"""Resolve an rsID within a selected species and versioned reference assembly."""

import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from http.client import HTTPException
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from pyfaidx import Fasta, FastaIndexingError

from genome_download import GenomeDownloadError, HttpTransport, find_assembly
from variant_input import VCF_DNA, open_vcf


EVA_URL = 'https://www.ebi.ac.uk/eva/webservices/identifiers/v1/clustered-variants/{number}/submitted'
_ACCESSION = re.compile(r'GC[AF]_[0-9]{9}\.[1-9][0-9]*', re.IGNORECASE | re.ASCII)
_ASSEMBLY_LABEL = re.compile(r'GC[AF]_[0-9]{9}(?:\.[0-9]+)?', re.IGNORECASE | re.ASCII)
_TAXON = re.compile(r'(?:(?:NCBITaxon|taxid):\s*)?([1-9][0-9]*)', re.IGNORECASE | re.ASCII)
_ATTRIBUTE = re.compile(r'(?:^|,)\s*([\w-]+)\s*=\s*(?:"((?:\\.|[^"\\])*)"|([^,]*))')


@dataclass(frozen=True)
class VariantMatch:
    taxid: int
    assembly: str
    contig: str
    position: int
    ref: str
    alts: tuple[str, ...]
    source_contig: str
    species: str = ''


@dataclass(frozen=True)
class LookupResult:
    rsid: str
    status: str
    matches: tuple[VariantMatch, ...] = ()
    detail: str = ''

    @property
    def match(self):
        return self.matches[0] if self.status == 'accepted' else None


@dataclass(frozen=True)
class _Record:
    contig: str
    position: int
    ref: str
    alts: tuple[str, ...]
    assemblies: tuple[str, ...]
    taxa: tuple[str, ...]


@dataclass
class _AssemblyReport:
    name: str = ''
    paired: str = ''
    organism: str = ''
    taxid: str = ''
    aliases: dict = field(default_factory=lambda: defaultdict(set))
    shared_contigs: set = field(default_factory=set)


def _read_report(contents, assembly):
    metadata, rows = {}, []
    for line in contents.splitlines():
        if line.startswith('# '):
            key, separator, value = line[2:].partition(':')
            if separator:
                metadata[key.casefold()] = value.strip()
        elif line and not line.startswith('#'):
            rows.append(line.split('\t'))
    genbank = metadata.get('genbank assembly accession', '').upper()
    refseq = metadata.get('refseq assembly accession', '').upper()
    if assembly not in (genbank, refseq):
        raise ValueError(f'Assembly report does not include {assembly}')
    bank_selected = assembly == genbank
    report = _AssemblyReport(metadata.get('assembly name', ''), refseq if bank_selected else genbank,
                             metadata.get('organism name', ''), metadata.get('taxid', ''))
    for row in rows:
        if len(row) < 10:
            continue
        name, role, chrom, _, bank_contig, relation, refseq_contig, _, _, ucsc_name = row[:10]
        own, other = (bank_contig, refseq_contig) if bank_selected else (refseq_contig, bank_contig)
        if own in ('na', ''):
            continue
        spellings = {own, name, ucsc_name}
        if role == 'assembled-molecule':
            spellings.add(chrom)
        spellings -= {'na', ''}
        if relation == '=' and other not in ('na', ''):
            spellings.add(other)
            report.shared_contigs.update(spellings)
        for alias in spellings:
            report.aliases[alias].update(spellings)
    return report


def _sidecar(fasta):
    reference = Path(fasta)
    if reference.suffix.lower() in ('.gz', '.bgz', '.bgzf'):
        reference = reference.with_suffix('')
    stem = reference.stem.removesuffix('_genomic')
    report = reference.with_name(stem + '_assembly_report.txt')
    return report if report.is_file() else None


def _build_matches(value, contig, assembly, report):
    accessions = _ASSEMBLY_LABEL.findall(value.upper())
    if not accessions:
        return bool(report.name) and value.casefold() == report.name.casefold()
    paired_contig = contig in report.shared_contigs or contig.removeprefix('chr') in report.shared_contigs
    return all(accession == assembly or (accession == report.paired and paired_contig)
               for accession in accessions)


def _taxon_matches(value, taxid, species, report):
    match = _TAXON.fullmatch(value)
    if match:
        return match[1] == str(taxid)
    labels = {species.casefold()}
    if report.taxid == str(taxid):
        organism = report.organism.split(' (', 1)[0].casefold()
        labels.update((organism, ' '.join(organism.split()[:2])))
    return bool(value) and value.casefold() in labels


def _named_reference(value):
    # A bare build label can be compared with the report. An arbitrary FASTA
    # filename or URI may lack assembly metadata altogether.
    return (re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.+-]*', value) is not None
            and Path(value).suffix.lower() not in ('.fasta', '.fa', '.fna', '.gz', '.bgz', '.bgzf'))


def _local_records(path, rsid, check_cancel):
    builds, taxa = [], []
    contigs = defaultdict(list)
    header = False
    with open_vcf(path) as stream:
        for row_number, line in enumerate(stream, 1):
            check_cancel()
            line = line.rstrip('\r\n')
            if line.startswith('##'):
                field, _, value = line[2:].partition('=')
                field, value = field.casefold(), value.strip().strip('"')
                if field == 'assembly' and ('://' not in value or _ASSEMBLY_LABEL.search(value)):
                    builds.append(value)
                elif field in ('reference', 'reference_url') and (
                    _ASSEMBLY_LABEL.search(value) or (field == 'reference' and _named_reference(value))
                ):
                    builds.append(value)
                elif field in ('taxid', 'species', 'taxonomy'):
                    taxa.append(value)
                elif field == 'contig' and value.startswith('<') and value.endswith('>'):
                    attrs = {key.casefold(): quoted if quoted else bare.strip()
                             for key, quoted, bare in _ATTRIBUTE.findall(value[1:-1])}
                    if attrs.get('id'):
                        contigs[attrs['id']].append(attrs)
                continue
            if line.startswith('#CHROM\t'):
                header = True
            if not line or line.startswith('#'):
                continue
            if not header:
                raise ValueError('Missing VCF #CHROM header')
            fields = line.split('\t', 8)
            if len(fields) < 3 or rsid not in fields[2].lower().split(';'):
                continue
            if len(fields) < 8 or not re.fullmatch(r'[1-9][0-9]*', fields[1]):
                raise ValueError(f'Invalid VCF record at line {row_number}')
            contig, position, _, ref, alt = fields[:5]
            attributes = contigs.get(contig, ())
            contig_builds = [x['assembly'] for x in attributes if 'assembly' in x]
            contig_taxa = [x[key] for x in attributes for key in ('taxid', 'species', 'taxonomy') if key in x]
            yield _Record(contig, int(position), ref.upper(), tuple(alt.upper().split(',')),
                          tuple(builds + contig_builds), tuple(taxa + contig_taxa))
    if not header:
        raise ValueError('Missing VCF #CHROM header')


def _fetch_json(url):
    request = Request(url, headers={'User-Agent': 'RFLP-Picker', 'Accept': 'application/json'})
    try:
        with urlopen(request, timeout=15, context=HttpTransport._ssl_context()) as answer:
            return json.load(answer)
    except HTTPError as error:
        if error.code == 404:
            return []
        raise


def _online_records(payload):
    if not isinstance(payload, list):
        raise ValueError('Unexpected EVA response type')
    for entry in payload:
        try:
            record = entry['data']
            contig = record['contig']
            value = record['start']
            ref, alt = record['referenceAllele'].upper(), record['alternateAllele'].upper()
            assembly = record['referenceSequenceAccession']
            taxon = str(record['taxonomyAccession'])
            if not isinstance(contig, str) or not contig or not _ACCESSION.fullmatch(assembly):
                raise ValueError('Invalid EVA submitted coordinate')
            if not re.fullmatch(r'[1-9][0-9]*', str(value)):
                raise ValueError('Invalid EVA submitted coordinate')
        except (AttributeError, KeyError, TypeError) as error:
            raise ValueError('Unexpected EVA submitted record') from error
        yield _Record(contig, int(value), ref, (alt,), (assembly,), (taxon,))


def _online_report(assembly):
    transport = HttpTransport()
    transport.timeout = 15
    entry = find_assembly(assembly, transport)
    return transport.text(f'{entry.url}{entry.stem}_assembly_report.txt')


def _fasta_contigs(contig, fasta, names, report):
    configured = names.get(contig)
    if configured and configured in fasta:
        return [configured]
    if contig in fasta:
        return [contig]
    alternatives = report.aliases.get(contig, set()) | report.aliases.get(contig.removeprefix('chr'), set()) | {contig}
    candidates = set(alternatives)
    for alias in alternatives:
        candidates.add(alias.removeprefix('chr') if alias.startswith('chr') else 'chr' + alias)
        if alias in names:
            candidates.add(names[alias])
    return sorted(candidate for candidate in candidates if candidate in fasta)


def _check_matches(records, fasta, names, report, rsid, assembly, taxid, species):
    grouped = {}
    absent = set()
    for record in records:
        choices = _fasta_contigs(record.contig, fasta, names, report)
        if not choices:
            choices = [record.contig]
            absent.add(record.contig)
        for contig in choices:
            grouped.setdefault((contig, record.position, record.ref), []).append(record)
    matches = []
    for (contig, position, ref), members in grouped.items():
        alts = tuple(dict.fromkeys(alt for record in members for alt in record.alts
                                   if alt and alt != ref and set(alt) <= VCF_DNA))
        matches.append(VariantMatch(taxid, assembly, contig, position, ref, alts, members[0].contig, species))
    matches = tuple(sorted(matches, key=lambda match: (match.contig, match.position, match.ref)))
    if len(matches) > 1:
        return LookupResult(rsid, 'ambiguous', matches)
    match = matches[0]
    if match.contig in absent:
        return LookupResult(rsid, 'contig_not_found', matches)
    if not (match.ref and set(match.ref) <= VCF_DNA and match.alts):
        return LookupResult(rsid, 'unsupported_allele', matches)
    start = match.position - 1
    end = start + len(match.ref)
    reference = fasta[match.contig]
    observed = str(reference[start:end]).upper() if end <= len(reference) else ''
    if observed != match.ref:
        return LookupResult(rsid, 'ref_mismatch', matches, observed)
    return LookupResult(rsid, 'accepted', matches)


def lookup_rsid(rsid, *, fasta, taxid, assembly, local_vcf=None, names=None, assembly_report=None,
                species='', fetch_json=None, check_cancel=None) -> LookupResult:
    """Find one REF-verified placement, or return the corresponding nonacceptance status.

    Without reference metadata in a local VCF, the selected build and taxon are
    taken from the caller. Paired GenBank/RefSeq assembly coordinates require
    equality at the contig level in an NCBI assembly report. A local lookup
    uses only local files. ``fetch_json`` substitutes the network for tests.
    """
    rsid, assembly, species = rsid.strip().lower(), assembly.strip().upper(), species.strip()
    query = re.fullmatch(r'rs([1-9][0-9]*)', rsid, re.ASCII)
    taxon = _TAXON.fullmatch(str(taxid).strip())
    if not query or not _ACCESSION.fullmatch(assembly) or not taxon:
        raise ValueError('Expected an rsID, a TaxID and a versioned GCA or GCF assembly accession')
    taxid = int(taxon[1])
    cancel = check_cancel if check_cancel is not None else lambda: None
    cancel()
    try:
        if local_vcf:
            records = list(_local_records(local_vcf, rsid, cancel))
        else:
            records = list(_online_records((fetch_json or _fetch_json)(EVA_URL.format(number=query[1]))))
        cancel()
        if not records:
            return LookupResult(rsid, 'no_match')
        report_path = assembly_report or _sidecar(fasta)
        report = _AssemblyReport()
        if report_path:
            report = _read_report(Path(report_path).read_text(encoding='utf-8-sig'), assembly)
            if local_vcf and report.taxid and report.taxid != str(taxid):
                return LookupResult(rsid, 'species_mismatch', detail=report.organism or f'TaxID {report.taxid}')
        elif not (local_vcf or fetch_json) and assembly.startswith('GCF_') and any(
            record.taxa == (str(taxid),) and record.assemblies[0].upper().startswith('GCA_')
            for record in records
        ):
            report = _read_report(_online_report(assembly), assembly)
        taxon_records = [record for record in records
                         if all(_taxon_matches(taxon, taxid, species, report) for taxon in record.taxa)]
        if not taxon_records:
            return LookupResult(rsid, 'species_mismatch')
        assembly_records = [record for record in taxon_records if all(
            _build_matches(value, record.contig, assembly, report) for value in record.assemblies)]
        if not assembly_records:
            return LookupResult(rsid, 'assembly_mismatch')
        cancel()
        species = species or (report.organism if report.taxid == str(taxid) else '')
        with Fasta(str(fasta), as_raw=True) as reference:
            return _check_matches(assembly_records, reference, names or {}, report, rsid, assembly, taxid, species)
    except (OSError, ValueError, FastaIndexingError, GenomeDownloadError, HTTPException) as error:
        return LookupResult(rsid, 'source_error', detail=str(error))
