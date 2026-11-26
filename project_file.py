"""Versioned, portable snapshots of completed RFLP Picker runs."""

import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath

from Bio import __version__ as biopython_version
from pyfaidx import Fasta

import rflp_core as core
from app_settings import _same_type
from assay_scoring import GelModel
from enzyme_catalog import CATALOG_VERSION
from gel_model import DEFAULT_DIGEST_PCT, DEFAULT_GEL_MARKER, GEL_MARKERS, digest_probability
from i18n import t
from reporting import _atomic_destination
from snp_summary import InputRecord, SnpSummary, build_summaries
from visualization import payload_for_scheme

FORMAT = 'rflp-picker-project'
VERSION = 1
PRODUCT_FIELDS = ('enzyme_product', 'enzyme_product_2')
_ACCESSION = re.compile(r'(?<![A-Z0-9])GC[AF]_\d+\.\d+(?!\d)', re.IGNORECASE)


def run_metadata(params: dict, report: dict) -> dict:
    genome_path = params.get('fasta', '')
    genome_name = PurePosixPath(genome_path.replace('\\', '/')).name
    accession = _ACCESSION.search(genome_name)
    return {
        'app_version': core.APP_VERSION,
        'scoring_model': report.get('scoring_model', core.MODEL_VERSION),
        'enzyme_catalog_version': CATALOG_VERSION,
        'restriction_dictionary_version': biopython_version,
        'assembly_accession': params.get('assembly_accession') or (accession[0].upper() if accession else ''),
        'run_timestamp': datetime.now().isoformat(timespec='seconds'),
        'genome_name': genome_name,
        'genome_path': genome_path,
    }


def capture_payloads(summaries: list[SnpSummary], params: dict) -> dict:
    """Store the drawing inputs of every selectable row, opening FASTA once."""
    payloads = {}
    fasta = None
    path = params.get('fasta', '')
    try:
        for summary in summaries:
            for scheme in summary.schemes:
                key = str(scheme['_index'])
                if key in payloads:
                    continue
                if (scheme.get('design') == 'allele_specific' or
                        not all(scheme.get(field) for field in ('mapped_id', 'amplicon_start', 'amplicon_end'))):
                    payloads[key] = None
                    continue
                if fasta is None:
                    fasta = Fasta(path, as_raw=True, sequence_always_upper=True)
                payload, error = payload_for_scheme(scheme, path, params, fasta)
                if payload is None:
                    raise ValueError(t('Не удалось сохранить схему {variant}: {error}',
                                       variant=scheme['variant'], error=error))
                payloads[key] = payload
    finally:
        if fasta is not None:
            fasta.close()
    return payloads


def cached_payload(payloads: dict, scheme: dict) -> tuple[dict | None, str]:
    payload = payloads.get(str(scheme.get('_index')))
    return payload, '' if payload else t('В проекте нет данных для визуализации этой схемы.')


def save_project(path, *, header, rows, report, params, metadata, has_ids, summaries, payloads) -> None:
    products = {}
    for summary in summaries:
        for scheme in summary.schemes:
            if any(field in scheme for field in PRODUCT_FIELDS):
                products[str(scheme['_index'])] = [scheme.get(field) for field in PRODUCT_FIELDS]
    document = {
        'format': FORMAT, 'version': VERSION, 'metadata': metadata, 'parameters': params,
        'results': {'header': header, 'rows': rows, 'report': report},
        'has_ids': has_ids, 'records': [asdict(summary.record) for summary in summaries],
        'choices': [summary.chosen for summary in summaries], 'products': products, 'payloads': payloads,
    }
    with _atomic_destination(path) as temporary:
        with temporary.open('w', encoding='utf-8') as stream:
            json.dump(document, stream, ensure_ascii=False, allow_nan=False, separators=(',', ':'))


@dataclass
class OpenProject:
    header: list
    rows: list
    report: dict
    params: dict
    metadata: dict
    has_ids: bool
    summaries: list[SnpSummary]
    payloads: dict


def _check(condition: bool) -> None:
    if not condition:
        raise ValueError(t('Повреждены данные проекта.'))


def _check_parameters(params: dict) -> None:
    defaults = dict(core.DEFAULT_PARAMS, fasta='', input_file='')
    for key, default in defaults.items():
        value = params.get(key, default)
        _check(_same_type(value, default))
        if isinstance(default, float):
            _check(math.isfinite(value))
    # These two values go directly to QSpinBox, whose input is a C++ int.
    _check(all(0 < params.get(key, defaults[key]) < 2**31 for key in ('prod_min', 'prod_max')))
    _check(not core.validate_params({**defaults, **params}, check_files=False))


def _check_payload(payload) -> None:
    if payload is None:
        return
    _check(isinstance(payload, dict) and isinstance(payload.get('gel'), dict))
    _check(all(isinstance(payload.get(key), str) for key in ('mapped_id', 'amplicon_seq', 'primer_left', 'primer_right')))
    _check(all(type(payload.get(key)) is int for key in ('genome_start', 'genome_end')))
    _check(1 <= payload['genome_start'] <= payload['genome_end'])
    _check(all(key in payload and (payload[key] is None or
               (type(payload[key]) in (int, float) and math.isfinite(payload[key])))
               for key in ('snp_offset', 'tm_left', 'tm_right')))
    _check(all(isinstance(payload.get(key, ''), str) for key in
               ('variant_str', 'enzyme_name', 'pattern', 'frags_ref', 'frags_alt')))
    enzymes, alleles = payload.get('enzyme_names', ()), payload.get('locus_products', ())
    _check(isinstance(enzymes, (list, tuple)) and all(isinstance(item, str) for item in enzymes))
    _check(isinstance(alleles, (list, tuple)) and (not alleles or len(alleles) >= 2))
    for allele in alleles:
        _check(isinstance(allele, dict) and all(isinstance(allele.get(key), str) and allele[key]
                                                for key in ('allele', 'sequence')))
        _check(type(allele.get('length')) is int and allele['length'] > 0)
        fragments = allele.get('fragments')
        _check(isinstance(fragments, list) and all(type(value) is int and value > 0 for value in fragments))
        _check(sum(fragments) == allele['length'])
    gel = dict(payload['gel'])
    _check(gel.pop('marker', DEFAULT_GEL_MARKER) in GEL_MARKERS)
    digest_probability(gel.pop('digest_pct', DEFAULT_DIGEST_PCT))
    if 'model_version' in gel:
        gel['version'] = gel.pop('model_version')
    GelModel(**gel)


def _restore(document: dict) -> OpenProject:
    try:
        metadata, params = document['metadata'], document['parameters']
        header, rows, report = (document['results'][name] for name in ('header', 'rows', 'report'))
        has_ids, records = document['has_ids'], document['records']
        choices, products, payloads = (document[name] for name in ('choices', 'products', 'payloads'))
        _check(all(isinstance(value, dict) for value in (metadata, params, report, products, payloads)))
        _check_parameters(params)
        _check(isinstance(header, list) and all(isinstance(name, str) for name in header))
        _check(isinstance(rows, list) and all(isinstance(row, list) and len(row) == len(header) for row in rows))
        _check(isinstance(has_ids, bool) and isinstance(records, list) and isinstance(choices, list))
        _check(len(records) == len(choices))
        row_lines = report.get('result_source_lines')
        if row_lines is not None:
            _check(isinstance(row_lines, list) and len(row_lines) == len(rows) and
                   all(type(line) is int and line > 0 for line in row_lines))
        for payload in payloads.values():
            _check_payload(payload)
        input_records = [InputRecord(**record) for record in records]
        _check(all(isinstance(record.line, int) and (record.pos is None or isinstance(record.pos, int))
                   and all(isinstance(value, str) for value in (record.snp_id, record.chrom, record.ref,
                                                               record.alt, record.raw)) for record in input_records))
        if row_lines is not None:
            _check(set(row_lines) <= {record.line for record in input_records})
        summaries = build_summaries(input_records, header, rows, report.get('rejected', []), report=report)
        for summary, index in zip(summaries, choices):
            _check(type(index) is int and 0 <= index < max(1, len(summary.schemes)))
            if summary.schemes:
                summary.choose(index)
            for scheme in summary.schemes:
                key = str(scheme['_index'])
                _check(key in payloads)
                if key in products:
                    selected = products[key]
                    _check(isinstance(selected, list) and len(selected) == len(PRODUCT_FIELDS))
                    _check(all(value is None or isinstance(value, str) for value in selected))
                    scheme.update(zip(PRODUCT_FIELDS, selected))
    except (AttributeError, IndexError, KeyError, OverflowError, TypeError) as exc:
        raise ValueError(t('Повреждены данные проекта.')) from exc
    return OpenProject(header, rows, report, params, metadata, has_ids, summaries, payloads)


def _reject_json_constant(_value):
    _check(False)


def open_project(path) -> OpenProject:
    with Path(path).open(encoding='utf-8') as stream:
        document = json.load(stream, parse_constant=_reject_json_constant)
    _check(isinstance(document, dict))
    if document.get('format') != FORMAT or type(document.get('version')) is not int or document['version'] != VERSION:
        raise ValueError(t('Неподдерживаемый формат проекта.'))
    return _restore(document)
