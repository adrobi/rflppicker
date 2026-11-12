"""Download a reference genome from NCBI by its assembly accession. No Qt.

The assembly folder of the NCBI FTP site holds the genome (``*_genomic.fna.gz``,
soft-masked), its MD5 and the assembly report with the chromosome names:

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/772/045/GCF_016772045.1_ARS-UI_Ramb_v2.0/

The result is ``genoms/<assembly name>.fa`` with its index and
``names/<assembly name>.txt`` (chromosome -> contig) next to the program.
"""
import gzip
import hashlib
import os
import re
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Callable, Iterator, List, Optional, Tuple

from i18n import N_, t

SITE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
CHUNK = 1 << 20
# Unpacked DNA takes about 3.3 times its archive; the archive is deleted only afterwards.
UNPACKED_RATIO = 4.5

# (title, species, accession, assembly name, size of the download in GB). The assemblies are
# fixed, not "the current reference": the coordinates of an SNP list belong to one assembly.
# Checked against the NCBI FTP site on 2026-11-13.
PRESETS = (
    (N_('Овца'), 'Ovis aries', 'GCF_016772045.2', 'ARS-UI_Ramb_v3.0', 0.83),
    (N_('Овца'), 'Ovis aries', 'GCF_016772045.1', 'ARS-UI_Ramb_v2.0', 0.83),
    (N_('Корова'), 'Bos taurus', 'GCF_002263795.3', 'ARS-UCD2.0', 0.85),
    (N_('Коза'), 'Capra hircus', 'GCF_001704415.2', 'ARS1.2', 0.89),
    (N_('Свинья'), 'Sus scrofa', 'GCF_000003025.6', 'Sscrofa11.1', 0.78),
    (N_('Курица'), 'Gallus gallus', 'GCF_016699485.2', 'bGalGal1.mat.broiler.GRCg7b', 0.33),
    (N_('Рис'), 'Oryza sativa', 'GCF_034140825.1', 'ASM3414082v1', 0.12),
    (N_('Пшеница'), 'Triticum aestivum', 'GCF_018294505.1', 'IWGSC_CS_RefSeq_v2.1', 4.34),
    (N_('Ячмень'), 'Hordeum vulgare', 'GCF_904849725.1', 'MorexV3_pseudomolecules_assembly', 1.27),
    (N_('Сорго'), 'Sorghum bicolor', 'GCF_000003195.3', 'Sorghum_bicolor_NCBIv3', 0.21),
)


class GenomeDownloadError(Exception):
    """The genome could not be downloaded; the text is shown to the user."""


class DownloadCancelled(Exception):
    """The caller cancelled the download; the unfinished archive is kept for a resume."""


@dataclass(frozen=True)
class Assembly:
    accession: str   # GCF_016772045.1
    name: str        # ARS-UI_Ramb_v2.0
    stem: str        # GCF_016772045.1_ARS-UI_Ramb_v2.0: the folder and the prefix of its files
    url: str         # the folder, with the closing slash


class HttpTransport:
    """The three kinds of request the download needs; replaced by a fake in the tests."""

    timeout = 60

    def _open(self, url: str, headers: Optional[dict] = None, method: str = 'GET'):
        # A generic client name: nothing about the user or the project leaves the computer.
        request = urllib.request.Request(url, headers={'User-Agent': 'RFLP-Picker', **(headers or {})},
                                         method=method)
        try:
            return urllib.request.urlopen(request, timeout=self.timeout)
        except urllib.error.HTTPError as error:
            raise GenomeDownloadError(t('NCBI ответил ошибкой {code} ({url})', code=error.code, url=url)) from error
        except (urllib.error.URLError, OSError) as error:
            raise GenomeDownloadError(t('Нет связи с NCBI: {error}', error=getattr(error, 'reason', error))) from error

    def text(self, url: str) -> str:
        with self._open(url) as response:
            return response.read().decode('utf-8', errors='replace')

    def size(self, url: str) -> int:
        with self._open(url, method='HEAD') as response:
            return int(response.headers.get('Content-Length') or 0)

    def stream(self, url: str, start: int) -> Iterator[bytes]:
        with self._open(url, {'Range': f'bytes={start}-'} if start else None) as response:
            if start and response.status != 206:
                raise GenomeDownloadError(t('NCBI не поддержал докачку; удалите незавершённый файл и повторите'))
            while True:
                try:
                    chunk = response.read(CHUNK)
                except OSError as error:
                    raise GenomeDownloadError(t('Соединение с NCBI прервано: {error}', error=error)) from error
                if not chunk:
                    return
                yield chunk


def accession_folder(accession: str) -> str:
    """The FTP folder that holds every version of an assembly."""
    match = re.fullmatch(r'(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+', accession.strip().upper())
    if match is None:
        raise GenomeDownloadError(t('Номер сборки должен иметь вид GCF_016772045.1, а не «{value}»', value=accession.strip()))
    return f'{SITE}{"/".join(match.groups())}/'


def find_assembly(accession: str, transport) -> Assembly:
    accession = accession.strip().upper()
    folder = accession_folder(accession)
    # The folder is named with the accession and the assembly name: GCF_..._ARS-UI_Ramb_v2.0/
    match = re.search(rf'href="({re.escape(accession)}_([^"/]+))/"', transport.text(folder))
    if match is None:
        raise GenomeDownloadError(t('Сборка {accession} не найдена на NCBI', accession=accession))
    stem, name = match.groups()
    return Assembly(accession, name, stem, f'{folder}{stem}/')


def chromosome_names(report: str, accession: str) -> List[Tuple[str, str]]:
    """(chromosome, contig) of the assembled molecules from an NCBI assembly report."""
    column = 6 if accession.upper().startswith('GCF') else 4     # RefSeq-Accn or GenBank-Accn
    names = []
    for line in report.splitlines():
        fields = line.split('\t')
        if line.startswith('#') or len(fields) < 7 or fields[1] != 'assembled-molecule':
            continue
        chromosome, contig = fields[2].strip(), fields[column].strip()
        if chromosome and contig and 'na' not in (chromosome, contig):
            names.append((chromosome, contig))
    return names


def _expected_md5(listing: str, file_name: str) -> str:
    for line in listing.splitlines():
        digest, _, path = line.strip().partition(' ')
        if path.strip().lstrip('./') == file_name:
            return digest.lower()
    raise GenomeDownloadError(t('На NCBI нет контрольной суммы файла {name}', name=file_name))


def _md5(path: str, cancelled: Callable[[], None]) -> str:
    digest = hashlib.md5()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(CHUNK), b''):
            cancelled()
            digest.update(chunk)
    return digest.hexdigest()


def download_genome(accession: str, base_dir: str, progress_cb=None, status_cb=None, cancel_cb=None,
                    transport=None, free_space=None) -> Tuple[str, str]:
    """Download, check, unpack and index an assembly; return (FASTA path, names path).

    ``progress_cb(done, total)`` counts the bytes of the archive, ``status_cb(text)``
    names the stage, ``cancel_cb()`` returning true raises DownloadCancelled. An
    unfinished archive is kept and continued by the next call.
    """
    transport = transport or HttpTransport()
    free_space = free_space or (lambda path: shutil.disk_usage(path).free)

    def cancelled() -> None:
        if cancel_cb and cancel_cb():
            raise DownloadCancelled()

    def status(text: str) -> None:
        if status_cb:
            status_cb(text)

    assembly = find_assembly(accession, transport)
    genomes, names_dir = os.path.join(base_dir, 'genoms'), os.path.join(base_dir, 'names')
    fasta = os.path.join(genomes, f'{assembly.name}.fa')
    names_path = os.path.join(names_dir, f'{assembly.name}.txt')
    if os.path.isfile(fasta):
        return fasta, names_path
    os.makedirs(genomes, exist_ok=True)

    archive_name = f'{assembly.stem}_genomic.fna.gz'
    archive_url = assembly.url + archive_name
    expected = _expected_md5(transport.text(assembly.url + 'md5checksums.txt'), archive_name)
    names = chromosome_names(transport.text(f'{assembly.url}{assembly.stem}_assembly_report.txt'), accession)
    total = transport.size(archive_url)
    part = fasta + '.gz.part'
    done = os.path.getsize(part) if os.path.isfile(part) else 0
    if done > total:          # not a part of this archive
        os.remove(part)
        done = 0
    needed = (total - done) + int(total * UNPACKED_RATIO)
    if free_space(genomes) < needed:
        raise GenomeDownloadError(t('Недостаточно места на диске: нужно около {size} ГБ в папке {folder}',
                                    size=f'{needed / 1e9:.1f}', folder=genomes))

    status(t('Загрузка с NCBI…'))
    if done < total:
        with open(part, 'ab') as output:
            for chunk in transport.stream(archive_url, done):
                cancelled()
                output.write(chunk)
                done += len(chunk)
                if progress_cb:
                    progress_cb(done, total)
    if progress_cb:
        progress_cb(done, total)

    status(t('Проверка контрольной суммы…'))
    if _md5(part, cancelled) != expected:
        os.remove(part)
        raise GenomeDownloadError(t('Контрольная сумма не совпала с NCBI: загрузка повреждена, повторите её'))

    status(t('Распаковка…'))
    unpacked = fasta + '.part'
    try:
        with gzip.open(part, 'rb') as source, open(unpacked, 'wb') as output:
            for chunk in iter(lambda: source.read(CHUNK), b''):
                cancelled()
                output.write(chunk)
    except BaseException:
        if os.path.exists(unpacked):
            os.remove(unpacked)
        raise
    os.replace(unpacked, fasta)
    os.remove(part)

    if names and not os.path.exists(names_path):      # a mapping edited by hand stays
        os.makedirs(names_dir, exist_ok=True)
        with open(names_path, 'w', encoding='utf-8', newline='\n') as output:
            output.writelines(f'{chromosome}\t{contig}\n' for chromosome, contig in names)

    status(t('Индексация…'))
    from pyfaidx import Fasta
    Fasta(fasta).close()
    return fasta, names_path
