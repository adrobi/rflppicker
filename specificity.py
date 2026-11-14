"""Where else in the genome do the primers anneal, and what would a pair amplify there? No Qt.

One pass over the whole FASTA looks for the 3' ends of all primers at once; no
index is kept and no other program is needed. The model is deliberately simple:

* a primer anneals where its last SEED bases match the genome exactly and the
  rest of it has at most ``max_mismatches`` mismatches (an N of the genome is a
  mismatch), on either strand;
* a pair gives a product where a site on the top strand and a site on the bottom
  strand of one contig face each other within ``max_size`` bp — the two primers
  in any combination, also one primer on both sides.

It does not judge how well a mismatched primer really anneals.
"""
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

import numpy as np
from pyfaidx import Fasta

SEED = 12                 # bases of the 3' end that must match exactly
BLOCK = 1 << 23           # genome bases read at once
MAX_SITES = 1000          # a primer with more sites is a repeat: its products are not counted
_COMPLEMENT = str.maketrans('ACGT', 'TGCA')
_CODE = np.full(256, 255, dtype=np.uint8)          # A, C, G, T -> 0..3 in either case; anything else 255
for _value, _base in enumerate('ACGT'):
    _CODE[ord(_base)] = _CODE[ord(_base.lower())] = _value


class ScanCancelled(Exception):
    """The caller cancelled the pass over the genome."""


class Site(NamedTuple):
    contig: str
    start: int            # 0-based start of the primer footprint on the top strand
    strand: str           # '+': the primer reads along the top strand; '-': along the bottom strand
    mismatches: int


@dataclass
class PrimerSites:
    sites: List[Site] = field(default_factory=list)
    saturated: bool = False    # more than the limit: a repeat, the list is cut


class Product(NamedTuple):
    contig: str
    start: int            # 0-based, top strand
    end: int              # exclusive
    primers: str          # 'F+R', 'R+F', 'F+F' or 'R+R': the forward-strand primer first

    @property
    def size(self) -> int:
        return self.end - self.start


def _codes(text: str) -> np.ndarray:
    return _CODE[np.frombuffer(text.encode('ascii'), dtype=np.uint8)]


def _seed_value(codes: np.ndarray) -> int:
    value = 0
    for code in codes:
        value = (value << 2) | int(code)
    return value


def _seed_values(codes: np.ndarray) -> np.ndarray:
    """The value of the SEED-mer starting at every position (meaningless where it spans a non-ACGT base)."""
    count = len(codes) - SEED + 1
    two = (codes[:-1] << 2) | (codes[1:] & 3)
    four = (two[:-2].astype(np.uint32) << 4) | two[2:]
    return ((four[:count] << 16) | (four[4:count + 4] << 8) | four[8:count + 8]) & 0xFFFFFF


def _scan_block(codes: np.ndarray, first: int, last: int, table: np.ndarray, patterns: dict, max_mismatches: int):
    """Sites whose seed starts in codes[first:last]: (primer index, footprint start in codes, strand, mismatches)."""
    found = []
    if len(codes) < SEED:
        return found
    values = _seed_values(codes)
    positions = np.flatnonzero(table[values[first:last]]) + first
    if not len(positions):
        return found
    hit_values = values[positions]
    order = np.argsort(hit_values, kind='stable')
    positions, hit_values = positions[order], hit_values[order]
    borders = np.flatnonzero(np.diff(hit_values)) + 1
    span = np.arange(SEED)
    for group in np.split(np.arange(len(positions)), borders):
        seeds = positions[group]
        # A seed over an N has a meaningless value: the seed itself must be real bases.
        seeds = seeds[(codes[seeds[:, None] + span] != 255).all(axis=1)]
        if not len(seeds):
            continue
        for index, strand, pattern in patterns[int(hit_values[group[0]])]:
            # '+': the seed is the end of the footprint; '-': the reverse complement starts with it.
            starts = seeds - (len(pattern) - SEED) if strand == '+' else seeds
            starts = starts[(starts >= 0) & (starts + len(pattern) <= len(codes))]
            if not len(starts):
                continue
            mismatches = (codes[starts[:, None] + np.arange(len(pattern))] != pattern).sum(axis=1)
            keep = mismatches <= max_mismatches
            found.extend((index, int(start), strand, int(count)) for start, count in zip(starts[keep], mismatches[keep]))
    return found


def find_sites(fasta_path: str, primers: Iterable[str], max_mismatches: int = 3, threads: Optional[int] = None,
               block: int = BLOCK, max_sites: int = MAX_SITES, progress_cb=None, cancel_cb=None) -> Dict[str, PrimerSites]:
    """The annealing sites of every primer in the genome, by primer in upper case.

    Primers shorter than the seed or with letters other than A, C, G, T are left
    out of the result. ``progress_cb(done, total)`` counts blocks of the genome.
    """
    names = []
    for primer in primers:
        primer = primer.upper()
        if len(primer) >= SEED and not set(primer) - set('ACGT') and primer not in names:
            names.append(primer)
    result = {primer: PrimerSites() for primer in names}
    if not names:
        return result

    table = np.zeros(4 ** SEED, dtype=bool)
    patterns: Dict[int, list] = {}
    longest = max(map(len, names))
    for index, primer in enumerate(names):
        forward = _codes(primer)
        reverse = _codes(primer.translate(_COMPLEMENT)[::-1])
        for strand, pattern, seed in (('+', forward, forward[-SEED:]), ('-', reverse, reverse[:SEED])):
            value = _seed_value(seed)
            table[value] = True
            patterns.setdefault(value, []).append((index, strand, pattern))

    with Fasta(fasta_path, as_raw=True) as index_only:
        lengths = [(name, len(index_only[name])) for name in index_only.keys()]
    jobs = [(name, start, min(length, start + block), length)
            for name, length in lengths for start in range(0, length, block)]
    local = threading.local()
    handles = []

    def work(job):
        if cancel_cb and cancel_cb():
            return None
        if not hasattr(local, 'fasta'):
            local.fasta = Fasta(fasta_path, as_raw=True)
            handles.append(local.fasta)
        name, start, end, length = job
        # The footprint of a site may reach beyond the block by the primer length.
        left, right = max(0, start - longest), min(length, end + longest)
        codes = _codes(str(local.fasta[name][left:right]))
        return name, left, _scan_block(codes, start - left, end - left, table, patterns, max_mismatches)

    workers = max(1, min(threads or 4, len(jobs)))
    try:
        with ThreadPoolExecutor(workers) as pool:
            for done, found in enumerate(pool.map(work, jobs), 1):
                if found is None or (cancel_cb and cancel_cb()):
                    raise ScanCancelled()
                name, left, sites = found
                for index, start, strand, mismatches in sorted(sites, key=lambda site: (site[1], site[2])):
                    entry = result[names[index]]
                    if len(entry.sites) >= max_sites:
                        entry.saturated = True
                    else:
                        entry.sites.append(Site(name, left + start, strand, mismatches))
                if progress_cb:
                    progress_cb(done, len(jobs))
    finally:
        for handle in handles:
            handle.close()
    return result


def pair_products(sites: Dict[str, PrimerSites], left: str, right: str, max_size: int,
                  intended: Optional[Tuple[str, int, int]] = None) -> List[Product]:
    """What the pair could amplify: (contig, start, end) of every product up to ``max_size`` bp.

    ``intended`` is the (contig, start, end) of the designed amplicon; it is left out.
    """
    forward, reverse = [], []        # (contig, coordinate, primer label)
    for label, primer in (('F', left.upper()), ('R', right.upper())):
        entry = sites.get(primer)
        if entry is None:
            continue
        for site in entry.sites:
            if site.strand == '+':
                forward.append((site.contig, site.start, label, len(primer)))
            else:
                reverse.append((site.contig, site.start + len(primer), label, len(primer)))
    reverse.sort()
    products = []
    for contig, start, label, length in sorted(forward):
        for other_contig, end, other_label, other_length in reverse:
            if other_contig != contig or end - start > max_size:
                continue
            # The two footprints must not overlap: there is nothing to extend into.
            if end - other_length < start + length:
                continue
            if intended is not None and (contig, start, end) == tuple(intended):
                continue
            products.append(Product(contig, start, end, f'{label}+{other_label}'))
    return products
