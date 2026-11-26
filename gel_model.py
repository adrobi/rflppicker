"""Virtual-gel reference settings and expected fragment mass for viewing."""

from math import isfinite

DEFAULT_AGAROSE_PCT = 2.0
MIN_AGAROSE_PCT = 0.5
MAX_AGAROSE_PCT = 4.0
DEFAULT_GEL_MARKER = '50-bp'
DEFAULT_DIGEST_PCT = 100.0

# Generic reference strips with size labels in base pairs, independent of vendor.
GEL_MARKERS = {
    '50-bp': (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1500, 2000),
    '100-bp': (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000),
    '1-kb': (500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000),
}


def digest_probability(digest_pct):
    if not isfinite(digest_pct) or not 0 <= digest_pct <= 100:
        raise ValueError('gel_digest_pct: expected a percentage from 0 to 100')
    return digest_pct / 100


def digest_mass(ordered_fragments, probability):
    """Yield fragment length and expected mass in units of bp per product copy.

    An interval survives if its internal cuts fail and its two flanking cuts
    (where present) succeed. Enumerating contiguous intervals also covers
    coincident sizes, which the lane sums. No subset enumeration is needed.
    """
    if probability == 1:
        for size in ordered_fragments:
            yield size, size
        return
    if probability == 0:
        if ordered_fragments:
            size = sum(ordered_fragments)
            yield size, size
        return
    last = len(ordered_fragments) - 1
    for left in range(last + 1):
        length = 0
        chance = probability if left > 0 else 1.0
        for right in range(left, last + 1):
            length += ordered_fragments[right]
            mass = length * chance * (probability if right < last else 1.0)
            if mass:
                yield length, mass
            chance *= 1 - probability
