"""Conservative, intensity-free scoring of diploid PCR-RFLP assays.

This is a deterministic screening model, not a calibrated gel simulator.
See docs/assay_design.md for equations, assumptions and ranking rules.
"""
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement
from math import isfinite, log2

from gel_model import DEFAULT_AGAROSE_PCT, MAX_AGAROSE_PCT, MIN_AGAROSE_PCT

MODEL_VERSION = "band-presence-maximin-v1"
AGAROSE_MODEL_VERSION = "band-presence-maximin-agarose-v2"


@dataclass(frozen=True)
class GelModel:
    min_visible: int = 80
    max_visible: int = 800
    resolution_bp: float = 10.0
    resolution_pct: float = 3.0
    stress_factor: float = 1.5
    version: str = MODEL_VERSION
    agarose_pct: float = DEFAULT_AGAROSE_PCT

    def __post_init__(self):
        values = (self.min_visible, self.max_visible, self.resolution_bp,
                  self.resolution_pct, self.stress_factor)
        if not all(isfinite(v) for v in values):
            raise ValueError("Параметры модели полос должны быть конечными числами")
        if not 0 < self.min_visible <= self.max_visible:
            raise ValueError("Некорректный диапазон видимых фрагментов")
        if self.resolution_bp <= 0 or not 0 <= self.resolution_pct <= 100:
            raise ValueError("Некорректное разрешение геля")
        if not 1 <= self.stress_factor <= 5:
            raise ValueError("Коэффициент ухудшения разрешения должен быть от 1 до 5")

        if self.version not in (MODEL_VERSION, AGAROSE_MODEL_VERSION):
            raise ValueError("Неизвестная версия модели геля")
        if not isfinite(self.agarose_pct) or not MIN_AGAROSE_PCT <= self.agarose_pct <= MAX_AGAROSE_PCT:
            raise ValueError("Концентрация агарозы должна быть от 0.5 до 4%")

    def distance(self, a, b, scale=1.0):
        length = max(a, b)
        threshold = max(self.resolution_bp, self.resolution_pct * length / 100.0)
        if self.version == AGAROSE_MODEL_VERSION and self.agarose_pct != DEFAULT_AGAROSE_PCT:
            # Engineering heuristic. Reference thresholds at 2% are unchanged.
            # Center the concentration adjustment at 600/c bp.
            reference_penalty = 1 + 0.5 * abs(log2(length / (600 / DEFAULT_AGAROSE_PCT)))
            penalty = 1 + 0.5 * abs(log2(length / (600 / self.agarose_pct)))
            threshold *= penalty / reference_penalty
        return abs(a - b) / (scale * threshold)


def band_clusters(lengths, model, scale):
    """Group unresolved consecutive sizes into bins shared across lanes."""
    clusters = []
    for size in sorted(set(lengths)):
        if clusters and model.distance(clusters[-1][-1], size, scale) <= 1:
            clusters[-1].append(size)
        else:
            clusters.append([size])
    return clusters


def _band_text(clusters, ids):
    return "; ".join(str(clusters[i][0]) if len(clusters[i]) == 1
                     else f"{clusters[i][0]}–{clusters[i][-1]}" for i in sorted(ids))


def _scenario(ref, alt, model, scale):
    clusters = band_clusters((*ref, *alt), model, scale)
    rset, aset = set(ref), set(alt)
    rids = {i for i, cluster in enumerate(clusters) if rset.intersection(cluster)}
    aids = {i for i, cluster in enumerate(clusters) if aset.intersection(cluster)}

    def directed_margin(source_ids, target_ids, target_sizes):
        # Missing entire lanes are not evidence of a reliable genotype call:
        # PCR failure cannot be distinguished from an invisible digest here.
        if not source_ids or not target_ids:
            return 0.0
        return max((min(model.distance(a, b, scale)
                        for a in clusters[i] for b in target_sizes)
                    for i in source_ids - target_ids), default=0.0)

    ref_unique = directed_margin(rids, aids, alt)
    alt_unique = directed_margin(aids, rids, ref)
    # RR vs RA requires an ALT-only band; RA vs AA requires a REF-only band.
    distances = {"REF/REF–REF/ALT": alt_unique,
                 "REF/ALT–ALT/ALT": ref_unique,
                 "REF/REF–ALT/ALT": max(ref_unique, alt_unique)}
    weakest = min(distances, key=distances.get)

    def lane(ids):
        return _band_text(clusters, ids)

    return {"margin": distances[weakest], "weakest": weakest,
            "bands_rr": lane(rids), "bands_ra": lane(rids | aids),
            "bands_aa": lane(aids), "band_count": len(rids | aids)}


def score_genotypes(frags_ref, frags_alt, model):
    """Return nominal and worst-resolution scores; score > 1 is resolvable.

    Equal-size copies have no extra weight. Short/large fragments are hidden.
    Consecutive unresolved lengths merge transitively (conservative bins).
    """
    if any(not isfinite(x) or x <= 0 for x in (*frags_ref, *frags_alt)):
        raise ValueError("Длины фрагментов должны быть положительными")
    ref = [x for x in frags_ref if model.min_visible <= x <= model.max_visible]
    alt = [x for x in frags_alt if model.min_visible <= x <= model.max_visible]
    nominal = _scenario(ref, alt, model, 1.0)
    stressed = _scenario(ref, alt, model, model.stress_factor)
    # Bins can only merge as resolution worsens, and distances decrease;
    # the upper endpoint is therefore the worst case over the whole interval.
    quality = ("robust" if stressed["margin"] > 1 else
               "nominal_only" if nominal["margin"] > 1 else "ambiguous")
    return {
        "genotype_quality": quality,
        "genotype_margin": nominal["margin"],
        "worst_margin": stressed["margin"],
        "weakest_genotypes": stressed["weakest"],
        "bands_ref_ref": nominal["bands_rr"],
        "bands_ref_alt": nominal["bands_ra"],
        "bands_alt_alt": nominal["bands_aa"],
        "worst_bands_ref_ref": stressed["bands_rr"],
        "worst_bands_ref_alt": stressed["bands_ra"],
        "worst_bands_alt_alt": stressed["bands_aa"],
        "visible_band_count": nominal["band_count"],
        "hidden_ref_count": len(frags_ref) - len(ref),
        "hidden_alt_count": len(frags_alt) - len(alt),
    }


def diploid_genotypes(count):
    names = ['REF', *(f'ALT{i}' for i in range(1, count))]
    return [(f'{names[i]}/{names[j]}', i, j) for i, j in combinations_with_replacement(range(count), 2)]


def _locus_scenario(lanes, model, scale, check_cancel):
    clusters = band_clusters((size for _, lengths in lanes for size in lengths), model, scale)
    cluster_ids = {size: i for i, cluster in enumerate(clusters) for size in cluster}
    ids = {label: {cluster_ids[size] for size in lengths} for label, lengths in lanes}

    def directed(source, target, target_lengths):
        if not source or not target:
            return 0.0
        return max((min(model.distance(a, b, scale) for a in clusters[i] for b in target_lengths)
                    for i in source - target), default=0.0)

    weakest, margin = '', float('inf')
    for (a, a_lengths), (b, b_lengths) in combinations(lanes, 2):
        check_cancel()
        value = max(directed(ids[a], ids[b], b_lengths), directed(ids[b], ids[a], a_lengths))
        if value < margin:
            weakest, margin = f'{a}–{b}', value
    return margin, weakest, {label: _band_text(clusters, lane) for label, lane in ids.items()}, len(clusters)


def score_locus(fragments, model, check_cancel=lambda: None):
    """Score all diploid genotype pairs using clusters shared across the locus."""
    visible = [{size for size in allele if model.min_visible <= size <= model.max_visible}
               for allele in fragments]
    lanes = [(label, visible[i] | visible[j]) for label, i, j in diploid_genotypes(len(fragments))]
    nominal, _, bands, count = _locus_scenario(lanes, model, 1.0, check_cancel)
    stressed, weakest, worst_bands, _ = _locus_scenario(lanes, model, model.stress_factor, check_cancel)
    quality = 'robust' if stressed > 1 else 'nominal_only' if nominal > 1 else 'ambiguous'
    return dict(genotype_quality=quality, genotype_margin=nominal, worst_margin=stressed,
                weakest_genotypes=weakest, genotype_bands=bands, worst_genotype_bands=worst_bands,
                visible_band_count=count)


def rank_key(record):
    """Lexicographic maximin ranking; no arbitrary weighted sum."""
    penalty = record.get("primer_pair_penalty")
    return (-record["worst_margin"], -record["genotype_margin"],
            penalty if penalty is not None else float("inf"),
            record["visible_band_count"], record["enzyme"],
            record.get("primer_pair_index") or 0)
