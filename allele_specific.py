"""Two-reaction allele-specific PCR with a shared flanking control pair."""

import json
from collections.abc import Mapping
from itertools import product
from math import ceil, floor


def interpret_genotype(*, ref_control, alt_control, ref_band, alt_band):
    """Return a diploid genotype from identified bands, or None when undetermined."""
    if not (ref_control and alt_control):
        return None
    return {(True, False): 'REF/REF', (True, True): 'REF/ALT', (False, True): 'ALT/ALT'}.get(
        (ref_band, alt_band))


def expected_bands(control_size, ref_size, alt_size):
    """Expected lengths in the REF and ALT tubes for each diploid genotype."""
    return {
        'REF/REF': {'REF': [ref_size, control_size], 'ALT': [control_size]},
        'REF/ALT': {'REF': [ref_size, control_size], 'ALT': [alt_size, control_size]},
        'ALT/ALT': {'REF': [control_size], 'ALT': [alt_size, control_size]},
    }


def _read_mapping(value):
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except ValueError:
            return None
    return value if isinstance(value, Mapping) else None


def reaction_table(result):
    """Return tube primers and expected bands from a result or exported row mapping."""
    side = result.get('aspcr_specific_side')
    control = _read_mapping(result.get('aspcr_control'))
    reactions = _read_mapping(result.get('aspcr_reactions'))
    genotypes = _read_mapping(result.get('aspcr_genotypes'))
    if side not in ('left', 'right') or not (control and reactions and genotypes):
        return []
    common = 'right' if side == 'left' else 'left'
    control_size = control['product_size']
    tubes = []
    for label in ('REF', 'ALT'):
        reaction = reactions[label]
        primers = []
        for role, primer_side, pair in (('allele_specific', side, reaction),
                                        ('common', common, control), ('control_only', side, control)):
            primers.append(dict(role=role, side=primer_side, sequence=pair[f'primers_{primer_side}'],
                                tm=pair[f'tm_{primer_side}']))
        tubes.append(dict(reaction=label, allele=reaction['allele'], primers=primers,
                          control_size=control_size, product_size=reaction['product_size'],
                          genotype_bands={genotype: list(lanes[label]) for genotype, lanes in genotypes.items()}))
    return tubes


def _nested_pairs(template, offset, side, control, *, pick_primers, parse_pairs,
                  min_size, max_size, tm_min, tm_max, excluded_regions, mask_repeats):
    common = 'right' if side == 'left' else 'left'
    target = offset + 1 if side == 'left' else offset - 1
    shared_primer = control[f'primers_{common}']
    pres = pick_primers(template, target, min_size, max_size, tm_min, tm_max,
                        num_return=5, lowercase_masking=mask_repeats, excluded_regions=excluded_regions,
                        **{f'fixed_{common}': shared_primer, f'force_{side}_end': offset})
    if not pres:
        return []
    try:
        pairs = list(parse_pairs(pres, template.upper(), target))
    except (KeyError, TypeError, ValueError):
        return []
    keep = []
    for pair in pairs:
        start, end = pair['start0'], pair['end0']
        left_len, right_len = pair['primer_left_len'], pair['primer_right_len']
        if side == 'left':
            nested = (start >= control['start0'] + control['primer_left_len']
                      and end == control['end0'] and start + left_len - 1 == offset)
        else:
            nested = (end <= control['end0'] - control['primer_right_len']
                      and start == control['start0'] and end - right_len == offset)
        if not nested or pair[f'primers_{common}'] != shared_primer:
            continue
        if pair[f'primers_{side}'] in (control['primers_left'], control['primers_right']):
            continue
        if not (min_size <= end - start <= max_size and tm_min <= pair[f'tm_{side}'] <= tm_max):
            continue
        if any(start < pos + length and pos < start + left_len
               or end - right_len < pos + length and pos < end for pos, length in excluded_regions):
            continue
        if mask_repeats and (template[start + left_len - 1].islower() or template[end - right_len].islower()):
            continue
        keep.append(dict(pair, product_size=end - start))
    return keep


def _score(control, ref, alt, model):
    control_size, ref_size, alt_size = [p['product_size'] for p in (control, ref, alt)]
    ref_margin = model.distance(ref_size, control_size)
    alt_margin = model.distance(alt_size, control_size)
    margin = min(ref_margin, alt_margin)
    worst_margin = min(model.distance(ref_size, control_size, model.stress_factor),
                       model.distance(alt_size, control_size, model.stress_factor))
    if margin <= 1:
        return None
    return dict(genotype_quality='robust' if worst_margin > 1 else 'nominal_only',
                genotype_margin=margin, worst_margin=worst_margin,
                weakest_genotypes='REF/ALT–ALT/ALT' if ref_margin <= alt_margin else 'REF/REF–REF/ALT',
                visible_band_count=4, hidden_ref_count=0, hidden_alt_count=0,
                aspcr_genotypes=expected_bands(control_size, ref_size, alt_size))


def design_assays(sequence, offset, ref, alt, *, control_pairs, pick_primers, parse_pairs,
                  prod_min, prod_max, tm_min, tm_max, model, excluded_regions=(),
                  mask_repeats=False, check_cancel=lambda: None):
    """Design a pair of diagnostic tubes for one SNP using Primer3 control pairs.

    Each tube contains the outer control pair and one inner primer that ends at
    the SNP. The inner primer shares its partner with the outer pair. Coordinates
    are window-relative, with an exclusive end. The requested PCR size range
    applies to the control product. Nested products must be visible and resolve
    from the control in each tube. Gel quality describes band separation only.
    """
    if ref not in {'A', 'C', 'G', 'T'} or alt not in {'A', 'C', 'G', 'T'} or ref == alt:
        return []
    if not 0 < offset < len(sequence) - 1 or sequence[offset].upper() != ref:
        return []
    if mask_repeats and sequence[offset].islower():
        return []
    if any(start <= offset < start + length for start, length in excluded_regions):
        return []
    template = sequence if mask_repeats else sequence.upper()
    alt_template = template[:offset] + alt + template[offset + 1:]
    assays = []
    for pair in control_pairs:
        check_cancel()
        start, end = pair['start0'], pair['end0']
        size = end - start
        if not (prod_min <= size <= prod_max and model.min_visible <= size <= model.max_visible):
            continue
        if set(template[start:end].upper()) - set('ACGT'):
            continue
        control = dict(pair, product_size=size)
        min_size, max_size = ceil(model.min_visible), min(size - 1, floor(model.max_visible))
        if min_size > max_size:
            continue
        for side in ('left', 'right'):
            options = dict(pick_primers=pick_primers, parse_pairs=parse_pairs,
                           min_size=min_size, max_size=max_size, tm_min=tm_min, tm_max=tm_max,
                           excluded_regions=excluded_regions, mask_repeats=mask_repeats)
            check_cancel()
            ref_pairs = _nested_pairs(template, offset, side, control, **options)
            if not ref_pairs:
                continue
            check_cancel()
            alt_pairs = _nested_pairs(alt_template, offset, side, control, **options)
            for ref_pair, alt_pair in product(ref_pairs, alt_pairs):
                score = _score(control, ref_pair, alt_pair, model)
                if score is not None:
                    assays.append(dict(score, aspcr_control=control, aspcr_specific_side=side,
                                       aspcr_reactions={'REF': dict(ref_pair, allele=ref),
                                                        'ALT': dict(alt_pair, allele=alt)}))
    return assays
