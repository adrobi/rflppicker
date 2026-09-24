"""Create SYNTHETIC smoke-test data, not biological validation data."""
import random
from pathlib import Path


def main():
    rng = random.Random(20260924)
    root = Path(__file__).parent
    fasta, variants = [], ['chrom\tpos\tref\talt']
    for i, site in enumerate(('GAATTC', 'GGATCC', 'AAGCTT', 'CTGCAG'), 1):
        seq = ''.join(rng.choice('ACGT') for _ in range(601))
        seq = seq[:298] + site + seq[304:]
        name = f'SYNTHETIC_{i}'
        ref = seq[300]
        alt = next(x for x in 'ACGT' if x != ref)
        fasta.extend([f'>{name}', seq])
        variants.append(f'{name}\t301\t{ref}\t{alt}')
    (root / 'synthetic.fa').write_text('\n'.join(fasta) + '\n', encoding='ascii')
    (root / 'synthetic.tsv').write_text('\n'.join(variants) + '\n', encoding='ascii')


if __name__ == '__main__':
    main()
