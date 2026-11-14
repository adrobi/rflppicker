import os
import random
import tempfile
import unittest

import specificity
from specificity import Product, Site

COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def revcomp(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


def mutate(text: str, positions, rng) -> str:
    bases = list(text)
    for position in positions:
        bases[position] = rng.choice([base for base in 'ACGT' if base != bases[position]])
    return ''.join(bases)


def brute_force(contigs: dict, primers, max_mismatches: int) -> dict:
    """Every site by plain comparison: the reference the scan must agree with."""
    seed = specificity.SEED
    found = {primer: [] for primer in primers}
    for primer in primers:
        for strand, pattern in (('+', primer), ('-', revcomp(primer))):
            exact = slice(len(pattern) - seed, None) if strand == '+' else slice(0, seed)
            for name, sequence in contigs.items():
                upper = sequence.upper()
                for start in range(len(upper) - len(pattern) + 1):
                    window = upper[start:start + len(pattern)]
                    if window[exact] != pattern[exact]:
                        continue
                    mismatches = sum(a != b for a, b in zip(window, pattern))
                    if mismatches <= max_mismatches:
                        found[primer].append(Site(name, start, strand, mismatches))
    return {primer: sorted(sites) for primer, sites in found.items()}


class FindSitesTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.addCleanup(self._dir.cleanup)

    def fasta(self, contigs: dict, width: int = 60) -> str:
        path = os.path.join(self._dir.name, f'g{len(os.listdir(self._dir.name))}.fa')
        with open(path, 'w') as stream:
            for name, sequence in contigs.items():
                stream.write(f'>{name}\n')
                stream.writelines(sequence[i:i + width] + '\n' for i in range(0, len(sequence), width))
        return path

    def scan(self, contigs, primers, **options):
        found = specificity.find_sites(self.fasta(contigs), primers, **options)
        return {primer: sorted(result.sites) for primer, result in found.items()}

    def test_planted_sites_on_both_strands_with_mismatches(self):
        rng = random.Random(1)
        primer = 'ACGTTGCATGCCATGGATCCA'                      # 21 nt
        genome = [rng.choice('ACGT') for _ in range(5000)]

        def plant(at, text):
            genome[at:at + len(text)] = text

        plant(100, primer)                                    # exact, + strand
        plant(700, revcomp(primer))                           # exact, - strand
        plant(1500, mutate(primer, (0, 4), rng))              # two mismatches in the 5' part
        plant(2300, mutate(primer, (0, 2, 4, 6), rng))        # four: too many
        plant(3100, mutate(primer, (len(primer) - 3,), rng))  # one inside the 3' seed: not a site
        plant(3900, revcomp(mutate(primer, (1, 3, 5), rng)))  # three, - strand
        sites = self.scan({'c1': ''.join(genome)}, [primer], max_mismatches=3)[primer]
        self.assertEqual(sites, [Site('c1', 100, '+', 0), Site('c1', 700, '-', 0), Site('c1', 1500, '+', 2),
                                 Site('c1', 3900, '-', 3)])

    def test_an_unknown_base_in_the_seed_is_not_a_site(self):
        rng = random.Random(4)
        primer = 'ACGTTGCATGCCATGGATTCA'                      # the 3' seed has several T
        genome = [rng.choice('ACGT') for _ in range(600)]
        for at, strand in ((100, '+'), (300, '-')):
            text = primer if strand == '+' else revcomp(primer)
            genome[at:at + len(text)] = text
        clean = ''.join(genome)
        self.assertEqual(len(self.scan({'c': clean}, [primer])[primer]), 2)
        for at in (100 + len(primer) - 3, 100 + len(primer) - 2, 300 + 1, 300 + 2):   # inside either seed
            broken = clean[:at] + 'N' + clean[at + 1:]
            with self.subTest(at=at):
                self.assertEqual(len(self.scan({'c': broken}, [primer])[primer]), 1)
        outside = clean[:102] + 'N' + clean[103:]             # in the 5' part: one mismatch, still a site
        self.assertEqual(self.scan({'c': outside}, [primer])[primer][0], Site('c', 100, '+', 1))

    def test_agrees_with_brute_force_across_block_borders(self):
        rng = random.Random(7)
        contigs = {f'c{i}': ''.join(rng.choice('ACGTacgt') for _ in range(length))
                   for i, length in enumerate((3000, 37, 1200, 11))}
        contigs['c0'] = contigs['c0'][:1000] + 'N' * 30 + contigs['c0'][1030:]
        primers = []
        for _ in range(40):                                   # primers cut from the genome, some then mutated
            name = rng.choice(['c0', 'c2'])
            start = rng.randrange(len(contigs[name]) - 30)
            primer = contigs[name][start:start + rng.randint(18, 27)].upper().replace('N', 'A')
            if rng.random() < 0.5:
                primer = revcomp(primer)
            if rng.random() < 0.5:
                primer = mutate(primer, rng.sample(range(len(primer) - specificity.SEED), 2), rng)
            primers.append(primer)
        # Repeats: the same 40 nt several times, also over a block border.
        repeat = contigs['c2'][100:140]
        contigs['c2'] = contigs['c2'][:500] + repeat + contigs['c2'][540:1000] + repeat + contigs['c2'][1040:]
        primers.append(repeat[5:27].upper())
        expected = brute_force(contigs, set(primers), 3)
        for block in (64, 257, 1 << 20):                      # block sizes around and above the contig sizes
            with self.subTest(block=block):
                self.assertEqual(self.scan(contigs, primers, max_mismatches=3, block=block, threads=3), expected)
        self.assertGreater(sum(map(len, expected.values())), len(primers))

    def test_lower_case_primer_letters_and_duplicates(self):
        rng = random.Random(3)
        genome = ''.join(rng.choice('ACGT') for _ in range(2000))
        primer = genome[300:322]
        found = specificity.find_sites(self.fasta({'c': genome}), [primer.lower(), primer, primer[:10] + primer[10:].lower()])
        self.assertEqual(list(found), [primer])               # one entry, in upper case
        self.assertEqual(found[primer].sites, [Site('c', 300, '+', 0)])
        self.assertFalse(found[primer].saturated)

    def test_unusable_primers_are_skipped(self):
        found = specificity.find_sites(self.fasta({'c': 'ACGT' * 100}), ['ACGTN' * 5, 'ACGTAC', ''])
        self.assertEqual(found, {})

    def test_a_repeat_saturates_the_primer(self):
        rng = random.Random(5)
        unit = ''.join(rng.choice('ACGT') for _ in range(30))
        spacer = ''.join(rng.choice('ACGT') for _ in range(20))
        genome = (unit + spacer) * 40
        found = specificity.find_sites(self.fasta({'c': genome}), [unit[2:24]], max_sites=10)[unit[2:24]]
        self.assertTrue(found.saturated)
        self.assertEqual(len(found.sites), 10)

    def test_progress_and_cancel(self):
        rng = random.Random(9)
        genome = ''.join(rng.choice('ACGT') for _ in range(4000))
        path = self.fasta({'c': genome})
        seen = []
        specificity.find_sites(path, [genome[50:72]], block=512, progress_cb=lambda done, total: seen.append((done, total)))
        self.assertEqual(seen[-1], (len(seen), len(seen)))
        self.assertGreater(len(seen), 4)
        with self.assertRaises(specificity.ScanCancelled):
            specificity.find_sites(path, [genome[50:72]], block=512, cancel_cb=lambda: True)


class ProductTests(unittest.TestCase):
    L, R = 'A' * 20, 'C' * 22                                 # only the lengths matter here

    def sites(self, left=(), right=()):
        return {self.L: specificity.PrimerSites(list(left)), self.R: specificity.PrimerSites(list(right))}

    def test_facing_sites_within_the_size_limit(self):
        found = self.sites(left=[Site('c', 100, '+', 0), Site('c', 5000, '+', 1)],
                           right=[Site('c', 478, '-', 0), Site('c', 9000, '-', 0), Site('d', 300, '-', 0)])
        products = specificity.pair_products(found, self.L, self.R, max_size=2000)
        self.assertEqual(products, [Product('c', 100, 500, 'F+R')])     # 478 + 22 = 500, a product of 400 bp
        self.assertEqual(products[0].size, 400)

    def test_one_primer_on_both_sides(self):
        found = self.sites(left=[Site('c', 100, '+', 0), Site('c', 600, '-', 2)])
        self.assertEqual(specificity.pair_products(found, self.L, self.R, max_size=2000),
                         [Product('c', 100, 620, 'F+F')])

    def test_sites_facing_away_or_too_far_give_nothing(self):
        found = self.sites(left=[Site('c', 1000, '-', 0)], right=[Site('c', 100, '+', 0), Site('c', 5000, '-', 0)])
        # R(+) at 100 and L(-) at 1000 face each other: 100..1020; R(-) at 5000 is too far from R(+).
        self.assertEqual(specificity.pair_products(found, self.L, self.R, max_size=2000),
                         [Product('c', 100, 1020, 'R+F')])
        self.assertEqual(specificity.pair_products(found, self.L, self.R, max_size=500), [])

    def test_the_intended_product_is_left_out(self):
        found = self.sites(left=[Site('c', 100, '+', 0), Site('c', 3000, '+', 0)],
                           right=[Site('c', 478, '-', 0), Site('c', 3378, '-', 0)])
        products = specificity.pair_products(found, self.L, self.R, max_size=2000, intended=('c', 100, 500))
        self.assertEqual(products, [Product('c', 3000, 3400, 'F+R')])

    def test_overlapping_sites_are_not_a_product(self):
        found = self.sites(left=[Site('c', 100, '+', 0)], right=[Site('c', 105, '-', 0)])
        self.assertEqual(specificity.pair_products(found, self.L, self.R, max_size=2000), [])

    def test_unknown_primer_has_no_products(self):
        self.assertEqual(specificity.pair_products({}, self.L, self.R, max_size=2000), [])


if __name__ == '__main__':
    unittest.main()
