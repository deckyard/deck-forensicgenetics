# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for popstats.j. Run with `jennifer test popstats_test.j`.
#
# Locus "L" is the kinship overlay's table: alleles 10/11/12/13 at
# .1/.2/.3/.4, N = 10000. Its power sums, used throughout:
#   S2 = .01 + .04 + .09 + .16           = .30
#   S3 = .001 + .008 + .027 + .064       = .100
#   S4 = .0001 + .0016 + .0081 + .0256   = .0354
#   S5 = .00001 + .00032 + .00243 + .01024 = .013
# Locus "M" is two alleles at .5/.5 (S2 = .5, S3 = .25, S4 = .125,
# S5 = .0625), the textbook SNP. Every expectation is written out from those
# sums or from the Levene-Haldane distribution by hand; the cross-checks
# recompute a statistic by a different route (genotype enumeration through the
# core, or the Levene formula through lgamma).
use testing;
use math;
use strings;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

# within asserts a Monte Carlo estimate is within `tol` of the exact value.
func within(actual as float, expected as float, tol as float) {
    testing.assertTrue(math.abs($actual - $expected) < $tol);
}

func fixture() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("test"),
        "L",
        10000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    return forensics.withLocus($out, "M", 500, {"1": 0.5, "2": 0.5});
}

func gl(locus as string, a as string, b as string) {
    return forensics.genotype($locus, $a, $b);
}

# --- frequency statistics ---------------------------------------------------

func testExpectedHeterozygosity() {
    # 1 - S2 = 1 - .30
    near(expectedHeterozygosity(fixture(), "L"), 0.7);
    near(expectedHeterozygosity(fixture(), "M"), 0.5);
}

func testUnbiasedHeterozygosity() {
    # 2N / (2N - 1) * (1 - S2) = 20000 / 19999 * .7
    near(unbiasedHeterozygosity(fixture(), "L"), 0.7 * 20000.0 / 19999.0);
    # M: 1000 / 999 * .5
    near(unbiasedHeterozygosity(fixture(), "M"), 0.5 * 1000.0 / 999.0);
}

func testMatchProbabilityAndDiscrimination() {
    # 2 S2^2 - S4 = 2 * .09 - .0354 = .1446; PD = .8554
    near(matchProbability(fixture(), "L"), 0.1446);
    near(powerOfDiscrimination(fixture(), "L"), 0.8554);
    # M: genotypes 11, 12, 22 at .25 / .5 / .25: .0625 + .25 + .0625 = .375
    near(matchProbability(fixture(), "M"), 0.375);
    near(powerOfDiscrimination(fixture(), "M"), 0.625);
}

func testMatchProbabilityEqualsGenotypeEnumeration() {
    # Cross-check: sum P(g)^2 over every genotype, each P(g) from the core's
    # Hardy-Weinberg genotypeFrequency at theta 0 (the 5/2N floor, .00025, is
    # below every allele, so it does not move anything).
    def fdb as forensics.FrequencyDb init fixture();
    def labels as list of string init ["10", "11", "12", "13"];
    def total as float init 0.0;
    for (def i in 0..len($labels)) {
        for (def j in $i..len($labels)) {
            def f as float init forensics.genotypeFrequency(
                $fdb,
                gl("L", $labels[$i], $labels[$j]),
                forensics.THETA_NONE);
            $total = $total + $f * $f;
        }
    }
    near(matchProbability($fdb, "L"), $total);
}

func testPolymorphismInformationContent() {
    # Botstein: 1 - S2 - sum_{i<j} 2 pi^2 pj^2. The pairs:
    #   2 * (.01*.04 + .01*.09 + .01*.16 + .04*.09 + .04*.16 + .09*.16)
    # = 2 * (.0004 + .0009 + .0016 + .0036 + .0064 + .0144) = 2 * .0273 = .0546
    # PIC = 1 - .3 - .0546 = .6454
    near(polymorphismInformationContent(fixture(), "L"), 0.6454);
    # M: 1 - .5 - 2 * .25 * .25 = .375 (Botstein's 3/8 for a biallelic 50/50)
    near(polymorphismInformationContent(fixture(), "M"), 0.375);
}

func testTypicalPaternityIndex() {
    # 1 / (2 S2) = 1 / .6
    near(typicalPaternityIndex(fixture(), "L"), 1.0 / 0.6);
    near(typicalPaternityIndex(fixture(), "M"), 1.0);
}

func testPowerOfExclusionTrio() {
    # 1 - 2 S2 + S3 + 2 S4 - 3 S5 - 2 S2^2 + 3 S2 S3
    # = 1 - .6 + .1 + .0708 - .039 - .18 + .09 = .4418
    near(powerOfExclusionTrio(fixture(), "L"), 0.4418);
    # M: 1 - 1 + .25 + .25 - .1875 - .5 + .375 = .1875 = 3/16
    near(powerOfExclusionTrio(fixture(), "M"), 0.1875);
}

func testPowerOfExclusionDuo() {
    # 1 - 4 S2 + 2 S2^2 + 4 S3 - 3 S4 = 1 - 1.2 + .18 + .4 - .1062 = .2738
    near(powerOfExclusionDuo(fixture(), "L"), 0.2738);
    # M: 1 - 2 + .5 + 1 - .375 = .125. By hand: a duo excludes only a man
    # sharing no allele with the child, which needs a homozygous child (.5)
    # and a man homozygous for the other allele (.25): .5 * .25 = .125.
    near(powerOfExclusionDuo(fixture(), "M"), 0.125);
}

func testMonomorphicLocusHasNoPower() {
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("mono"),
        "X",
        100,
        {"7": 1.0});
    def row as LocusStats init locusStats($fdb, "X");
    near($row.expectedHeterozygosity, 0.0);
    near($row.matchProbability, 1.0);
    near($row.pic, 0.0);
    near($row.typicalPaternityIndex, 0.5);
    near($row.powerOfExclusionTrio, 0.0);
    near($row.powerOfExclusionDuo, 0.0);
}

func testRoundedTableIsRescaled() {
    # .099/.198/.297/.396 sums to .99; rescaled it is exactly .1/.2/.3/.4.
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("rounded"),
        "L",
        10000,
        {"10": 0.099, "11": 0.198, "12": 0.297, "13": 0.396});
    near(expectedHeterozygosity($fdb, "L"), 0.7);
    near(matchProbability($fdb, "L"), 0.1446);
}

func testLocusStatsGathersEveryStatistic() {
    def row as LocusStats init locusStats(fixture(), "L");
    testing.assertEqual($row.locus, "L");
    testing.assertEqual($row.n, 10000);
    testing.assertEqual($row.alleleCount, 4);
    near($row.expectedHeterozygosity, 0.7);
    near($row.matchProbability, 0.1446);
    near($row.powerOfDiscrimination, 0.8554);
    near($row.pic, 0.6454);
    near($row.typicalPaternityIndex, 1.0 / 0.6);
    near($row.powerOfExclusionTrio, 0.4418);
    near($row.powerOfExclusionDuo, 0.2738);
}

func testStatsTableFollowsLocusOrder() {
    def table as list of LocusStats init statsTable(fixture());
    testing.assertEqual(len($table), 2);
    testing.assertEqual($table[0].locus, "L");
    testing.assertEqual($table[1].locus, "M");
    near($table[1].powerOfExclusionTrio, 0.1875);
}

# --- combining loci ---------------------------------------------------------

func testCombinedStats() {
    def all as CombinedStats init combinedStats(statsTable(fixture()));
    testing.assertEqual($all.loci, 2);
    # MP: .1446 * .375 = .054225; PD = .945775
    near($all.matchProbability, 0.054225);
    near($all.powerOfDiscrimination, 0.945775);
    # trio: 1 - (1 - .4418)(1 - .1875) = 1 - .5582 * .8125 = 1 - .4535375
    near($all.powerOfExclusionTrio, 0.5464625);
    # duo: 1 - (1 - .2738)(1 - .125) = 1 - .7262 * .875 = 1 - .635425
    near($all.powerOfExclusionDuo, 0.364575);
    # TPI: (1 / .6) * 1
    near($all.typicalPaternityIndex, 1.0 / 0.6);
    near($all.log10MatchProbability, math.log10(0.054225));
    near($all.log10NonExclusionTrio, math.log10(0.4535375));
    near($all.log10NonExclusionDuo, math.log10(0.635425));
    near($all.log10TypicalPaternityIndex, math.log10(1.0 / 0.6));
}

func testCombinedLogDoesNotUnderflow() {
    # 400 copies of locus L: MP = .1446^400 ~ 10^-336, below the smallest
    # double, so the product is 0 while the log stays 400 * log10(.1446).
    def row as LocusStats init locusStats(fixture(), "L");
    def table as list of LocusStats init [];
    for (def i in 0..400) {
        $table[] = $row;
    }
    def all as CombinedStats init combinedStats($table);
    near($all.log10MatchProbability, 400.0 * math.log10(0.1446));
    testing.assertTrue($all.matchProbability < 1e-300);
    near($all.powerOfDiscrimination, 1.0);
}

func testCombinedPower() {
    # 1 - (1 - .4418)(1 - .1875), as in testCombinedStats
    near(combinedPower([0.4418, 0.1875]), 0.5464625);
    near(log10CombinedComplement([0.4418, 0.1875]), math.log10(0.4535375));
    near(combinedPower([0.0]), 0.0);
}

# --- observed genotypes -----------------------------------------------------

# sampleA, five people at L: 10/10, 10/11, 11/12, 12/12, 10/12.
#   heterozygotes 3 of 5 -> H_obs = .6, homozygosity .4
#   copies 10: 2 + 1 + 1 = 4, 11: 1 + 1 = 2, 12: 1 + 2 + 1 = 4 (sum 10 = 2n)
func sampleA() {
    return [
        gl("L", "10", "10"),
        gl("L", "10", "11"),
        gl("L", "11", "12"),
        gl("L", "12", "12"),
        gl("L", "10", "12")
    ];
}

func testObservedHeterozygosity() {
    near(observedHeterozygosity(sampleA()), 0.6);
}

func testAlleleCountsAndFrequencies() {
    def counts as map of string to int init alleleCounts(sampleA());
    testing.assertEqual($counts, {"10": 4, "11": 2, "12": 4});
    def freqs as map of string to float init estimateFrequencies(sampleA());
    near($freqs["10"], 0.4);
    near($freqs["11"], 0.2);
    near($freqs["12"], 0.4);
}

func testWithEstimatedLocus() {
    def fdb as forensics.FrequencyDb init withEstimatedLocus(forensics.db("sample"), sampleA());
    testing.assertEqual(forensics.sampleSize($fdb, "L"), 5);
    # 1 - (.16 + .04 + .16) = .64; Nei's correction 10/9: .64 * 10 / 9
    near(expectedHeterozygosity($fdb, "L"), 0.64);
    near(unbiasedHeterozygosity($fdb, "L"), 0.64 * 10.0 / 9.0);
}

func testObservedMatchProbability() {
    # five distinct genotypes, each 1/5: 5 * (1/5)^2 = .2
    near(observedMatchProbability(sampleA()), 0.2);
    # AA, AB, AB, BB: (1/4)^2 + (2/4)^2 + (1/4)^2 = .375, the HWE value for
    # .5/.5 because this sample is in exact HWE proportions
    def hwe as list of forensics.Genotype init [
        gl("M", "1", "1"),
        gl("M", "1", "2"),
        gl("M", "1", "2"),
        gl("M", "2", "2")
    ];
    near(observedMatchProbability($hwe), matchProbability(fixture(), "M"));
}

func testObservedTypicalPaternityIndex() {
    # homozygosity .4: 1 / (2 * .4) = 1.25
    near(observedTypicalPaternityIndex(sampleA()), 1.25);
}

func testObservedPowerOfExclusion() {
    # H = .6, h = .4: .36 * (1 - 2 * .6 * .16) = .36 * .808 = .29088
    near(observedPowerOfExclusion(sampleA()), 0.29088);
    # An HWE .5/.5 sample (H = .5): .25 * (1 - 2 * .5 * .25) = .1875, the
    # exact trio value of locus M.
    def hwe as list of forensics.Genotype init [
        gl("M", "1", "1"),
        gl("M", "1", "2"),
        gl("M", "1", "2"),
        gl("M", "2", "2")
    ];
    near(observedPowerOfExclusion($hwe), powerOfExclusionTrio(fixture(), "M"));
}

# --- Hardy-Weinberg ---------------------------------------------------------

# Levene-Haldane for n = 4 people, allele counts A 4, B 4 (2n = 8):
#   P(H) = 4! 4! 4! 2^H / (8! nAA! nAB! nBB!),  4!^3 = 13824, 8! = 40320
#   H = 0 (AA 2, BB 2):  13824 / (40320 * 2 * 2)           = 3/35
#   H = 2 (1, 2, 1):     13824 * 4 / (40320 * 2)           = 24/35
#   H = 4 (0, 4, 0):     13824 * 16 / (40320 * 24)         = 8/35
# p(H = 4) = 8/35 + 3/35 = 11/35; p(H = 0) = 3/35; p(H = 2) = 1.
func allHet() {
    return [gl("M", "1", "2"), gl("M", "1", "2"), gl("M", "1", "2"), gl("M", "1", "2")];
}

func testHweExactBiallelic() {
    near(hweExactTest(allHet()), 11.0 / 35.0);
    near(
        hweExactTest([gl("M", "1", "1"), gl("M", "1", "1"), gl("M", "2", "2"), gl("M", "2", "2")]),
        3.0 / 35.0);
    near(
        hweExactTest([gl("M", "1", "1"), gl("M", "1", "2"), gl("M", "2", "1"), gl("M", "2", "2")]),
        1.0);
    # n = 3, A 2, B 4: P(H = 0) = 3! 2! 4! / (6! 1! 2!) = 288 / 1440 = .2,
    # P(H = 2) = 288 * 4 / (720 * 2) = .8; observing H = 0 gives .2
    near(hweExactTest([gl("M", "1", "1"), gl("M", "2", "2"), gl("M", "2", "2")]), 0.2);
}

func testHweExactMonomorphic() {
    near(hweExactTest([gl("M", "1", "1"), gl("M", "1", "1")]), 1.0);
}

# leveneP sums, by the Levene formula through lgamma, every heterozygote count
# no more likely than `het` for n people with `rare` copies of the minor
# allele - a second route to the p-value that shares nothing with the
# Wigginton recurrence.
func leveneLog(n as int, rare as int, h as int) {
    def homr as int init ($rare - $h) // 2;
    def homc as int init $n - $h - $homr;
    return math.lgamma(convert.toFloat($n + 1)) + math.lgamma(convert.toFloat($rare + 1)) +
        math.lgamma(convert.toFloat(2 * $n - $rare + 1)) +
        convert.toFloat($h) *
        math.ln(2.0) -
        math.lgamma(convert.toFloat(2 * $n + 1)) -
        math.lgamma(convert.toFloat($homr + 1)) -
        math.lgamma(convert.toFloat($h + 1)) - math.lgamma(convert.toFloat($homc + 1));
}

func leveneP(n as int, rare as int, het as int) {
    def obs as float init math.exp(leveneLog($n, $rare, $het));
    def p as float init 0.0;
    def h as int init $rare % 2;
    while ($h <= $rare) {
        def x as float init math.exp(leveneLog($n, $rare, $h));
        if ($x <= $obs * 1.0000001) {
            $p = $p + $x;
        }
        $h = $h + 2;
    }
    return $p;
}

func testHweExactAgreesWithLeveneFormula() {
    # 100 people, 60 copies of the minor allele, 20 heterozygotes
    # (20 homozygous minor, 60 homozygous major): a heterozygote deficit.
    def g as list of forensics.Genotype init [];
    for (def i in 0..20) {
        $g[] = gl("M", "1", "1");
        $g[] = gl("M", "1", "2");
    }
    for (def i in 0..60) {
        $g[] = gl("M", "2", "2");
    }
    def p as float init hweExactTest($g);
    testing.assertTrue(math.abs($p - leveneP(100, 60, 20)) < 1e-9);
    testing.assertTrue($p < 0.001);
}

# Three people, alleles a, b, c twice each (2n = 6). Levene:
#   P = 3! 2! 2! 2! 2^H / (6! Π n_ij!) = (1/15) 2^H / Π n_ij!
#   aa bb cc            H 0              1/15
#   aa bc bc (x3 ways)  H 2, bc twice    4 / 2 / 15 = 2/15 each
#   ab ac bc            H 3              8/15
# Total 15/15. p(aa bb cc) = 1/15, p(aa bc bc) = 1/15 + 6/15 = 7/15,
# p(ab ac bc) = 1.
func threeHoms() {
    return [gl("T", "a", "a"), gl("T", "b", "b"), gl("T", "c", "c")];
}

func testHwePermutationMultiallelic() {
    # 10000 permutations: SE ~ sqrt(p(1-p)/10000) = .0025 at 1/15, .005 at
    # 7/15; tolerances are 4 SE. The seed fixes the draw, so these are stable.
    within(hwePermutationTest(threeHoms(), 10000, 12345), 1.0 / 15.0, 0.01);
    within(
        hwePermutationTest([gl("T", "a", "a"), gl("T", "b", "c"), gl("T", "b", "c")], 10000, 777),
        7.0 / 15.0,
        0.02);
    # The most likely table: every shuffle is no more likely, so r = B and
    # (B + 1) / (B + 1) = 1 exactly.
    near(
        hwePermutationTest([gl("T", "a", "b"), gl("T", "a", "c"), gl("T", "b", "c")], 500, 9),
        1.0);
}

func testHwePermutationEstimatesTheExactTest() {
    # 11/35 = .3143, SE .0046 at 10000 permutations
    within(hwePermutationTest(allHet(), 10000, 2026), hweExactTest(allHet()), 0.02);
}

func testHwePermutationIsReproducible() {
    near(hwePermutationTest(threeHoms(), 2000, 42), hwePermutationTest(threeHoms(), 2000, 42));
    near(hwePermutationTest([gl("M", "1", "1")], 10, 1), 1.0);
}

func testGeneratorIsParkMiller() {
    # Park & Miller (1988), multiplier 48271, modulus 2^31 - 1. From seed 1:
    #   48271
    #   48271^2 = 2330089441; - 2147483647 = 182605794
    # ISO C++ [rand.predef], std::minstd_rand (a = 48271, seed 1): the 10000th
    # state is 399268537.
    testing.assertEqual(nextState(1), 48271);
    testing.assertEqual(nextState(48271), 182605794);
    def x as int init 1;
    for (def i in 0..10000) {
        $x = nextState($x);
    }
    testing.assertEqual($x, 399268537);
}

func testRedrawRejectsAboveTheLimit() {
    # From 1914720637 the generator gives 2078669041, then 407355683 (both
    # by 48271 x mod 2^31 - 1). Under a limit of 1e9 the first is rejected.
    testing.assertEqual(redraw(1914720637, 1000000000), 407355683);
    testing.assertEqual(redraw(1, 48271), 48271);
}

# --- F-statistics -----------------------------------------------------------

# Every expectation below is Weir & Cockerham's (1984) a, b, c worked by hand
# per allele: n̄, n_c, p̄, s², h̄, then
#   pq = p̄ (1 - p̄) - (r - 1) s² / r
#   a  = n̄ / n_c · [s² - (pq - h̄ / 4) / (n̄ - 1)]
#   b  = n̄ / (n̄ - 1) · [pq - (2 n̄ - 1) h̄ / (4 n̄)]
#   c  = h̄ / 2

# copies is `count` people of genotype a/b at `locus`.
func copies(locus as string, a as string, b as string, count as int) {
    def out as list of forensics.Genotype init [];
    for (def i in 0..$count) {
        $out[] = gl($locus, $a, $b);
    }
    return $out;
}

func joined(parts as list of list of forensics.Genotype) {
    def out as list of forensics.Genotype init [];
    for (def i in 0..len($parts)) {
        for (def j in 0..len($parts[$i])) {
            $out[] = $parts[$i][$j];
        }
    }
    return $out;
}

# Two populations of two: {AA, AB} and {AB, BB}. r = 2, n̄ = 2,
# n_c = (4 - 8/4) / 1 = 2. Allele A: p = 3/4, 1/4, p̄ = 1/2,
# s² = (2 (1/4)² + 2 (1/4)²) / (1 · 2) = 1/8, h = 1/2, 1/2, h̄ = 1/2;
#   pq = 1/4 - 1/16 = 3/16
#   a  = 2/2 · [1/8 - (3/16 - 1/8) / 1] = 1/16
#   b  = 2/1 · [3/16 - (3/4)(1/2) / 2] = 0
#   c  = 1/4
# B mirrors A: a = 1/8, b = 0, c = 1/2 in all;
# fst = (1/8) / (5/8) = 1/5, fit = 1/5, fis = 0.
func twoByTwo() {
    return [[gl("A", "a", "a"), gl("A", "a", "b")], [gl("A", "a", "b"), gl("A", "b", "b")]];
}

func testFStatisticsTwoPopulationsByHand() {
    def s as FStatistics init fStatistics(twoByTwo());
    testing.assertEqual($s.locus, "A");
    testing.assertEqual($s.loci, 1);
    testing.assertEqual($s.populations, 2);
    near($s.a, 0.125);
    near($s.b, 0.0);
    near($s.c, 0.5);
    near(fst($s), 0.2);
    near(fit($s), 0.2);
    near(fis($s), 0.0);
}

# Unbalanced, three alleles: {aa, ab, bc} (n 3) and {cc, bc} (n 2).
# n̄ = 5/2, n_c = (5 - 13/5) / 1 = 12/5, n̄/n_c = 25/24, n̄/(n̄-1) = 5/3,
# (2n̄ - 1)/(4n̄) = 2/5.
# a: p = 1/2, 0; p̄ = 3/10; s² = (3 (1/5)² + 2 (3/10)²) / (5/2) = 3/25;
#    h̄ = (1 + 0) / 5 = 1/5; pq = 21/100 - 3/50 = 3/20
#    a = 25/24 (3/25 - (3/20 - 1/20)(2/3)) = 25/24 · 4/75 = 1/18
#    b = 5/3 (3/20 - 2/25) = 7/60;  c = 1/10
# b: p = 1/3, 1/4; p̄ = 3/10; s² = (3 (1/30)² + 2 (1/20)²) / (5/2) = 1/300;
#    h̄ = (2 + 1) / 5 = 3/5; pq = 21/100 - 1/600 = 5/24
#    a = 25/24 (1/300 - (5/24 - 3/20)(2/3)) = 25/24 · (-32/900) = -1/27
#    b = 5/3 (5/24 - 6/25) = -19/360;  c = 3/10
# c: p = 1/6, 3/4; p̄ = 2/5; s² = (3 (7/30)² + 2 (7/20)²) / (5/2) = 49/300;
#    h̄ = (1 + 1) / 5 = 2/5; pq = 6/25 - 49/600 = 19/120
#    a = 25/24 (49/300 - (19/120 - 1/10)(2/3)) = 25/24 · 112/900 = 7/54
#    b = 5/3 (19/120 - 4/25) = -1/360;  c = 1/5
# Sums: a = 4/27, b = 11/180, c = 3/5; over 540: 80 + 33 + 324 = 437, so
# fst = 80/437, fit = 113/437, fis = 33/357 = 11/119.
func testFStatisticsUnbalancedByHand() {
    def s as FStatistics init fStatistics([
        [gl("U", "a", "a"), gl("U", "a", "b"), gl("U", "b", "c")],
        [gl("U", "c", "c"), gl("U", "b", "c")]
    ]);
    near($s.a, 4.0 / 27.0);
    near($s.b, 11.0 / 180.0);
    near($s.c, 0.6);
    near(fst($s), 80.0 / 437.0);
    near(fit($s), 113.0 / 437.0);
    near(fis($s), 11.0 / 119.0);
}

# Identical samples: s² = 0, so per allele a = -(pq - h̄/4) / (n̄ - 1) · n̄/n_c
# and a + b + c = pq; with n_c = n̄ = n, θ = -(Σ pq - Σ h̄/4) / ((n - 1) Σ pq):
# a negative bias of order 1/n, as an unbiased-components ratio should show
# when the populations do not differ.
# Two samples of {3 aa, 4 ab, 3 bb}: p̄ = 1/2, h̄ = 2/5 per allele,
#   n = 10:  θ = -(1/2 - 1/5) / (9 · 1/2) = -1/15
#   n = 100 (30/40/30): θ = -(3/10) / (99/2) = -1/165
# b at n = 10: 10/9 (1/4 - (19/40)(2/5)) = 1/15 per allele, c = 1/5:
#   fis = (2/15) / (2/15 + 2/5) = 1/4, fit = (-1/30 + 2/15) / (1/2) = 1/5.
func identical(count as int) {
    def one as list of forensics.Genotype init joined([
        copies("I", "a", "a", 3 * $count),
        copies("I", "a", "b", 4 * $count),
        copies("I", "b", "b", 3 * $count)
    ]);
    return [$one, $one];
}

func testIdenticalPopulationsGiveThetaNearZero() {
    def small as FStatistics init fStatistics(identical(1));
    near(fst($small), -1.0 / 15.0);
    near(fis($small), 0.25);
    near(fit($small), 0.2);
    near(fst(fStatistics(identical(10))), -1.0 / 165.0);
}

# Populations fixed for different alleles: per allele p̄ = 1/2, s² = 1/2,
# h̄ = 0, pq = 1/4 - 1/4 = 0; a = 2/2 (1/2 - 0) = 1/2, b = c = 0.
# fst = fit = 1; fis = 0/0, undefined.
func fixedApart() {
    return fStatistics([
        [gl("X", "a", "a"), gl("X", "a", "a")],
        [gl("X", "b", "b"), gl("X", "b", "b")]
    ]);
}

func testFixedPopulations() {
    near(fst(fixedApart()), 1.0);
    near(fit(fixedApart()), 1.0);
    testing.assertThrows("fisOfFixed", "popstats");
}

func fisOfFixed() {
    fis(fixedApart());
}

# One sample {aa, ab, bb, bb}, n = 4, r = 1 (s² terms vanish, a = 0):
# a: p = 3/8, h = 1/4; b = 4/3 (15/64 - (7/16)(1/4)) = 4/3 · 8/64 = 1/6
# b: p = 5/8, h = 1/4; same b = 1/6; c = 1/8 each.
# f = (1/3) / (1/3 + 1/4) = 4/7. By the closed form, Ho = 1/4,
# He = 1 - 9/64 - 25/64 = 15/32: 1 - 2·3·(1/4) / (8·15/32 - 1/4) = 1 - 3/7.
func testSinglePopulationInbreeding() {
    def one as list of forensics.Genotype init [
        gl("S", "a", "a"),
        gl("S", "a", "b"),
        gl("S", "b", "b"),
        gl("S", "b", "b")
    ];
    def s as FStatistics init fStatistics([$one]);
    testing.assertEqual($s.populations, 1);
    near($s.a, 0.0);
    near($s.b, 1.0 / 3.0);
    near($s.c, 0.25);
    near(fis($s), 4.0 / 7.0);
    near(inbreedingCoefficient($one), 4.0 / 7.0);
    testing.assertThrows("fstOfOne", "popstats");
    testing.assertThrows("fitOfOne", "popstats");
}

func fstOfOne() {
    fst(fStatistics([[gl("S", "a", "a"), gl("S", "a", "b")]]));
}

func fitOfOne() {
    fit(fStatistics([[gl("S", "a", "a"), gl("S", "a", "b")]]));
}

func testInbreedingCoefficientClosedForm() {
    # 1 - 2 (n - 1) Ho / (2 n He - Ho), with Ho and He from the module's own
    # counting functions, on sampleA (3 alleles, n = 5).
    def g as list of forensics.Genotype init sampleA();
    def ho as float init observedHeterozygosity($g);
    def he as float init 0.0;
    def freqs as map of string to float init estimateFrequencies($g);
    for (def a in $freqs) {
        $he = $he + $freqs[$a] * $freqs[$a];
    }
    $he = 1.0 - $he;
    def n as float init convert.toFloat(len($g));
    near(inbreedingCoefficient($g), 1.0 - 2.0 * ($n - 1.0) * $ho / (2.0 * $n * $he - $ho));
    # an all-heterozygous sample is a heterozygote excess: f < 0
    testing.assertTrue(inbreedingCoefficient(allHet()) < 0.0);
}

# Genepop documentation, section "ANOVA estimators: single- and multilocus
# definitions" (genepop R package 1.2.17, inst/doc/GenepopS.tex), example file
# Fmulti.txt (3 samples, 3 loci; data from A.J. Gharrett). Genotype counts
# read from that file, its published Weir & Cockerham estimates to 4 places:
#   locus 1   FIS -0.0483  FST  0.5712  FIT 0.5505
#   locus 2   FIS -0.1161  FST  0.8560  FIT 0.8393
#   locus 3   FIS  0.0051  FST -0.0023  FIT 0.0028
#   multilocus (1984 a, b, c method) FIS -0.0286  FST 0.5606  FIT 0.5480
func fmulti() {
    def loc1 as list of list of forensics.Genotype init [
        joined([copies("Loc1", "01", "01", 64), copies("Loc1", "01", "02", 16)]),
        joined([copies("Loc1", "01", "02", 8), copies("Loc1", "02", "02", 32)]),
        joined([
            copies("Loc1", "01", "01", 10),
            copies("Loc1", "01", "02", 20),
            copies("Loc1", "02", "02", 10)
        ])
    ];
    def loc2 as list of list of forensics.Genotype init [
        joined([copies("Loc2", "01", "01", 32), copies("Loc2", "01", "02", 8)]),
        joined([copies("Loc2", "01", "02", 10), copies("Loc2", "02", "02", 30)]),
        copies("Loc2", "03", "03", 40)
    ];
    def loc3 as list of list of forensics.Genotype init [
        joined([
            copies("Loc3", "01", "01", 15),
            copies("Loc3", "01", "02", 30),
            copies("Loc3", "02", "02", 15)
        ]),
        joined([
            copies("Loc3", "01", "01", 8),
            copies("Loc3", "01", "02", 20),
            copies("Loc3", "02", "02", 12)
        ]),
        joined([
            copies("Loc3", "01", "01", 12),
            copies("Loc3", "01", "02", 20),
            copies("Loc3", "02", "02", 8)
        ])
    ];
    return [fStatistics($loc1), fStatistics($loc2), fStatistics($loc3)];
}

func fourPlaces(actual as float, expected as float) {
    within($actual, $expected, 0.00005001);
}

func testGenepopFmultiExample() {
    def t as list of FStatistics init fmulti();
    def want as list of list of float init [
        [-0.0483, 0.5712, 0.5505],
        [-0.1161, 0.8560, 0.8393],
        [0.0051, -0.0023, 0.0028]
    ];
    for (def i in 0..3) {
        fourPlaces(fis($t[$i]), $want[$i][0]);
        fourPlaces(fst($t[$i]), $want[$i][1]);
        fourPlaces(fit($t[$i]), $want[$i][2]);
    }
    def all as FStatistics init combineFStatistics($t);
    testing.assertEqual($all.locus, "");
    testing.assertEqual($all.loci, 3);
    fourPlaces(fis($all), -0.0286);
    fourPlaces(fst($all), 0.5606);
    fourPlaces(fit($all), 0.5480);
}

func testCombineSumsComponents() {
    def x as FStatistics init fStatistics(twoByTwo());
    def y as FStatistics init fixedApart();
    def both as FStatistics init combineFStatistics([$x, $y]);
    # a = 1/8 + 1, b = 0, c = 1/2: fst = (9/8) / (13/8) = 9/13
    near($both.a, 1.125);
    near($both.c, 0.5);
    near(fst($both), 9.0 / 13.0);
    near(fis($both), 0.0);
}

# Locus A is twoByTwo; locus B: population 0 {12, 22, 11}, population 1
# {12, 12} plus a third profile untyped at B. The table must equal
# fStatistics on those lists, and population 1's untyped profile must drop out.
func profiles2() {
    def pop0 as list of forensics.Profile init [
        forensics.profile("p0", [gl("A", "a", "a"), gl("B", "1", "2")]),
        forensics.profile("p1", [gl("A", "a", "b"), gl("B", "2", "2")]),
        forensics.profile("p2", [gl("B", "1", "1")])
    ];
    def pop1 as list of forensics.Profile init [
        forensics.profile("q0", [gl("A", "a", "b"), gl("B", "1", "2")]),
        forensics.profile("q1", [gl("A", "b", "b"), gl("B", "1", "2")]),
        forensics.profile("q2", [gl("C", "5", "6")])
    ];
    return [$pop0, $pop1];
}

func testFStatisticsTableFromProfiles() {
    def t as list of FStatistics init fStatisticsTable(profiles2(), ["A", "B"]);
    testing.assertEqual(len($t), 2);
    testing.assertEqual($t[0].locus, "A");
    near(fst($t[0]), 0.2);
    def b as FStatistics init fStatistics([
        [gl("B", "1", "2"), gl("B", "2", "2"), gl("B", "1", "1")],
        [gl("B", "1", "2"), gl("B", "1", "2")]
    ]);
    near($t[1].a, $b.a);
    near($t[1].b, $b.b);
    near($t[1].c, $b.c);
    # B by hand: n̄ = 5/2, n_c = 12/5, p̄ = 1/2 for both alleles, p = 1/2 in
    # each sample so s² = 0, h̄ = (1 + 2)/5 = 3/5, pq = 1/4:
    #   a = 25/24 (0 - (1/4 - 3/20)(2/3)) = -5/72,  b = 5/3 (1/4 - 6/25) = 1/60,
    #   c = 3/10;  x2 alleles: fst = (-5/36) / (-5/36 + 1/30 + 3/5) = -25/89
    near(fst($t[1]), -25.0 / 89.0);
}

# --- linkage disequilibrium -------------------------------------------------

# exactLd enumerates every permutation of the second locus's genotypes (n!)
# and returns the share whose table is no more probable than the observed
# one: the exact p-value the permutation test estimates. Table probability
# given margins is ∝ 1 / Π n_ij!.
func exactLd(first as list of forensics.Genotype, second as list of forensics.Genotype) {
    def rows as list of int init genotypeClasses(encode($first));
    def cols as list of int init genotypeClasses(encode($second));
    def n as int init len($rows);
    def observed as float init tableLogP($rows, $cols);
    def perm as list of int init [];
    for (def i in 0..$n) {
        $perm[] = $i;
    }
    def total as int init 0;
    def hits as int init 0;
    def more as bool init true;
    while ($more) {
        def shuffled as list of int init [];
        for (def i in 0..$n) {
            $shuffled[] = $cols[$perm[$i]];
        }
        $total = $total + 1;
        if (tableLogP($rows, $shuffled) <= $observed + 0.0000001) {
            $hits = $hits + 1;
        }
        # next lexicographic permutation
        def i as int init $n - 2;
        while ($i >= 0 and $perm[$i] >= $perm[$i + 1]) {
            $i = $i - 1;
        }
        if ($i < 0) {
            $more = false;
        } else {
            def j as int init $n - 1;
            while ($perm[$j] <= $perm[$i]) {
                $j = $j - 1;
            }
            def t as int init $perm[$i];
            $perm[$i] = $perm[$j];
            $perm[$j] = $t;
            def lo as int init $i + 1;
            def hi as int init $n - 1;
            while ($lo < $hi) {
                $t = $perm[$lo];
                $perm[$lo] = $perm[$hi];
                $perm[$hi] = $t;
                $lo = $lo + 1;
                $hi = $hi - 1;
            }
        }
    }
    return convert.toFloat($hits) / convert.toFloat($total);
}

func tableLogP(rows as list of int, cols as list of int) {
    def counts as map of string to int init {};
    for (def i in 0..len($rows)) {
        def key as string init "{$rows[$i]}/{$cols[$i]}";
        if (maps.has($counts, $key)) {
            $counts[$key] = $counts[$key] + 1;
        } else {
            $counts[$key] = 1;
        }
    }
    def out as float init 0.0;
    for (def key in $counts) {
        $out = $out - math.lgamma(convert.toFloat($counts[$key] + 1));
    }
    return $out;
}

# Three people: rows {11, 11, 22}, columns {33, 33, 44}. The 3! shuffles give
# two tables: the 44 with the 22 (n_11 = 2, n_22 = 1; weight 1/(2! 1!) = 1/2,
# 2 of 6 shuffles, probability 1/3) or with an 11 (all cells 1; weight 1,
# probability 2/3). Observed diagonal: p = 1/3. Observed off-diagonal: the
# most probable table, so p = 1 exactly.
func locusPair(a as list of string, b as list of string) {
    def out as list of list of forensics.Genotype init [[], []];
    for (def i in 0..len($a)) {
        def x as list of string init strings.split($a[$i], "/");
        def y as list of string init strings.split($b[$i], "/");
        $out[0][] = gl("L1", $x[0], $x[1]);
        $out[1][] = gl("L2", $y[0], $y[1]);
    }
    return $out;
}

func testLdThreePeopleByHand() {
    def diag as list of list of forensics.Genotype init locusPair(
        ["1/1", "1/1", "2/2"],
        ["3/3", "3/3", "4/4"]);
    # SE sqrt((1/3)(2/3)/10000) = .0047; tolerance about 4 SE
    within(ldPermutationTest($diag[0], $diag[1], 10000, 2024), 1.0 / 3.0, 0.02);
    near(exactLd($diag[0], $diag[1]), 1.0 / 3.0);
    def off as list of list of forensics.Genotype init locusPair(
        ["1/1", "1/1", "2/2"],
        ["3/3", "4/4", "3/3"]);
    near(ldPermutationTest($off[0], $off[1], 500, 7), 1.0);
}

# Four people, rows (2, 2) and columns (2, 2), perfectly associated: n_11 is
# hypergeometric, P(k) = C(2,k) C(2,2-k) / C(4,2) = 1/6, 4/6, 1/6 for
# k = 0, 1, 2. Observed k = 2 (probability 1/6): p = 1/6 + 1/6 = 1/3.
func testLdFourPeopleByHand() {
    def t as list of list of forensics.Genotype init locusPair(
        ["1/2", "1/2", "1/1", "1/1"],
        ["5/6", "5/6", "6/6", "6/6"]);
    near(exactLd($t[0], $t[1]), 1.0 / 3.0);
    within(ldPermutationTest($t[0], $t[1], 10000, 99), 1.0 / 3.0, 0.02);
}

func testLdEstimatesTheExactEnumeration() {
    # Seven people, three genotypes at each locus: 5040 shuffles enumerated.
    def t as list of list of forensics.Genotype init locusPair(
        [
            "1/1",
            "1/1",
            "1/2",
            "1/2",
            "2/2",
            "2/2",
            "1/2"
        ],
        ["7/7", "7/7", "7/8", "7/8", "8/8", "7/8", "8/8"]);
    def exact as float init exactLd($t[0], $t[1]);
    # SE at most sqrt(.25 / 20000) = .0035
    within(ldPermutationTest($t[0], $t[1], 20000, 31337), $exact, 0.015);
}

func testLdMonomorphicAndReproducible() {
    def t as list of list of forensics.Genotype init locusPair(
        ["1/1", "1/2", "2/2"],
        ["3/3", "3/3", "3/3"]);
    near(ldPermutationTest($t[0], $t[1], 100, 5), 1.0);
    def u as list of list of forensics.Genotype init locusPair(
        ["1/1", "1/2", "2/2", "1/1", "2/2"],
        ["3/3", "3/4", "4/4", "3/4", "4/4"]);
    near(ldPermutationTest($u[0], $u[1], 1000, 8), ldPermutationTest($u[0], $u[1], 1000, 8));
}

# Five profiles at three loci; Z is untyped in the last. Pairs in name order:
# (X, Y) over 5 people, (X, Z) and (Y, Z) over 4. The first pair starts from
# the seed, so it equals ldPermutationTest with that seed; the others follow
# on the same stream and are checked against exact enumeration.
func ldProfiles() {
    return [
        forensics.profile("1", [gl("Y", "3", "3"), gl("X", "1", "1"), gl("Z", "5", "5")]),
        forensics.profile("2", [gl("Y", "3", "3"), gl("X", "1", "1"), gl("Z", "5", "5")]),
        forensics.profile("3", [gl("Y", "4", "4"), gl("X", "2", "2"), gl("Z", "6", "6")]),
        forensics.profile("4", [gl("Y", "4", "4"), gl("X", "2", "2"), gl("Z", "6", "6")]),
        forensics.profile("5", [gl("Y", "3", "4"), gl("X", "1", "2")])
    ];
}

func testLdPairwise() {
    def r as list of LdTest init ldPermutationTests(ldProfiles(), 10000, 4242);
    testing.assertEqual(len($r), 3);
    testing.assertEqual($r[0].first, "X");
    testing.assertEqual($r[0].second, "Y");
    testing.assertEqual($r[0].n, 5);
    testing.assertEqual($r[1].second, "Z");
    testing.assertEqual($r[1].n, 4);
    testing.assertEqual($r[2].first, "Y");
    def xy as list of list of forensics.Genotype init locusPair(
        ["1/1", "1/1", "2/2", "2/2", "1/2"],
        ["3/3", "3/3", "4/4", "4/4", "3/4"]);
    near($r[0].pValue, ldPermutationTest($xy[0], $xy[1], 10000, 4242));
    # (X, Z) over four people is testLdFourPeopleByHand's table up to labels,
    # perfectly associated: p = 1/3.
    within($r[1].pValue, 1.0 / 3.0, 0.02);
    within($r[2].pValue, 1.0 / 3.0, 0.02);
    # a pair typed together in fewer than two people has no table: p = 1
    def sparse as list of LdTest init ldPermutationTests(
        [
            forensics.profile("a", [gl("X", "1", "1")]),
            forensics.profile("b", [gl("Y", "1", "2")])
        ],
        100,
        1);
    testing.assertEqual($sparse[0].n, 0);
    near($sparse[0].pValue, 1.0);
}

func noPopulations() {
    fStatistics([]);
}

func emptyPopulation() {
    fStatistics([[gl("A", "a", "a"), gl("A", "a", "b")], []]);
}

func populationsMixLoci() {
    fStatistics([[gl("A", "a", "a"), gl("A", "a", "b")], [gl("B", "a", "b")]]);
}

func singletons() {
    fStatistics([[gl("A", "a", "a")], [gl("A", "a", "b")]]);
}

func monomorphicFst() {
    fst(fStatistics([[gl("A", "a", "a"), gl("A", "a", "a")], [gl("A", "a", "a")]]));
}

func tableNoLoci() {
    fStatisticsTable(profiles2(), []);
}

func tableUntypedPopulation() {
    fStatisticsTable(profiles2(), ["C"]);
}

func combineEmpty() {
    combineFStatistics([]);
}

func combineMismatch() {
    combineFStatistics([fixedApart(), fStatistics([[gl("S", "a", "a"), gl("S", "a", "b")]])]);
}

func ldLengths() {
    ldPermutationTest([gl("L1", "1", "1")], [gl("L2", "1", "1"), gl("L2", "1", "2")], 10, 1);
}

func ldSameLocus() {
    ldPermutationTest([gl("L1", "1", "1")], [gl("L1", "1", "2")], 10, 1);
}

func ldZeroPermutations() {
    ldPermutationTest([gl("L1", "1", "1")], [gl("L2", "1", "2")], 0, 1);
}

func ldBadSeed() {
    ldPermutationTests(ldProfiles(), 10, 0);
}

func ldNoProfiles() {
    ldPermutationTests([], 10, 1);
}

func ldOneLocus() {
    ldPermutationTests([forensics.profile("a", [gl("X", "1", "1")])], 10, 1);
}

func testPopulationGeneticsErrors() {
    testing.assertThrows("noPopulations", "popstats");
    testing.assertThrows("emptyPopulation", "popstats");
    testing.assertThrows("populationsMixLoci", "popstats");
    testing.assertThrows("singletons", "popstats");
    testing.assertThrows("monomorphicFst", "popstats");
    testing.assertThrows("tableNoLoci", "popstats");
    testing.assertThrows("tableUntypedPopulation", "popstats");
    testing.assertThrows("combineEmpty", "popstats");
    testing.assertThrows("combineMismatch", "popstats");
    testing.assertThrows("ldLengths", "popstats");
    testing.assertThrows("ldSameLocus", "popstats");
    testing.assertThrows("ldZeroPermutations", "popstats");
    testing.assertThrows("ldBadSeed", "popstats");
    testing.assertThrows("ldNoProfiles", "popstats");
    testing.assertThrows("ldOneLocus", "popstats");
}

# --- errors -----------------------------------------------------------------

func missingLocus() {
    expectedHeterozygosity(fixture(), "SE33");
}

func partialTable() {
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("p"),
        "L",
        100,
        {"10": 0.5, "11": 0.4});
    matchProbability($fdb, "L");
}

func emptyLocus() {
    def fdb as forensics.FrequencyDb init forensics.withLocus(forensics.db("p"), "L", 100, {});
    powerOfDiscrimination($fdb, "L");
}

func noSampleSize() {
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("p"),
        "L",
        0,
        {"10": 0.5, "11": 0.5});
    unbiasedHeterozygosity($fdb, "L");
}

func emptyDb() {
    statsTable(forensics.db("empty"));
}

func emptyCombined() {
    combinedStats([]);
}

func emptyPowers() {
    combinedPower([]);
}

func powerAboveOne() {
    combinedPower([1.5]);
}

func emptyLogPowers() {
    log10CombinedComplement([]);
}

func powerOfOneLog() {
    log10CombinedComplement([0.5, 1.0]);
}

func emptySample() {
    observedHeterozygosity([]);
}

func mixedLoci() {
    alleleCounts([gl("L", "10", "11"), gl("M", "1", "2")]);
}

func emptyAllele() {
    estimateFrequencies([forensics.Genotype{locus: "L", a: "", b: "10"}]);
}

func noHomozygote() {
    observedTypicalPaternityIndex([gl("L", "10", "11")]);
}

func exactMultiallelic() {
    hweExactTest(threeHoms());
}

func zeroPermutations() {
    hwePermutationTest(threeHoms(), 0, 1);
}

func tooManyPermutations() {
    hwePermutationTest(threeHoms(), MAX_PERMUTATIONS + 1, 1);
}

func zeroSeed() {
    hwePermutationTest(threeHoms(), 10, 0);
}

func seedAtModulus() {
    hwePermutationTest(threeHoms(), 10, 2147483647);
}

func testErrorsCarryThePopstatsKind() {
    testing.assertThrows("missingLocus", "popstats");
    testing.assertThrows("partialTable", "popstats");
    testing.assertThrows("emptyLocus", "popstats");
    testing.assertThrows("noSampleSize", "popstats");
    testing.assertThrows("emptyDb", "popstats");
    testing.assertThrows("emptyCombined", "popstats");
    testing.assertThrows("emptyPowers", "popstats");
    testing.assertThrows("powerAboveOne", "popstats");
    testing.assertThrows("emptyLogPowers", "popstats");
    testing.assertThrows("powerOfOneLog", "popstats");
    testing.assertThrows("emptySample", "popstats");
    testing.assertThrows("mixedLoci", "popstats");
    testing.assertThrows("emptyAllele", "popstats");
    testing.assertThrows("noHomozygote", "popstats");
    testing.assertThrows("exactMultiallelic", "popstats");
    testing.assertThrows("zeroPermutations", "popstats");
    testing.assertThrows("tooManyPermutations", "popstats");
    testing.assertThrows("zeroSeed", "popstats");
    testing.assertThrows("seedAtModulus", "popstats");
}
