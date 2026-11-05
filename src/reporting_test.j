# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# White-box test overlay for reporting.j. Run with `jennifer test reporting_test.j`.
#
# Locus "L": alleles 10/11/12/13 at .1/.2/.3/.4 with N = 100, so 2N = 200 and
# the counts are 20/40/60/80. Size-bias expectations are written out from
# FSR-GUI-0012 11.3: (x + 2)/(2N + 4) per heterozygote allele, (x + 4)/(2N + 4)
# for a homozygote. Posterior quantiles are checked against distributions with
# a closed-form or independently computed quantile (Beta(2, 2), the arcsine
# Beta(1/2, 1/2), and Beta(1/3, 8/3) inverted through math.regBetaI), each
# tolerance justified by the Monte Carlo error of a nearest-rank quantile,
# sd ≈ sqrt(q (1 - q) / M) / f(x_q). Verbal bands are the ENFSI (2015) table.
use testing;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

func within(actual as float, expected as float, tol as float) {
    testing.assertTrue(math.abs($actual - $expected) < $tol);
}

func gl(locus as string, a as string, b as string) {
    return forensics.genotype($locus, $a, $b);
}

func fixture() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("test"),
        "L",
        100,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    # M: two alleles at .5 with N = 1, counts 1 and 1.
    return forensics.withLocus($out, "M", 1, {"1": 0.5, "2": 0.5});
}

# --- size-bias correction ---------------------------------------------------

func testSizeBiasAlleleFrequencyHeterozygote() {
    # 10/11: 10 -> (20 + 2) / (200 + 4) = 22/204; 11 -> 42/204
    near(sizeBiasAlleleFrequency(fixture(), gl("L", "10", "11"), "10"), 22.0 / 204.0);
    near(sizeBiasAlleleFrequency(fixture(), gl("L", "10", "11"), "11"), 42.0 / 204.0);
    # an allele the two people do not carry only sees the larger denominator: 60/204
    near(sizeBiasAlleleFrequency(fixture(), gl("L", "10", "11"), "12"), 60.0 / 204.0);
}

func testSizeBiasAlleleFrequencyHomozygoteAndUnlisted() {
    # 12/12: (60 + 4) / 204
    near(sizeBiasAlleleFrequency(fixture(), gl("L", "12", "12"), "12"), 64.0 / 204.0);
    # 10/14, 14 unlisted (x = 0): (0 + 2) / 204, no 5/2N floor
    near(sizeBiasAlleleFrequency(fixture(), gl("L", "10", "14"), "14"), 2.0 / 204.0);
}

func testSizeBiasGenotypeFrequencyThetaZero() {
    # 2 * 22/204 * 42/204 = 1848 / 41616
    near(sizeBiasGenotypeFrequency(fixture(), gl("L", "10", "11"), 0.0), 1848.0 / 41616.0);
    # (64/204)^2 = 4096 / 41616
    near(sizeBiasGenotypeFrequency(fixture(), gl("L", "12", "12"), 0.0), 4096.0 / 41616.0);
}

func testSizeBiasGenotypeFrequencyWithThetaMatchesCore() {
    # The core's NRC II 4.2 on a table holding the corrected frequencies
    # (N = 10^6 puts its floor, 2.5e-6, below both) is a separate code path.
    def corrected as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("corrected"),
        "L",
        1000000,
        {"10": 22.0 / 204.0, "11": 42.0 / 204.0, "12": 64.0 / 204.0});
    for (def t in [0.0, 0.01, 0.03]) {
        near(
            sizeBiasGenotypeFrequency(fixture(), gl("L", "10", "11"), $t),
            forensics.genotypeFrequency($corrected, gl("L", "10", "11"), $t));
        near(
            sizeBiasGenotypeFrequency(fixture(), gl("L", "12", "12"), $t),
            forensics.genotypeFrequency($corrected, gl("L", "12", "12"), $t));
    }
    # Hand value at theta .01, 10/11:
    # 2 [.01 + .99 * 22/204] [.01 + .99 * 42/204] / (1.01 * 1.02)
    def u as float init 0.01 + 0.99 * 22.0 / 204.0;
    def v as float init 0.01 + 0.99 * 42.0 / 204.0;
    near(
        sizeBiasGenotypeFrequency(fixture(), gl("L", "10", "11"), 0.01),
        2.0 * $u * $v / (1.01 * 1.02));
}

func testSizeBiasProfile() {
    # L 10/11 and M 1/1. M: 2N = 2, (1 + 4) / (2 + 4) = 5/6, so (5/6)^2 = 25/36.
    # RMP = 1848/41616 * 25/36
    def p as forensics.Profile init forensics.profile(
        "S",
        [gl("L", "10", "11"), gl("M", "1", "1")]);
    def rmp as float init 1848.0 / 41616.0 * 25.0 / 36.0;
    near(sizeBiasRandomMatchProbability(fixture(), $p, 0.0), $rmp);
    near(log10SizeBiasRandomMatchProbability(fixture(), $p, 0.0), math.log10($rmp));
    near(sizeBiasLikelihoodRatio(fixture(), $p, 0.0), 1.0 / $rmp);
}

func testSizeBiasStaysAProbability() {
    # A table summing to 1.05 (accepted by withLocus), N = 1: counts 1.2 and .9.
    # 1/1: (1.2 + 4) / 6 = .86667, squared .75111. 1/2: 2 * 3.2/6 * 2.9/6 = .51556.
    def big as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("b"),
        "X",
        1,
        {"1": 0.6, "2": 0.45});
    near(sizeBiasGenotypeFrequency($big, gl("X", "1", "1"), 0.0), 5.2 * 5.2 / 36.0);
    near(sizeBiasGenotypeFrequency($big, gl("X", "1", "2"), 0.0), 2.0 * 3.2 * 2.9 / 36.0);
}

func testSizeBiasVanishesForLargeN() {
    # N = 10^8: 10 -> (2e7 + 2)/(2e8 + 4), 11 -> (4e7 + 2)/(2e8 + 4).
    # Relative shift: (2 - .4)/(2e8 * .1) = 8e-8 and (2 - .8)/(2e8 * .2) = 3e-8,
    # so the genotype moves by about 1.1e-7 relative; 1e-6 bounds it.
    def large as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("large"),
        "L",
        100000000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    def g as forensics.Genotype init gl("L", "10", "11");
    def plain as float init forensics.genotypeFrequency($large, $g, 0.0);
    def corrected as float init sizeBiasGenotypeFrequency($large, $g, 0.0);
    testing.assertTrue($corrected > $plain);
    testing.assertTrue(math.abs($corrected / $plain - 1.0) < 0.000001);
}

func testSizeBiasRaisesTheMatchProbability() {
    # A heterozygote allele with p < 1/2 rises: (x + 2)/(2N + 4) > x/2N iff p < 1/2.
    def g as forensics.Genotype init gl("L", "10", "11");
    testing.assertTrue(
        sizeBiasGenotypeFrequency(fixture(), $g, 0.0) >
            forensics.genotypeFrequency(fixture(), $g, 0.0));
}

# --- generator and samplers -------------------------------------------------

func testGeneratorCheckValue() {
    # ISO C++ requires std::minstd_rand (48271, 2^31 - 1) to reach 399268537
    # on its 10000th call from seed 1.
    def state as int init 1;
    for (def i in 0..10000) {
        $state = nextState($state);
    }
    testing.assertEqual($state, 399268537);
}

func testLnGammasThreadsTheStream() {
    # One call drawing [3, .5, 3, .5] equals four one-deviate calls passing
    # the Stream on: generator state and the polar method's spare normal
    # carry over between calls, so splitting a draw changes no deviate.
    def whole as Stream init lnGammas(stream(4242), [3.0, 0.5, 3.0, 0.5]);
    def st as Stream init stream(4242);
    for (def j in 0..4) {
        def a as float init 3.0;
        if ($j % 2 == 1) {
            $a = 0.5;
        }
        $st = lnGammas($st, [$a]);
        testing.assertEqual(len($st.values), 1);
        testing.assertTrue($st.values[0] == $whole.values[$j]);
    }
    testing.assertEqual($st.state, $whole.state);
    testing.assertTrue($st.hasSpare == $whole.hasSpare);
}

# gammaMoments draws m Gamma(alpha, 1) deviates, one lnGammas call each, and
# returns their mean and variance.
func gammaMoments(alpha as float, m as int, seed as int) {
    def st as Stream init stream($seed);
    def s as float init 0.0;
    def s2 as float init 0.0;
    for (def i in 0..$m) {
        $st = lnGammas($st, [$alpha]);
        def x as float init math.exp($st.values[0]);
        $s = $s + $x;
        $s2 = $s2 + $x * $x;
    }
    def mean as float init $s / convert.toFloat($m);
    return [$mean, $s2 / convert.toFloat($m) - $mean * $mean];
}

func testGammaMoments() {
    # Gamma(a, 1) has mean a and variance a. M = 20000: SE of the mean
    # sqrt(a/M) is .0122 at a = 3, .0050 at a = .5 (the U^(1/a) branch),
    # .0016 at a = .05; tolerances are 4 SE. The variance's SE is
    # sqrt((mu4 - a^2)/M), mu4 = 3a^2 + 6a: .058 at a = 3, .0122 at a = .5;
    # tolerances 4 SE. Both depend on the normal deviates' first two moments.
    def m as list of float init gammaMoments(3.0, 20000, 7);
    within($m[0], 3.0, 0.05);
    within($m[1], 3.0, 0.24);
    $m = gammaMoments(0.5, 20000, 11);
    within($m[0], 0.5, 0.02);
    within($m[1], 0.5, 0.05);
    within(gammaMoments(0.05, 20000, 13)[0], 0.05, 0.0065);
}

func testArcsineQuantileFromGammaPair() {
    # X/(X + Y), X, Y ~ Gamma(1/2): Beta(1/2, 1/2), the arcsine law, CDF
    # (2/pi) asin(sqrt(x)); its .9 quantile is sin^2(.45 pi) ≈ .97553, density
    # there 1/(pi sqrt(x(1-x))) ≈ 2.060. M = 20000: sd sqrt(.09/M)/2.06 ≈
    # .00103; tolerance .005 (about 5 sd). The median is 1/2: sd
    # sqrt(.25/M) * pi/2 ≈ .0056 (density 2/pi), tolerance .025.
    def st as Stream init stream(99);
    def xs as list of float init [];
    for (def i in 0..20000) {
        $st = lnGammas($st, [0.5, 0.5]);
        def a as float init math.exp($st.values[0]);
        def b as float init math.exp($st.values[1]);
        $xs[] = $a / ($a + $b);
    }
    def sorted as list of float init lists.sort($xs);
    def s as float init math.sin(0.45 * math.PI);
    within(quantileOf($sorted, 0.9), $s * $s, 0.005);
    within(quantileOf($sorted, 0.5), 0.5, 0.025);
}

func testQuantileOfNearestRank() {
    # ceil(q M)-th smallest of [1..4]: q .25 -> 1st, .26 -> 2nd, .5 -> 2nd, .99 -> 4th
    def xs as list of float init [1.0, 2.0, 3.0, 4.0];
    near(quantileOf($xs, 0.25), 1.0);
    near(quantileOf($xs, 0.26), 2.0);
    near(quantileOf($xs, 0.5), 2.0);
    near(quantileOf($xs, 0.99), 4.0);
    near(quantileOf($xs, 0.000001), 1.0);
}

func testQuantileOfExactProductKeepsItsRank() {
    # [1..100]: q M is a whole number, so the rank is q M itself. In floating
    # point .07 * 100 = 7.000000000000001, .14 * 100 = 14.000000000000002,
    # .28 and .56 likewise; a bare ceil gave ranks 8, 15, 29 and 57.
    def xs as list of float init [];
    for (def i in 1..101) {
        $xs[] = convert.toFloat($i);
    }
    near(quantileOf($xs, 0.07), 7.0);
    near(quantileOf($xs, 0.14), 14.0);
    near(quantileOf($xs, 0.28), 28.0);
    near(quantileOf($xs, 0.56), 56.0);
    # a true fraction still rounds up: .071 * 100 = 7.1 -> 8th
    near(quantileOf($xs, 0.071), 8.0);
}

# --- Bayesian point estimate -----------------------------------------------
#
# FSR-GUI-0012 11.4 on locus L: k = 4 listed alleles, 2N = 200, counts
# 20/40/60/80. Partial table P: N = 50 (2N = 100), alleles 1/2 at .5/.3 (counts
# 50/30, 20 unlisted), k = 2.

func partial() {
    return forensics.withLocus(fixture(), "P", 50, {"1": 0.5, "2": 0.3});
}

func testBayesianFrequencyByHand() {
    # uniform: (x + 1) / (2N + k) = 21/204; unlisted 14: 1/204
    near(bayesianFrequency(fixture(), "L", "10", PRIOR_UNIFORM), 21.0 / 204.0);
    near(bayesianFrequency(fixture(), "L", "14", PRIOR_UNIFORM), 1.0 / 204.0);
    # 11.4.1: (x + 1/k) / (2N + 1) = 40.25/201; unlisted 14 (x = 0, k still 4): .25/201
    near(bayesianFrequency(fixture(), "L", "11", PRIOR_ONE_OVER_K), 40.25 / 201.0);
    near(bayesianFrequency(fixture(), "L", "14", PRIOR_ONE_OVER_K), 0.25 / 201.0);
    # 11.4.2: (x + 1/(k + 2)) / (2N + 1) = (60 + 1/6)/201; unlisted (1/6)/201
    def sixth as float init 1.0 / 6.0;
    near(
        bayesianFrequency(fixture(), "L", "12", PRIOR_ONE_OVER_K_PLUS_TWO),
        (60.0 + $sixth) / 201.0);
    near(bayesianFrequency(fixture(), "L", "14", PRIOR_ONE_OVER_K_PLUS_TWO), $sixth / 201.0);
}

func testBayesianFrequencyPartialTable() {
    # P: 1 -> uniform 51/102, 1/k (50 + .5)/101, 1/(k + 2) (50 + .25)/101
    near(bayesianFrequency(partial(), "P", "1", PRIOR_UNIFORM), 51.0 / 102.0);
    near(bayesianFrequency(partial(), "P", "1", PRIOR_ONE_OVER_K), 50.5 / 101.0);
    near(bayesianFrequency(partial(), "P", "1", PRIOR_ONE_OVER_K_PLUS_TWO), 50.25 / 101.0);
    # unlisted 3: .25/101 under 1/(k + 2)
    near(bayesianFrequency(partial(), "P", "3", PRIOR_ONE_OVER_K_PLUS_TWO), 0.25 / 101.0);
}

func testBayesianRandomMatchProbabilityByHand() {
    # theta 0, 1/k: L 10/11 2 (20.25/201)(40.25/201); P 1/1 (50.5/101)^2
    def p as forensics.Profile init forensics.profile(
        "S",
        [gl("L", "10", "11"), gl("P", "1", "1")]);
    def a as float init 2.0 * (20.25 / 201.0) * (40.25 / 201.0);
    def b as float init (50.5 / 101.0) * (50.5 / 101.0);
    near(bayesianRandomMatchProbability(partial(), $p, 0.0, PRIOR_ONE_OVER_K), $a * $b);
    near(
        log10BayesianRandomMatchProbability(partial(), $p, 0.0, PRIOR_ONE_OVER_K),
        math.log10($a) + math.log10($b));
    near(bayesianLikelihoodRatio(partial(), $p, 0.0, PRIOR_ONE_OVER_K), 1.0 / ($a * $b));
}

func testBayesianThetaIsTheCoresFormula() {
    # theta .02: NRC II 4.2 at the estimates equals the core's genotypeFrequency
    # on a table holding those estimates with n = 0 (no floor).
    def sixth as float init 1.0 / 6.0;
    def pa as float init (20.0 + $sixth) / 201.0;
    def pb as float init $sixth / 201.0;
    def est as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("est"),
        "L",
        0,
        {"10": $pa, "14": $pb});
    def het as forensics.Profile init forensics.profile("S", [gl("L", "10", "14")]);
    near(
        bayesianRandomMatchProbability(fixture(), $het, 0.02, PRIOR_ONE_OVER_K_PLUS_TWO),
        forensics.genotypeFrequency($est, gl("L", "10", "14"), 0.02));
    def hom as forensics.Profile init forensics.profile("H", [gl("L", "10", "10")]);
    near(
        bayesianRandomMatchProbability(fixture(), $hom, 0.02, PRIOR_ONE_OVER_K_PLUS_TWO),
        forensics.genotypeFrequency($est, gl("L", "10", "10"), 0.02));
}

# posteriorMeans lists the Dirichlet posterior means of a genotype's allele(s)
# and the total parameter mass.
func posteriorMeans(fdb as forensics.FrequencyDb, g as forensics.Genotype, prior as int) {
    def lp as LocusPrior init locusPrior($fdb, $g, $prior);
    def total as float init $lp.alphaA + $lp.alphaRest;
    if (not $lp.homozygous) {
        $total = $total + $lp.alphaB;
    }
    return [$lp.alphaA / $total, $lp.alphaB / $total, $total];
}

func testPosteriorOneOverKPlusTwoMean() {
    # 1/(k + 2): k + 2 = 6 categories, alpha 1/6, prior mass 1, so the
    # parameters sum to 2N + 1 = 201 and each mean is (x + 1/6)/201, whether
    # the profile alleles are listed (10/11), one is not (10/14), or a
    # homozygote is not (14/14). Partial P (1/2): 2N + 1 = 101, (x + .25)/101;
    # the 20 unlisted copies and two slots' .5 are the rest class.
    def sixth as float init 1.0 / 6.0;
    def m as list of float init posteriorMeans(
        fixture(),
        gl("L", "10", "11"),
        PRIOR_ONE_OVER_K_PLUS_TWO);
    near($m[0], (20.0 + $sixth) / 201.0);
    near($m[1], (40.0 + $sixth) / 201.0);
    near($m[2], 201.0);
    $m = posteriorMeans(fixture(), gl("L", "10", "14"), PRIOR_ONE_OVER_K_PLUS_TWO);
    near($m[0], (20.0 + $sixth) / 201.0);
    near($m[1], $sixth / 201.0);
    near($m[2], 201.0);
    $m = posteriorMeans(fixture(), gl("L", "14", "14"), PRIOR_ONE_OVER_K_PLUS_TWO);
    near($m[0], $sixth / 201.0);
    near($m[2], 201.0);
    $m = posteriorMeans(partial(), gl("P", "1", "2"), PRIOR_ONE_OVER_K_PLUS_TWO);
    near($m[0], 50.25 / 101.0);
    near($m[1], 30.25 / 101.0);
    near($m[2], 101.0);
    # each mean equals the point estimate
    near($m[0], bayesianFrequency(partial(), "P", "1", PRIOR_ONE_OVER_K_PLUS_TWO));
}

func testPosteriorOneOverKUsesKPrime() {
    # 1/k counts the profile's unlisted alleles: all listed (10/11) k' = 4,
    # mean (x + 1/4)/201 = 11.4.1; 10/14 k' = 5, (x + 1/5)/201 (neither);
    # 14/15 k' = 6, (x + 1/6)/201 = 11.4.2.
    def m as list of float init posteriorMeans(fixture(), gl("L", "10", "11"), PRIOR_ONE_OVER_K);
    near($m[0], bayesianFrequency(fixture(), "L", "10", PRIOR_ONE_OVER_K));
    $m = posteriorMeans(fixture(), gl("L", "10", "14"), PRIOR_ONE_OVER_K);
    near($m[0], 20.2 / 201.0);
    near($m[1], 0.2 / 201.0);
    near($m[2], 201.0);
    $m = posteriorMeans(fixture(), gl("L", "14", "15"), PRIOR_ONE_OVER_K);
    near($m[0], bayesianFrequency(fixture(), "L", "14", PRIOR_ONE_OVER_K_PLUS_TWO));
    near($m[2], 201.0);
}

func testPosteriorScalesATableSummingAboveOne() {
    # S: N = 10, a .5 / b .55, listed copies 10 + 11 = 21 > 2N = 20; scaled
    # by 20/21 to 200/21 and 220/21. 1/(k + 2), k = 2: alpha 1/4, mass 1, so
    # the parameters sum to 2N + 1 = 21 (unscaled: 22). a/a: alpha_A =
    # 200/21 + 1/4 = 821/84, mean 821/84/21 = 821/1764.
    def s as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("s"),
        "S",
        10,
        {"a": 0.5, "b": 0.55});
    def m as list of float init posteriorMeans($s, gl("S", "a", "a"), PRIOR_ONE_OVER_K_PLUS_TWO);
    near($m[2], 21.0);
    near($m[0], 821.0 / 1764.0);
    # a/b, uniform prior (k' = 2, alpha 1, no rest): 20 + 2 = 22 in all,
    # a's mean (200/21 + 1)/22 = 221/462.
    $m = posteriorMeans($s, gl("S", "a", "b"), PRIOR_UNIFORM);
    near($m[2], 22.0);
    near($m[0], 221.0 / 462.0);
}

func testPosteriorQuantilesOneOverKPlusTwo() {
    # M (counts 1/1, k = 2), profile 1/1, alpha 1/4: p1 ~ Beta(1 + 1/4,
    # 1 + 3/4), RMP p1^2. Quantiles of Beta(1.25, 1.75) by inverting
    # math.regBetaI. Monte Carlo sd of 2 log10 x_q at M = 20000 (numerically,
    # outside Jennifer): .0059 at .5, .0032 at .9; tolerances .03 and .016.
    def q as list of float init log10RmpQuantiles(
        fixture(),
        homozygoteM(),
        0.0,
        PRIOR_ONE_OVER_K_PLUS_TWO,
        [0.5, 0.9],
        20000,
        606);
    within($q[0], 2.0 * math.log10(betaQuantile(1.25, 1.75, 0.5)), 0.03);
    within($q[1], 2.0 * math.log10(betaQuantile(1.25, 1.75, 0.9)), 0.016);
    def lr as float init likelihoodRatioLowerBound(
        fixture(),
        homozygoteM(),
        0.0,
        PRIOR_ONE_OVER_K_PLUS_TWO,
        0.9,
        20000,
        606);
    near($lr, math.pow(10.0, -$q[1]));
}

# --- posterior quantiles ----------------------------------------------------

func homozygoteM() {
    return forensics.profile("H", [gl("M", "1", "1")]);
}

# beta22Quantile solves 3x^2 - 2x^3 = q (the Beta(2, 2) CDF) by bisection.
func beta22Quantile(q as float) {
    def lo as float init 0.0;
    def hi as float init 1.0;
    for (def i in 0..60) {
        def mid as float init ($lo + $hi) / 2.0;
        if (3.0 * $mid * $mid - 2.0 * $mid * $mid * $mid < $q) {
            $lo = $mid;
        } else {
            $hi = $mid;
        }
    }
    return $lo;
}

func testPosteriorQuantilesBeta22() {
    # M: counts 1/1, uniform prior: p1 ~ Beta(1 + 1, 1 + 1) = Beta(2, 2), and
    # 1/1's RMP is p1^2, so its log10 quantile is 2 log10 of p1's.
    # Median: p1 = 1/2 by symmetry, 2 log10 .5 = -.60206.
    # .9: 3x^2 - 2x^3 = .9 at x ≈ .8042, density 6x(1-x) ≈ .945.
    # M = 20000: sd of the p1 quantile sqrt(.09/M)/.945 ≈ .00225, times
    # d(2 log10 x)/dx = 2/(x ln 10) ≈ 1.08 -> .0024; the median's
    # sqrt(.25/M)/1.5 * 2/(.5 ln 10) ≈ .0041. Tolerances .012 and .02 (~5 sd).
    def q as list of float init log10RmpQuantiles(
        fixture(),
        homozygoteM(),
        0.0,
        PRIOR_UNIFORM,
        [0.5, 0.9],
        20000,
        2024);
    within($q[0], 2.0 * math.log10(0.5), 0.02);
    within($q[1], 2.0 * math.log10(beta22Quantile(0.9)), 0.012);
}

# betaQuantile inverts math.regBetaI (the interpreter's regularised incomplete
# beta, independent of the sampler) by bisection.
func betaQuantile(a as float, b as float, q as float) {
    def lo as float init 0.0;
    def hi as float init 1.0;
    for (def i in 0..60) {
        def mid as float init ($lo + $hi) / 2.0;
        if (math.regBetaI($mid, $a, $b) < $q) {
            $lo = $mid;
        } else {
            $hi = $mid;
        }
    }
    return $lo;
}

func testPosteriorQuantileOneOverKUnlistedAllele() {
    # M with profile 3/3, allele 3 unlisted: k = 3, alpha = 1/3 each.
    # p3 ~ Beta(0 + 1/3, (1 + 1/3) + (1 + 1/3)) = Beta(1/3, 8/3); RMP = p3^2.
    # Quartiles (numerically, outside Jennifer) x ≈ .00477 and .1554. The sd
    # of log10 x_q^2 is 2 sqrt(q(1-q)/M) / (f(x_q) x_q ln 10): .032 at .25 and
    # .013 at .75 for M = 20000 (the lower tail is steep in log). Tolerances
    # .15 and .06, about 4.5 sd.
    def p as forensics.Profile init forensics.profile("U", [gl("M", "3", "3")]);
    def q as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.0,
        PRIOR_ONE_OVER_K,
        [0.25, 0.75],
        20000,
        77);
    within($q[0], 2.0 * math.log10(betaQuantile(1.0 / 3.0, 8.0 / 3.0, 0.25)), 0.15);
    within($q[1], 2.0 * math.log10(betaQuantile(1.0 / 3.0, 8.0 / 3.0, 0.75)), 0.06);
}

func largeDb() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("large"),
        "L",
        1000000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    return forensics.withLocus($out, "K", 1000000, {"12": 0.3, "13": 0.7});
}

func testMedianNearPointEstimateForLargeN() {
    # N = 10^6. Plug-in RMP of L 10/11, K 12/12: 2 * .1 * .2 * .3^2 = .0036,
    # log10 -2.44370. Posterior sd of ln p is about sqrt((1-p)/(2N p)):
    # .00212 (p .1), .00141 (.2), and ln p^2 at .3: 2 * .00108 = .00216;
    # together ≈ .0033 in ln, .00145 in log10. The median's Monte Carlo sd at
    # M = 2000 is 1.2533 * .00145 / sqrt(2000) ≈ .00004; the prior and the
    # mean-median gap shift it by O(1/N) and O(sd^2) ≈ 1e-5. Tolerance .0003.
    def p as forensics.Profile init forensics.profile(
        "S",
        [gl("L", "10", "11"), gl("K", "12", "12")]);
    def q as list of float init log10RmpQuantiles(
        largeDb(),
        $p,
        0.0,
        PRIOR_UNIFORM,
        [0.5],
        2000,
        5);
    within($q[0], math.log10(0.0036), 0.0003);
    # With theta .01 the reference is the core's NRC II 4.2 plug-in value;
    # theta damps the sampling spread, so the same tolerance holds.
    def t as list of float init log10RmpQuantiles(
        largeDb(),
        $p,
        0.01,
        PRIOR_ONE_OVER_K,
        [0.5],
        2000,
        5);
    within($t[0], forensics.log10RandomMatchProbability(largeDb(), $p, 0.01), 0.0003);
}

func testIntervalShrinksWithN() {
    # The .01-.99 spread at N = 100 must exceed the one at N = 10^6.
    def p as forensics.Profile init forensics.profile("S", [gl("L", "10", "11")]);
    def small as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.0,
        PRIOR_UNIFORM,
        [0.01, 0.99],
        2000,
        3);
    def big as list of float init log10RmpQuantiles(
        largeDb(),
        $p,
        0.0,
        PRIOR_UNIFORM,
        [0.01, 0.99],
        2000,
        3);
    testing.assertTrue($small[1] -$small[0] > 10.0 * ($big[1] -$big[0]));
    testing.assertTrue($small[0] < $small[1]);
}

func testSeedFixesTheAnswer() {
    def p as forensics.Profile init forensics.profile("S", [gl("L", "10", "11")]);
    def a as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.01,
        PRIOR_UNIFORM,
        [0.99],
        500,
        42);
    def b as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.01,
        PRIOR_UNIFORM,
        [0.99],
        500,
        42);
    def c as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.01,
        PRIOR_UNIFORM,
        [0.99],
        500,
        43);
    testing.assertEqual($a[0], $b[0]);
    testing.assertTrue($a[0] != $c[0]);
}

func testLikelihoodRatioLowerBound() {
    # Same simulation as log10RmpQuantiles at that level: bound = 10^-q.
    def p as forensics.Profile init forensics.profile("S", [gl("L", "10", "11")]);
    def q as list of float init log10RmpQuantiles(
        fixture(),
        $p,
        0.01,
        PRIOR_UNIFORM,
        [0.99],
        1000,
        8);
    def lr as float init likelihoodRatioLowerBound(
        fixture(),
        $p,
        0.01,
        PRIOR_UNIFORM,
        0.99,
        1000,
        8);
    near($lr, math.pow(10.0, -$q[0]));
    # At N = 100 the posterior sd of ln(2 p10 p11) is about
    # sqrt(.9/20 + .8/40) ≈ .25, so the 99% bound sits about 2.3 sd (a factor
    # near 1.8) below the core's plug-in LR at theta .01.
    testing.assertTrue($lr < forensics.matchLikelihoodRatio(fixture(), $p, 0.01));
}

# --- verbal equivalents -----------------------------------------------------

func band(lr as float) {
    return verbalEquivalent($lr).label;
}

func testVerbalBandsFirstProposition() {
    # ENFSI (2015) p. 17, upper bound inclusive.
    testing.assertEqual(band(1.5), "weak support");
    testing.assertEqual(band(2.0), "weak support");
    testing.assertEqual(band(10.0), "weak support");
    testing.assertEqual(band(10.000001), "moderate support");
    testing.assertEqual(band(100.0), "moderate support");
    testing.assertEqual(band(100.00001), "moderately strong support");
    testing.assertEqual(band(1000.0), "moderately strong support");
    testing.assertEqual(band(1000.0001), "strong support");
    testing.assertEqual(band(10000.0), "strong support");
    testing.assertEqual(band(10000.001), "very strong support");
    testing.assertEqual(band(1000000.0), "very strong support");
    testing.assertEqual(band(1000000.1), "extremely strong support");
    testing.assertEqual(band(1e300), "extremely strong support");
    testing.assertEqual(verbalEquivalent(50.0).supports, SUPPORTS_FIRST);
}

func testVerbalBandsAlternative() {
    # 1/LR on the same scale; the mirror of (1, 10] is [.1, 1).
    testing.assertEqual(band(0.5), "weak support");
    testing.assertEqual(band(0.1), "weak support");
    testing.assertEqual(band(0.0999999), "moderate support");
    testing.assertEqual(band(0.01), "moderate support");
    testing.assertEqual(band(0.001), "moderately strong support");
    testing.assertEqual(band(0.0001), "strong support");
    testing.assertEqual(band(0.000001), "very strong support");
    testing.assertEqual(band(0.00000099), "extremely strong support");
    testing.assertEqual(verbalEquivalent(0.5).supports, SUPPORTS_ALTERNATIVE);
}

func testVerbalNeutral() {
    def v as VerbalEquivalent init verbalEquivalent(1.0);
    testing.assertEqual($v.label, "no support");
    testing.assertEqual($v.supports, SUPPORTS_NEITHER);
    testing.assertEqual(
        $v.statement,
        "The forensic findings do not support one proposition over the other.");
}

func testVerbalStatementsAndBounds() {
    def w as VerbalEquivalent init verbalEquivalent(5.0);
    testing.assertEqual(
        $w.statement,
        "The forensic findings provide weak support for the first proposition relative to the alternative.");
    near($w.lower, 1.0);
    near($w.upper, 10.0);
    testing.assertFalse($w.lowerInclusive);
    testing.assertTrue($w.upperInclusive);
    testing.assertTrue($w.hasUpper);
    def s as VerbalEquivalent init verbalEquivalent(500.0);
    testing.assertEqual(
        $s.statement,
        "The forensic findings provide moderately strong support for the first proposition rather than the alternative.");
    near($s.lower, 100.0);
    near($s.upper, 1000.0);
    def v as VerbalEquivalent init verbalEquivalent(50000.0);
    near($v.lower, 10000.0);
    near($v.upper, 1000000.0);
    def top as VerbalEquivalent init verbalEquivalent(1e9);
    near($top.lower, 1000000.0);
    testing.assertFalse($top.hasUpper);
    def alt as VerbalEquivalent init verbalEquivalent(0.05);
    testing.assertEqual(
        $alt.statement,
        "The forensic findings provide moderate support for the alternative rather than the first proposition.");
    near($alt.lower, 0.01);
    near($alt.upper, 0.1);
    testing.assertTrue($alt.lowerInclusive);
    testing.assertFalse($alt.upperInclusive);
    def bottom as VerbalEquivalent init verbalEquivalent(1e-9);
    near($bottom.lower, 0.0);
    near($bottom.upper, 0.000001);
    testing.assertFalse($bottom.lowerInclusive);
}

func testVerbalLog10() {
    testing.assertEqual(verbalEquivalentLog10(0.0).supports, SUPPORTS_NEITHER);
    testing.assertEqual(verbalEquivalentLog10(1.0).label, "weak support");
    testing.assertEqual(verbalEquivalentLog10(1.0000001).label, "moderate support");
    testing.assertEqual(verbalEquivalentLog10(6.0).label, "very strong support");
    testing.assertEqual(verbalEquivalentLog10(6.0000001).label, "extremely strong support");
    # Past float range: log10 LR = 400.
    testing.assertEqual(verbalEquivalentLog10(400.0).label, "extremely strong support");
    testing.assertEqual(verbalEquivalentLog10(-1.0).label, "weak support");
    testing.assertEqual(verbalEquivalentLog10(-1.0).supports, SUPPORTS_ALTERNATIVE);
    testing.assertEqual(verbalEquivalentLog10(-4.5).label, "very strong support");
    # Agrees with verbalEquivalent away from the edges.
    for (def x in [0.3, 3.0, 30.0, 300.0, 3000.0, 30000.0, 3e7, 0.03, 3e-5, 3e-8]) {
        testing.assertEqual(
            verbalEquivalentLog10(math.log10($x)).label,
            verbalEquivalent($x).label);
    }
}

# --- errors -----------------------------------------------------------------

func zeroLr() {
    verbalEquivalent(0.0);
}

func negativeLr() {
    verbalEquivalent(-2.0);
}

func noSampleSize() {
    def d as forensics.FrequencyDb init forensics.withLocus(forensics.db("d"), "L", 0, {"10": 0.5});
    sizeBiasAlleleFrequency($d, gl("L", "10", "10"), "10");
}

func missingLocus() {
    sizeBiasGenotypeFrequency(fixture(), gl("Z", "1", "2"), 0.0);
}

func badThetaSizeBias() {
    sizeBiasGenotypeFrequency(fixture(), gl("L", "10", "11"), 1.0);
}

func emptyProfileSizeBias() {
    sizeBiasRandomMatchProbability(fixture(), forensics.profile("E", []), 0.0);
}

func emptyProfileLog10() {
    log10SizeBiasRandomMatchProbability(fixture(), forensics.profile("E", []), 0.0);
}

func emptyProfileLR() {
    sizeBiasLikelihoodRatio(fixture(), forensics.profile("E", []), 0.0);
}

func badPrior() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, 4, [0.5], 10, 1);
}

func bayesianBadPrior() {
    bayesianFrequency(fixture(), "L", "10", 0);
}

func bayesianNoSampleSize() {
    def db as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("t"),
        "L",
        0,
        {"10": 0.5, "11": 0.5});
    bayesianFrequency($db, "L", "10", PRIOR_ONE_OVER_K);
}

func bayesianMissingLocus() {
    bayesianFrequency(fixture(), "Z", "1", PRIOR_ONE_OVER_K);
}

func bayesianNoAlleles() {
    def db as forensics.FrequencyDb init forensics.withLocus(forensics.db("t"), "E", 10, {});
    bayesianFrequency($db, "E", "1", PRIOR_ONE_OVER_K);
}

func bayesianEmptyProfile() {
    bayesianRandomMatchProbability(fixture(), forensics.profile("E", []), 0.0, PRIOR_UNIFORM);
}

func bayesianEmptyProfileLog10() {
    log10BayesianRandomMatchProbability(fixture(), forensics.profile("E", []), 0.0, PRIOR_UNIFORM);
}

func bayesianBadTheta() {
    bayesianLikelihoodRatio(
        fixture(),
        forensics.profile("S", [gl("L", "10", "11")]),
        1.0,
        PRIOR_UNIFORM);
}

func bayesianBadPriorRmp() {
    bayesianRandomMatchProbability(
        fixture(),
        forensics.profile("S", [gl("L", "10", "11")]),
        0.0,
        7);
}

func bayesianSubnormalLr() {
    def c as SubnormalCase init subnormalProfile();
    bayesianLikelihoodRatio($c.db, $c.profile, 0.0, PRIOR_ONE_OVER_K);
}

func noQuantiles() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [], 10, 1);
}

func quantileOne() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [1.0], 10, 1);
}

func quantileZero() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [0.0], 10, 1);
}

func zeroDraws() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [0.5], 0, 1);
}

func tooManyDraws() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [0.5], MAX_DRAWS + 1, 1);
}

func zeroSeed() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [0.5], 10, 0);
}

func seedAtModulus() {
    log10RmpQuantiles(fixture(), homozygoteM(), 0.0, PRIOR_UNIFORM, [0.5], 10, 2147483647);
}

func badThetaPosterior() {
    likelihoodRatioLowerBound(fixture(), homozygoteM(), -0.1, PRIOR_UNIFORM, 0.99, 10, 1);
}

func posteriorMissingLocus() {
    likelihoodRatioLowerBound(
        fixture(),
        forensics.profile("Z", [gl("Z", "1", "2")]),
        0.0,
        PRIOR_UNIFORM,
        0.99,
        10,
        1);
}

# subnormalProfile has 41 loci typed X/Y, alleles the table (A, B at .5,
# N = 10000) lacks: each corrected frequency is (0 + 2) / 20004, each locus
# 2 (2/20004)^2 ≈ 2.0e-8, the product about 10^-315.7, below the normal range.
func subnormalProfile() {
    def db as forensics.FrequencyDb init forensics.db("t");
    def gs as list of forensics.Genotype init [];
    for (def i in 0..41) {
        $db = forensics.withLocus($db, "L{$i}", 10000, {"A": 0.5, "B": 0.5});
        $gs[] = gl("L{$i}", "X", "Y");
    }
    return SubnormalCase{db: $db, profile: forensics.profile("S", $gs)};
}

def struct SubnormalCase {
    db as forensics.FrequencyDb,
    profile as forensics.Profile
};

func subnormalLr() {
    def c as SubnormalCase init subnormalProfile();
    sizeBiasLikelihoodRatio($c.db, $c.profile, 0.0);
}

func testSizeBiasLikelihoodRatioPastFloatRange() {
    # 1 / 10^-315.7 overflows; the reporting error names the log10 variant,
    # which gives 41 log10(2 (2/20004)^2).
    testing.assertThrows("subnormalLr", "reporting");
    def c as SubnormalCase init subnormalProfile();
    def rmp as float init sizeBiasRandomMatchProbability($c.db, $c.profile, 0.0);
    testing.assertTrue($rmp > 0.0 and $rmp < 1e-307);
    def f as float init 2.0 / 20004.0;
    near(
        log10SizeBiasRandomMatchProbability($c.db, $c.profile, 0.0),
        41.0 * math.log10(2.0 * $f * $f));
}

func testErrorsCarryTheReportingKind() {
    testing.assertThrows("zeroLr", "reporting");
    testing.assertThrows("negativeLr", "reporting");
    testing.assertThrows("noSampleSize", "reporting");
    testing.assertThrows("missingLocus", "reporting");
    testing.assertThrows("badThetaSizeBias", "reporting");
    testing.assertThrows("emptyProfileSizeBias", "reporting");
    testing.assertThrows("emptyProfileLog10", "reporting");
    testing.assertThrows("emptyProfileLR", "reporting");
    testing.assertThrows("badPrior", "reporting");
    testing.assertThrows("noQuantiles", "reporting");
    testing.assertThrows("quantileOne", "reporting");
    testing.assertThrows("quantileZero", "reporting");
    testing.assertThrows("zeroDraws", "reporting");
    testing.assertThrows("tooManyDraws", "reporting");
    testing.assertThrows("zeroSeed", "reporting");
    testing.assertThrows("seedAtModulus", "reporting");
    testing.assertThrows("badThetaPosterior", "reporting");
    testing.assertThrows("posteriorMissingLocus", "reporting");
    testing.assertThrows("bayesianBadPrior", "reporting");
    testing.assertThrows("bayesianNoSampleSize", "reporting");
    testing.assertThrows("bayesianMissingLocus", "reporting");
    testing.assertThrows("bayesianNoAlleles", "reporting");
    testing.assertThrows("bayesianEmptyProfile", "reporting");
    testing.assertThrows("bayesianEmptyProfileLog10", "reporting");
    testing.assertThrows("bayesianBadTheta", "reporting");
    testing.assertThrows("bayesianBadPriorRmp", "reporting");
    testing.assertThrows("bayesianSubnormalLr", "reporting");
}
