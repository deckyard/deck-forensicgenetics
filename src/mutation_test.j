# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for mutation.j. Run with `jennifer test mutation_test.j`.
#
# Locus "L" is kinship_test.j's: alleles 10/11/12/13 at .1/.2/.3/.4, N = 10000,
# so the 5/2N floor (0.00025) never bites. Matrix entries are worked by hand
# from the FamiliasLocus definitions; parentage results from the joint
# genotype sums written out beside each test. With no mutation and no silent
# allele every function is compared with its kinship.j original for exact
# equality.
use testing;
use math;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

func fixture() {
    return forensics.withLocus(
        forensics.db("test"),
        "L",
        10000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

func twoLocus() {
    def out as forensics.FrequencyDb init fixture();
    return forensics.withLocus($out, "M", 10000, {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

func gt(a as string, b as string) {
    return forensics.genotype("L", $a, $b);
}

func allelesL(silent as float) {
    return locusAlleles(fixture(), "L", [], $silent);
}

func matrixL(model as MutationSpec) {
    return mutationMatrix($model, allelesL(0.0));
}

# onLocus builds a one-locus table and its LocusAlleles.
func onLocus(locus as string, freqs as map of string to float) {
    def db as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("t"),
        $locus,
        10000,
        $freqs);
    return locusAlleles($db, $locus, [], 0.0);
}

func at(m as MutationMatrix, source as string, target as string) {
    return transitionProbability($m, $source, $target);
}

func assertRowsSumToOne(m as MutationMatrix) {
    for (def row in $m.rows) {
        def sum as float init 0.0;
        for (def x in $row) {
            $sum = $sum + $x;
        }
        near($sum, 1.0);
    }
}

func assertSameMatrix(a as MutationMatrix, b as MutationMatrix) {
    testing.assertEqual($a.alleles, $b.alleles);
    for (def i in 0..len($a.rows)) {
        for (def j in 0..len($a.rows)) {
            near($a.rows[$i][$j], $b.rows[$i][$j]);
        }
    }
}

# --- matrices by hand -------------------------------------------------------

func testLocusAllelesScalesDetectableFrequencies() {
    def la as LocusAlleles init locusAlleles(fixture(), "L", ["10", "14"], 0.05);
    # table order, then the unlisted typed 14 at the floor 5/20000, then silent
    testing.assertEqual($la.alleles, ["10", "11", "12", "13", "14", SILENT_ALLELE]);
    near($la.frequencies[0], 0.95 * 0.1);
    near($la.frequencies[4], 0.95 * 0.00025);
    near($la.frequencies[5], 0.05);
    testing.assertEqual(allelesL(0.0).alleles, ["10", "11", "12", "13"]);
}

func testEqualMatrixByHand() {
    # rate .003 over n = 4: 1 - .003 on the diagonal, .003 / 3 elsewhere
    def m as MutationMatrix init matrixL(equalModel(0.003));
    near(at($m, "10", "10"), 0.997);
    near(at($m, "10", "13"), 0.001);
    near(at($m, "12", "11"), 0.001);
    assertRowsSumToOne($m);
}

func testProportionalMatrixByHand() {
    # sum p(1-p) = 1 - (.01 + .04 + .09 + .16) = .7; alpha = .007 / 1 / .7 = .01
    def m as MutationMatrix init matrixL(proportionalModel(0.007));
    near(at($m, "10", "10"), 1.0 - 0.01 + 0.01 * 0.1);
    near(at($m, "10", "11"), 0.01 * 0.2);
    near(at($m, "12", "10"), 0.01 * 0.1);
    near(at($m, "13", "13"), 1.0 - 0.01 + 0.01 * 0.4);
    assertRowsSumToOne($m);
    # expected rate sum p_i alpha (1 - p_i) = .01 * .7
    near(overallRate($m, allelesL(0.0)), 0.007);
    testing.assertTrue(isStationary($m, allelesL(0.0)));
}

func testStepwiseMatrixByHand() {
    # rate .004, range .5. From 10: weights .5, .25, .125 to 11, 12, 13,
    # total .875, so .004 * .5 / .875 = .016 / 7 and so on. From 11: 10 and 12
    # one step (.5 each), 13 two (.25), total 1.25.
    def m as MutationMatrix init matrixL(stepwiseModel(0.004, 0.5));
    near(at($m, "10", "10"), 0.996);
    near(at($m, "10", "11"), 0.016 / 7.0);
    near(at($m, "10", "12"), 0.008 / 7.0);
    near(at($m, "10", "13"), 0.004 / 7.0);
    near(at($m, "11", "10"), 0.004 * 0.5 / 1.25);
    near(at($m, "11", "13"), 0.004 * 0.25 / 1.25);
    assertRowsSumToOne($m);
}

func testExtendedStepwiseMatrixByHand() {
    # 10, 11, 12 are one class, 11.3 another. rate .004, range .5, rate2 .001.
    # From 10: diagonal 1 - .004 - .001; 11 and 12 weights .5, .25 scaled to
    # .004; 11.3 takes all of rate2. From 11.3, alone in its class: diagonal
    # 1 - .001, each of the three others .001 / 3.
    def la as LocusAlleles init onLocus("V", {"10": 0.25, "11": 0.25, "11.3": 0.25, "12": 0.25});
    def m as MutationMatrix init mutationMatrix(extendedStepwiseModel(0.004, 0.5, 0.001), $la);
    near(at($m, "10", "10"), 0.995);
    near(at($m, "10", "11"), 0.004 * 0.5 / 0.75);
    near(at($m, "10", "12"), 0.004 * 0.25 / 0.75);
    near(at($m, "10", "11.3"), 0.001);
    near(at($m, "11.3", "11.3"), 0.999);
    near(at($m, "11.3", "12"), 0.001 / 3.0);
    assertRowsSumToOne($m);
}

func testMicroClassFollowsFamilias() {
    testing.assertEqual(microClass(15.0), 0);
    testing.assertEqual(microClass(15.2), 2);
    testing.assertEqual(microClass(15.3), 3);
    # 15.7 is nearer 16: FamiliasLocus's round() gives class -3
    testing.assertEqual(microClass(15.7), -3);
}

func testSilentAlleleNeverMutates() {
    def la as LocusAlleles init allelesL(0.05);
    def m as MutationMatrix init mutationMatrix(equalModel(0.003), $la);
    # Equal spreads over the four detectable alleles only
    near(at($m, "10", "11"), 0.001);
    near(at($m, "10", SILENT_ALLELE), 0.0);
    near(at($m, SILENT_ALLELE, SILENT_ALLELE), 1.0);
    near(at($m, SILENT_ALLELE, "10"), 0.0);
    assertRowsSumToOne($m);
}

func testProportionalWithSilentKeepsTheExpectedRate() {
    # detectable sum .95: alpha = .007 / .95 / .7, expected rate over all
    # alleles (silent included) .95 * alpha * .7 = .007
    def la as LocusAlleles init allelesL(0.05);
    def m as MutationMatrix init mutationMatrix(proportionalModel(0.007), $la);
    near(at($m, "10", "11"), 0.007 / 0.95 / 0.7 * 0.2);
    near(overallRate($m, $la), 0.007);
    testing.assertTrue(isStationary($m, $la));
    assertRowsSumToOne($m);
}

func testPerLocusRatesOverride() {
    def model as MutationSpec init withLocusRate(equalModel(0.003), "L", 0.006);
    near(at(matrixL($model), "10", "11"), 0.002);
    def many as MutationSpec init withLocusRates(equalModel(0.003), {"L": 0.009, "M": 0.0});
    near(at(matrixL($many), "10", "11"), 0.003);
    near(at(mutationMatrix($many, locusAlleles(twoLocus(), "M", [], 0.0)), "10", "10"), 1.0);
}

# --- models reduce to each other --------------------------------------------

func testEveryModelAtRateZeroIsTheIdentity() {
    def id as list of list of float init identityRows(4);
    for (def model in [
        NO_MUTATION,
        equalModel(0.0),
        proportionalModel(0.0),
        stepwiseModel(0.0, 0.5),
        stationaryModel(equalModel(0.0))
    ]) {
        testing.assertEqual(matrixL($model).rows, $id);
    }
}

func testStepwiseWithRangeOneIsEqual() {
    # every weight range^k = 1, so each of the 3 others gets rate / 3
    assertSameMatrix(matrixL(stepwiseModel(0.003, 1.0)), matrixL(equalModel(0.003)));
}

func testExtendedStepwiseWithoutMicrovariantsIsStepwise() {
    testing.assertEqual(
        matrixL(extendedStepwiseModel(0.004, 0.3, 0.002)).rows,
        matrixL(stepwiseModel(0.004, 0.3)).rows);
}

func testProportionalOnUniformFrequenciesIsEqual() {
    # p = 1/4: alpha = rate / (4 * 1/4 * 3/4) = rate * 4/3; off-diagonal
    # alpha / 4 = rate / 3
    def la as LocusAlleles init onLocus("U", {"10": 0.25, "11": 0.25, "12": 0.25, "13": 0.25});
    assertSameMatrix(
        mutationMatrix(proportionalModel(0.006), $la),
        mutationMatrix(equalModel(0.006), $la));
    testing.assertTrue(isStationary(mutationMatrix(equalModel(0.006), $la), $la));
}

func testProportionalIsUnchangedByPM() {
    def la as LocusAlleles init allelesL(0.0);
    assertSameMatrix(
        mutationMatrix(stationaryModel(proportionalModel(0.007)), $la),
        mutationMatrix(proportionalModel(0.007), $la));
}

# --- stationarity -----------------------------------------------------------

func testPMOnTwoAllelesByHand() {
    # Equal, rate .1, p = (.25, .75): M = [[.9, .1], [.1, .9]] has stationary
    # v = (.5, .5); R = c = .1; d_i = R v_i / (c p_i) = (2, 2/3).
    # S = [[1 - .2, .2], [.1 * 2/3, 1 - .1 * 2/3]]; .25 * .2 = .75 * (.2/3).
    def la as LocusAlleles init onLocus("T", {"10": 0.25, "11": 0.75});
    def plain as MutationMatrix init mutationMatrix(equalModel(0.1), $la);
    testing.assertFalse(isStationary($plain, $la));
    def s as MutationMatrix init mutationMatrix(stationaryModel(equalModel(0.1)), $la);
    near(at($s, "10", "11"), 0.2);
    near(at($s, "10", "10"), 0.8);
    near(at($s, "11", "10"), 0.2 / 3.0);
    near(at($s, "11", "11"), 1.0 - 0.2 / 3.0);
    testing.assertTrue(isStationary($s, $la));
    near(overallRate($s, $la), overallRate($plain, $la));
    near(overallRate($s, $la), 0.1);
}

func testPMMakesStepwiseStationaryAtTheSameRate() {
    for (def silent in [0.0, 0.05]) {
        def la as LocusAlleles init allelesL($silent);
        def plain as MutationMatrix init mutationMatrix(stepwiseModel(0.004, 0.5), $la);
        testing.assertFalse(isStationary($plain, $la));
        def s as MutationMatrix init mutationMatrix(
            stationaryModel(stepwiseModel(0.004, 0.5)),
            $la);
        testing.assertTrue(isStationary($s, $la));
        assertRowsSumToOne($s);
        near(overallRate($s, $la), overallRate($plain, $la));
    }
}

func pmDoesNotExist() {
    # d_0 (1 - M_00) = 2 * .6 > 1: the scaled diagonal would be negative
    mutationMatrix(stationaryModel(equalModel(0.6)), onLocus("T", {"10": 0.25, "11": 0.75}));
}

func pmWithoutUniqueStationaryDistribution() {
    # rate2 0: 10.1 and 10 never exchange, so the chain has two closed classes
    mutationMatrix(
        stationaryModel(stepwiseModel(0.01, 0.5)),
        onLocus("W", {"10": 0.3, "11": 0.3, "10.1": 0.2, "11.1": 0.2}));
}

func testPMFailuresThrow() {
    testing.assertThrows("pmDoesNotExist", "mutation");
    testing.assertThrows("pmWithoutUniqueStationaryDistribution", "mutation");
}

# --- validation -------------------------------------------------------------

func badEqualRate() {
    equalModel(1.0);
}

func badProportionalRate() {
    proportionalModel(-0.1);
}

func badRange() {
    stepwiseModel(0.002, 0.0);
}

func badRate2() {
    extendedStepwiseModel(0.002, 0.5, 1.0);
}

func ratesOverOne() {
    extendedStepwiseModel(0.6, 0.5, 0.5);
}

func badLocusRate() {
    withLocusRate(equalModel(0.002), "L", 1.5);
}

func badLocusRates() {
    withLocusRates(equalModel(0.002), {"L": -0.5});
}

func locusRateOverOneWithRate2() {
    withLocusRate(extendedStepwiseModel(0.1, 0.5, 0.5), "L", 0.6);
}

func stepwiseOnLetters() {
    mutationMatrix(stepwiseModel(0.002, 0.5), onLocus("AMEL", {"X": 0.5, "Y": 0.5}));
}

func stepwiseOnTwoDecimals() {
    mutationMatrix(stepwiseModel(0.002, 0.5), onLocus("T", {"15": 0.5, "15.25": 0.5}));
}

func proportionalTooHigh() {
    # max rate (1 - sum p^2) / (1 - min p) = (1 - .66) / .9 < .5: the .8
    # allele's diagonal 1 - alpha (1 - .8) goes negative
    mutationMatrix(proportionalModel(0.5), onLocus("T", {"10": 0.8, "11": 0.1, "12": 0.1}));
}

func unknownModel() {
    def m as MutationSpec init NO_MUTATION;
    $m.model = "random";
    $m.rate = 0.1;
    matrixL($m);
}

func typedSilent() {
    locusAlleles(fixture(), "L", [SILENT_ALLELE], 0.1);
}

func silentOutOfRange() {
    locusAlleles(fixture(), "L", [], 1.0);
}

func uncoveredLocus() {
    locusAlleles(fixture(), "Q", [], 0.0);
}

func transitionOfUnknownAllele() {
    transitionProbability(matrixL(NO_MUTATION), "10", "99");
}

func rateOverOtherAlleles() {
    overallRate(matrixL(NO_MUTATION), allelesL(0.05));
}

func stationaryOverOtherAlleles() {
    isStationary(matrixL(NO_MUTATION), allelesL(0.05));
}

func testModelValidation() {
    for (def f in [
        "badEqualRate",
        "badProportionalRate",
        "badRange",
        "badRate2",
        "ratesOverOne",
        "badLocusRate",
        "badLocusRates",
        "locusRateOverOneWithRate2",
        "stepwiseOnLetters",
        "stepwiseOnTwoDecimals",
        "proportionalTooHigh",
        "unknownModel",
        "typedSilent",
        "silentOutOfRange",
        "uncoveredLocus",
        "transitionOfUnknownAllele",
        "rateOverOtherAlleles",
        "stationaryOverOtherAlleles"
    ]) {
        testing.assertThrows($f, "mutation");
    }
}

# --- trio -------------------------------------------------------------------

func testMendelianTrioIsExactlyKinship() {
    def trios as list of list of forensics.Genotype init [
        [gt("10", "11"), gt("11", "12"), gt("12", "13")],
        [gt("10", "10"), gt("10", "11"), gt("11", "11")],
        [gt("10", "11"), gt("10", "11"), gt("10", "11")],
        [gt("12", "13"), gt("12", "12"), gt("12", "13")],
        [gt("10", "11"), gt("11", "13"), gt("12", "12")]
    ];
    for (def t in $trios) {
        for (def theta in [0.0, 0.01, 0.1]) {
            def want as float init kinship.paternityIndexAtTheta(
                fixture(),
                $t[0],
                $t[1],
                $t[2],
                kinship.NO_MUTATION,
                $theta);
            testing.assertEqual(
                paternityIndexAt(fixture(), $t[0], $t[1], $t[2], MENDELIAN, $theta),
                $want);
            def zero as ParentageModel init parentageModel(
                stepwiseModel(0.0, 0.5),
                equalModel(0.0));
            testing.assertEqual(
                paternityIndexAt(fixture(), $t[0], $t[1], $t[2], $zero, $theta),
                $want);
        }
    }
}

func testTrioPaternalMutationByHand() {
    # mother 10/11, child 11/13, alleged father 12/12: the 13 needs a 12 -> 13
    # change. X = T_M(11) T_F(13) = .5 * M[12][13]; Y = T_M(11) p13 = .5 * .4.
    def m as forensics.Genotype init gt("10", "11");
    def c as forensics.Genotype init gt("11", "13");
    def f as forensics.Genotype init gt("12", "12");
    def equal as ParentageModel init parentageModel(equalModel(0.003), NO_MUTATION);
    # Equal: M[12][13] = .001
    near(paternityIndexAt(fixture(), $m, $c, $f, $equal, 0.0), 0.5 * 0.001 / (0.5 * 0.4));
    # Stepwise .004, range .5: from 12, 11 and 13 one step (.5), 10 two (.25);
    # M[12][13] = .004 * .5 / 1.25 = .0016
    def step as ParentageModel init parentageModel(stepwiseModel(0.004, 0.5), NO_MUTATION);
    near(paternityIndexAt(fixture(), $m, $c, $f, $step, 0.0), 0.5 * 0.0016 / (0.5 * 0.4));
}

func testTrioMaternalMutationByHand() {
    # mother 10/10, child 11/12, alleged father 12/12: the child's 11 must be
    # a maternal mutation. Maternal Equal .003: T_M(11) = T_M(12) = .001;
    # T_F(12) = 1, T_F(11) = 0.
    # X = T_M(11) T_F(12) + T_M(12) T_F(11) = .001
    # Y = T_M(11) p12 + T_M(12) p11 = .001 * .3 + .001 * .2 = .0005
    def mat as ParentageModel init parentageModel(NO_MUTATION, equalModel(0.003));
    near(
        paternityIndexAt(fixture(), gt("10", "10"), gt("11", "12"), gt("12", "12"), $mat, 0.0),
        2.0);
}

func testTrioBothParentsMutateByHand() {
    # mother 10/11, child 11/13, alleged father 12/12, Equal .003 for both.
    # T_M(11) = .5 (M[10][11] + M[11][11]) = .5 (.001 + .997) = .499
    # T_M(13) = .5 (.001 + .001) = .001; T_F(13) = T_F(11) = .001
    # X = .499 * .001 + .001 * .001 = .0005
    # Y = .499 * p13 + .001 * p11 = .1996 + .0002 = .1998
    def both as ParentageModel init parentageModel(equalModel(0.003), equalModel(0.003));
    near(
        paternityIndexAt(fixture(), gt("10", "11"), gt("11", "13"), gt("12", "12"), $both, 0.0),
        0.0005 / 0.1998);
}

func maternalExclusion() {
    paternityIndexAt(fixture(), gt("10", "10"), gt("11", "12"), gt("12", "12"), MENDELIAN, 0.0);
}

func trioAtMixedLoci() {
    paternityIndexAt(
        fixture(),
        gt("10", "11"),
        forensics.genotype("M", "10", "11"),
        gt("10", "11"),
        MENDELIAN,
        0.0);
}

func trioBadTheta() {
    paternityIndexAt(fixture(), gt("10", "11"), gt("10", "11"), gt("10", "11"), MENDELIAN, 1.0);
}

func trioTypedSilent() {
    paternityIndexAt(
        fixture(),
        gt("10", "11"),
        forensics.Genotype{locus: "L", a: "10", b: SILENT_ALLELE},
        gt("10", "11"),
        MENDELIAN,
        0.0);
}

func testTrioErrors() {
    testing.assertThrows("maternalExclusion", "mutation");
    testing.assertThrows("trioAtMixedLoci", "mutation");
    testing.assertThrows("trioBadTheta", "mutation");
    testing.assertThrows("trioTypedSilent", "mutation");
}

func testTrioOppositeHomozygotesWithSilentAlleleByHand() {
    # mother 10/12, child shows 10, alleged father shows 11, q = .05:
    # pi10 = .95 * .1 = .095, pi11 = .95 * .2 = .19.
    # Child 10/silent, father 11/silent the only route: X = 2 pi11 q * (1/2 * 1/2).
    # Y = P(father shows 11) * T_M(10) (pi10 + q) = pi11 (pi11 + 2q) * .5 (pi10 + q).
    # PI = q / ((pi11 + 2q)(pi10 + q)) = .05 / (.29 * .145)
    def m as forensics.Genotype init gt("10", "12");
    def c as forensics.Genotype init gt("10", "10");
    def f as forensics.Genotype init gt("11", "11");
    testing.assertEqual(paternityIndexAt(fixture(), $m, $c, $f, MENDELIAN, 0.0), 0.0);
    def silent as ParentageModel init withSilentFrequency(MENDELIAN, "L", 0.05);
    near(paternityIndexAt(fixture(), $m, $c, $f, $silent, 0.0), 0.05 / (0.29 * 0.145));
}

func testTrioSilentAlleleWithThetaByHand() {
    # As above at theta t = .1, D(x | S) = (m t + (1-t) pi_x) / (1 + (n-1) t).
    # Mother 10/12; father 11/11 or 11/s. Pair weights o(Gm) o(Gf) Ptheta(10,12,f1,f2)
    # share F = 2 * pi10 * (1-t) pi12 * (1-t) pi11 / (1+t), then
    #   11/11: a = (t + (1-t) pi11) / (1+2t)     11/s: b = 2 (1-t) q / (1+2t)
    # X: only father 11/s and child 10/s, T_M(10) T_F(s) = 1/4, so X = F b / 4.
    # Y: child 10/10 or 10/s from mother 10/12 and a random man drawn given
    # the four parental alleles S: .5 D(10|S) + .5 D(s|S), with
    #   D(10|S) = (t + (1-t) pi10) / (1+3t),
    #   D(s|S) = (1-t) q / (1+3t) for S = 10,12,11,11, (t + (1-t) q) / (1+3t) for 10,12,11,s.
    # F cancels.
    def t as float init 0.1;
    def q as float init 0.05;
    def p10 as float init 0.095;
    def p11 as float init 0.19;
    def a as float init ($t + (1.0 - $t) * $p11) / (1.0 + 2.0 * $t);
    def b as float init 2.0 * (1.0 - $t) * $q / (1.0 + 2.0 * $t);
    def d10 as float init ($t + (1.0 - $t) * $p10) / (1.0 + 3.0 * $t);
    def x as float init $b * 0.25;
    def y as float init $a * 0.5 * ($d10 + (1.0 - $t) * $q / (1.0 + 3.0 * $t)) +
        $b * 0.5 * ($d10 + ($t + (1.0 - $t) * $q) / (1.0 + 3.0 * $t));
    def silent as ParentageModel init withSilentFrequency(MENDELIAN, "L", $q);
    near(
        paternityIndexAt(fixture(), gt("10", "12"), gt("10", "10"), gt("11", "11"), $silent, $t),
        $x / $y);
}

func testCombinedTrioIsTheProductOverLoci() {
    def mother as forensics.Profile init forensics.profile(
        "m",
        [gt("10", "11"), forensics.genotype("M", "10", "12")]);
    def child as forensics.Profile init forensics.profile(
        "c",
        [gt("11", "13"), forensics.genotype("M", "10", "10")]);
    def father as forensics.Profile init forensics.profile(
        "f",
        [gt("12", "12"), forensics.genotype("M", "11", "11")]);
    def model as ParentageModel init withSilentFrequency(
        parentageModel(withLocusRate(equalModel(0.003), "M", 0.0), NO_MUTATION),
        "M",
        0.05);
    # L: .5 * .001 / (.5 * .4) as in the paternal Equal case; M (no mutation
    # there): the silent-allele case, .05 / (.29 * .145)
    near(
        paternityIndex($mother, $child, $father, twoLocus(), $model, 0.0),
        0.0025 * 0.05 / (0.29 * 0.145));
}

func trioNoSharedLocus() {
    paternityIndex(
        forensics.profile("m", [gt("10", "11")]),
        forensics.profile("c", [forensics.genotype("M", "10", "11")]),
        forensics.profile("f", [gt("10", "11")]),
        twoLocus(),
        MENDELIAN,
        0.0);
}

func trioUncoveredLocus() {
    def g as forensics.Genotype init forensics.genotype("M", "10", "11");
    paternityIndex(
        forensics.profile("m", [$g]),
        forensics.profile("c", [$g]),
        forensics.profile("f", [$g]),
        fixture(),
        MENDELIAN,
        0.0);
}

func trioProfilesBadTheta() {
    def p as forensics.Profile init forensics.profile("p", [gt("10", "11")]);
    paternityIndex($p, $p, $p, fixture(), MENDELIAN, -0.1);
}

func testCombinedTrioErrors() {
    testing.assertThrows("trioNoSharedLocus", "mutation");
    testing.assertThrows("trioUncoveredLocus", "mutation");
    testing.assertThrows("trioProfilesBadTheta", "mutation");
}

# --- duo --------------------------------------------------------------------

func testMendelianDuoIsExactlyKinship() {
    def pairs as list of list of forensics.Genotype init [
        [gt("10", "11"), gt("11", "12")],
        [gt("10", "10"), gt("10", "11")],
        [gt("12", "12"), gt("12", "12")],
        [gt("11", "13"), gt("12", "12")]
    ];
    for (def p in $pairs) {
        for (def theta in [0.0, 0.01, 0.1]) {
            testing.assertEqual(
                duoPaternityIndexAt(fixture(), $p[0], $p[1], MENDELIAN, $theta),
                kinship.duoPaternityIndexAt(fixture(), $p[0], $p[1], kinship.NO_MUTATION, $theta));
        }
    }
}

func testDuoPaternalMutationByHand() {
    # child 11/13, alleged father 12/12, Equal .003 (M[12][x] = .001):
    # X = T_F(11) p13 + T_F(13) p11 = .001 * .4 + .001 * .2
    # Y = 2 p11 p13 = 2 * .2 * .4. The maternal model plays no part.
    def model as ParentageModel init parentageModel(equalModel(0.003), equalModel(0.5));
    near(
        duoPaternityIndexAt(fixture(), gt("11", "13"), gt("12", "12"), $model, 0.0),
        (0.001 * 0.4 + 0.001 * 0.2) / (2.0 * 0.2 * 0.4));
}

func testDuoOppositeHomozygotesWithSilentAlleleByHand() {
    # child shows 10, alleged father shows 11, q = .05, pi10 = .095, pi11 = .19.
    # Paternity: father 11/s (2 pi11 q) passes s (1/2), mother passes 10 (pi10):
    #   pi10 pi11 q. Unrelated: pi10 (pi10 + 2q) * pi11 (pi11 + 2q).
    # PI = q / ((pi10 + 2q)(pi11 + 2q)) = .05 / (.195 * .29)
    def c as forensics.Genotype init gt("10", "10");
    def f as forensics.Genotype init gt("11", "11");
    testing.assertEqual(duoPaternityIndexAt(fixture(), $c, $f, MENDELIAN, 0.0), 0.0);
    def silent as ParentageModel init withSilentFrequency(MENDELIAN, "L", 0.05);
    near(duoPaternityIndexAt(fixture(), $c, $f, $silent, 0.0), 0.05 / (0.195 * 0.29));
    # the pairwise parent/child LR with the same silent allele is the same
    # quantity by a different code path
    near(
        kinshipLRAt(fixture(), $c, $f, kinship.PARENT_CHILD, kinship.UNRELATED, 0.05, 0.0),
        0.05 / (0.195 * 0.29));
}

func testDuoSilentAlleleWithThetaByHand() {
    # child shows 10 (10/10 or 10/s), father shows 11 (11/11 or 11/s), t = .1.
    # w(11/11) = pi11 (t + (1-t) pi11), w(11/s) = 2 pi11 (1-t) q.
    # X: only father 11/s, child 10/s: T_F(s) D(10 | 11,s) = .5 (1-t) pi10 / (1+t).
    # Y for father genotype f: D(10|f) [D(10|f,10) + 2 D(s|f,10)], with
    #   D(10|f) = (1-t) pi10 / (1+t), D(10|f,10) = (t + (1-t) pi10) / (1+2t),
    #   D(s|f,10) = (1-t) q / (1+2t) for f = 11/11, (t + (1-t) q) / (1+2t) for 11/s.
    def t as float init 0.1;
    def q as float init 0.05;
    def p10 as float init 0.095;
    def p11 as float init 0.19;
    def w1 as float init $p11 * ($t + (1.0 - $t) * $p11);
    def w2 as float init 2.0 * $p11 * (1.0 - $t) * $q;
    def d10 as float init (1.0 - $t) * $p10 / (1.0 + $t);
    def again as float init ($t + (1.0 - $t) * $p10) / (1.0 + 2.0 * $t);
    def x as float init $w2 * 0.5 * $d10;
    def y as float init $w1 * $d10 * ($again + 2.0 * (1.0 - $t) * $q / (1.0 + 2.0 * $t)) +
        $w2 * $d10 * ($again + 2.0 * ($t + (1.0 - $t) * $q) / (1.0 + 2.0 * $t));
    def silent as ParentageModel init withSilentFrequency(MENDELIAN, "L", $q);
    near(duoPaternityIndexAt(fixture(), gt("10", "10"), gt("11", "11"), $silent, $t), $x / $y);
    near(
        kinshipLRAt(
            fixture(),
            gt("10", "10"),
            gt("11", "11"),
            kinship.PARENT_CHILD,
            kinship.UNRELATED,
            $q,
            $t),
        $x / $y);
}

func testCombinedDuoIsTheProductOverLoci() {
    def child as forensics.Profile init forensics.profile(
        "c",
        [gt("11", "13"), forensics.genotype("M", "10", "10")]);
    def father as forensics.Profile init forensics.profile(
        "f",
        [gt("12", "12"), forensics.genotype("M", "11", "11")]);
    def model as ParentageModel init withSilentFrequency(
        parentageModel(withLocusRate(equalModel(0.003), "M", 0.0), NO_MUTATION),
        "M",
        0.05);
    near(
        duoPaternityIndex($child, $father, twoLocus(), $model, 0.0),
        (0.001 * 0.4 + 0.001 * 0.2) / (2.0 * 0.2 * 0.4) * 0.05 / (0.195 * 0.29));
}

func duoMixedLoci() {
    duoPaternityIndexAt(
        fixture(),
        gt("10", "11"),
        forensics.genotype("M", "10", "11"),
        MENDELIAN,
        0.0);
}

func duoBadTheta() {
    duoPaternityIndexAt(fixture(), gt("10", "11"), gt("10", "11"), MENDELIAN, 1.0);
}

func duoNoSharedLocus() {
    duoPaternityIndex(
        forensics.profile("c", [gt("10", "11")]),
        forensics.profile("f", [forensics.genotype("M", "10", "11")]),
        twoLocus(),
        MENDELIAN,
        0.0);
}

func duoProfilesBadTheta() {
    def p as forensics.Profile init forensics.profile("p", [gt("10", "11")]);
    duoPaternityIndex($p, $p, fixture(), MENDELIAN, 2.0);
}

func testDuoErrors() {
    testing.assertThrows("duoMixedLoci", "mutation");
    testing.assertThrows("duoBadTheta", "mutation");
    testing.assertThrows("duoNoSharedLocus", "mutation");
    testing.assertThrows("duoProfilesBadTheta", "mutation");
}

# --- pairwise kinship -------------------------------------------------------

func testKinshipWithoutSilentIsExactlyKinship() {
    def a as forensics.Profile init forensics.profile(
        "a",
        [gt("10", "11"), forensics.genotype("M", "12", "12")]);
    def b as forensics.Profile init forensics.profile(
        "b",
        [gt("10", "12"), forensics.genotype("M", "12", "13")]);
    for (def rel in [kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, kinship.PARENT_CHILD]) {
        for (def theta in [0.0, 0.03]) {
            testing.assertEqual(
                kinshipLR($a, $b, $rel, kinship.UNRELATED, twoLocus(), {}, $theta),
                kinship.kinshipLRAgainstTheta($a, $b, $rel, kinship.UNRELATED, twoLocus(), $theta));
            testing.assertEqual(
                kinshipLRAt(
                    fixture(),
                    gt("10", "10"),
                    gt("10", "11"),
                    $rel,
                    kinship.UNRELATED,
                    0.0,
                    $theta),
                kinship.kinshipLRAtTheta(
                    fixture(),
                    gt("10", "10"),
                    gt("10", "11"),
                    $rel,
                    kinship.UNRELATED,
                    $theta));
        }
    }
}

func testGeneralJointAgreesWithKinshipWithoutSilent() {
    # The silent-allele path itself, at q = 0, against kinship.j's joint form.
    def la as LocusAlleles init allelesL(0.0);
    def idx as map of string to int init indexOfAllele($la.alleles);
    for (def theta in [0.0, 0.05]) {
        for (def pair in [[gt("10", "10"), gt("10", "11")], [gt("11", "12"), gt("11", "12")]]) {
            near(
                phenotypeJointTheta(
                    $la,
                    $idx,
                    $pair[0],
                    $pair[1],
                    kinship.FULL_SIBLINGS,
                    0.0,
                    $theta),
                kinship.jointGenotypeProbabilityTheta(
                    fixture(),
                    $pair[0],
                    $pair[1],
                    kinship.FULL_SIBLINGS,
                    $theta));
        }
    }
}

func testPhenotypeJointSumsToTheMarginal() {
    # A shows 10. Summed over every phenotype B can show, P(A, B) is P(A shows
    # 10) = pi10^2 + 2 pi10 q = .095^2 + 2 * .095 * .05 = .018525 at theta 0;
    # at theta .1, pi10 (t + (1-t) pi10) + 2 pi10 (1-t) q
    # = .095 * .1855 + 2 * .095 * .9 * .05 = .0261725.
    def la as LocusAlleles init allelesL(0.05);
    def idx as map of string to int init indexOfAllele($la.alleles);
    def labels as list of string init ["10", "11", "12", "13"];
    def want as list of float init [0.018525, 0.0261725];
    def thetas as list of float init [0.0, 0.1];
    for (def k in 0..2) {
        for (def rel in [kinship.FULL_SIBLINGS, kinship.PARENT_CHILD, kinship.UNRELATED]) {
            def total as float init 0.0;
            for (def i in 0..4) {
                for (def j in $i..4) {
                    $total = $total +
                        phenotypeJointTheta(
                        $la,
                        $idx,
                        gt("10", "10"),
                        gt($labels[$i], $labels[$j]),
                        $rel,
                        0.05,
                        $thetas[$k]);
                }
            }
            # B = silent/silent shows no peak, a phenotype no typing reports;
            # its share comes straight from trueJoint
            def ss as forensics.Genotype init forensics.Genotype{
                locus: "L",
                a: SILENT_ALLELE,
                b: SILENT_ALLELE
            };
            for (def ta in expand(gt("10", "10"), 0.05)) {
                $total = $total + trueJoint($la, $idx, $ta, $ss, $rel, $thetas[$k]);
            }
            near($total, $want[$k]);
        }
    }
}

func testSilentKinshipIsContinuousAtZero() {
    def lr0 as float init kinshipLRAt(
        fixture(),
        gt("10", "10"),
        gt("10", "11"),
        kinship.FULL_SIBLINGS,
        kinship.UNRELATED,
        0.0,
        0.0);
    def lr as float init kinshipLRAt(
        fixture(),
        gt("10", "10"),
        gt("10", "11"),
        kinship.FULL_SIBLINGS,
        kinship.UNRELATED,
        0.000000001,
        0.0);
    testing.assertTrue(math.abs($lr - $lr0) < 0.000001 * $lr0);
}

func kinshipZeroAlternative() {
    kinshipLRAt(
        fixture(),
        gt("10", "10"),
        gt("11", "11"),
        kinship.UNRELATED,
        kinship.PARENT_CHILD,
        0.0,
        0.0);
}

func kinshipZeroAlternativeWithSilent() {
    kinshipLRAt(
        fixture(),
        gt("10", "11"),
        gt("12", "13"),
        kinship.UNRELATED,
        kinship.IDENTICAL,
        0.05,
        0.0);
}

func kinshipMixedLoci() {
    kinshipLRAt(
        fixture(),
        gt("10", "11"),
        forensics.genotype("M", "10", "11"),
        kinship.FULL_SIBLINGS,
        kinship.UNRELATED,
        0.05,
        0.0);
}

func kinshipBadSilent() {
    kinshipLRAt(
        fixture(),
        gt("10", "11"),
        gt("10", "11"),
        kinship.FULL_SIBLINGS,
        kinship.UNRELATED,
        -0.1,
        0.0);
}

func kinshipNoSharedLocus() {
    kinshipLR(
        forensics.profile("a", [gt("10", "11")]),
        forensics.profile("b", [forensics.genotype("M", "10", "11")]),
        kinship.FULL_SIBLINGS,
        kinship.UNRELATED,
        twoLocus(),
        {},
        0.0);
}

func kinshipProfilesBadTheta() {
    def p as forensics.Profile init forensics.profile("p", [gt("10", "11")]);
    kinshipLR($p, $p, kinship.FULL_SIBLINGS, kinship.UNRELATED, fixture(), {}, 1.0);
}

func testKinshipErrors() {
    # at q = 0 the call goes to kinship.j and keeps its error kind
    testing.assertThrows("kinshipZeroAlternative", "kinship");
    testing.assertThrows("kinshipZeroAlternativeWithSilent", "mutation");
    testing.assertThrows("kinshipMixedLoci", "mutation");
    testing.assertThrows("kinshipBadSilent", "mutation");
    testing.assertThrows("kinshipNoSharedLocus", "mutation");
    testing.assertThrows("kinshipProfilesBadTheta", "mutation");
}

# --- published counts -------------------------------------------------------

func testAabbCountsAreInternallyConsistent() {
    # Every row is printed with its parts and its totals: mutations
    # maternal + paternal + either = total, meioses maternal + paternal = total.
    def counts as list of MutationCount init aabb2003Counts();
    testing.assertEqual(len($counts), 18);
    for (def c in $counts) {
        testing.assertEqual(
            $c.maternalMutations + $c.paternalMutations + $c.unassignedMutations,
            $c.totalMutations);
        testing.assertEqual($c.maternalMeioses + $c.paternalMeioses, $c.totalMeioses);
    }
}

func testAabbTotalRatesRoundToThePublishedColumn() {
    # STRBase's "Mutation Rate" column in percent, two decimals, page order
    def published as list of int init [
        16,
        28,
        1,
        1,
        17,
        12,
        11,
        10,
        14,
        14,
        11,
        22,
        19,
        14,
        16,
        12,
        11,
        64
    ];
    def rates as map of string to float init aabb2003Rates("total");
    def counts as list of MutationCount init aabb2003Counts();
    for (def i in 0..len($counts)) {
        testing.assertEqual(math.round($rates[$counts[$i].locus] * 10000.0), $published[$i]);
    }
    near($rates["CSF1PO"], 1487.0 / 947425.0);
    near(aabb2003Rates("paternal")["FGA"], 2210.0 / 692776.0);
    near(aabb2003Rates("maternal")["SE33"], 0.0);
    near(aabb2003Rates("maternal")["D21S11"], 464.0 / 435388.0);
}

func badRateKind() {
    aabb2003Rates("female");
}

func testAabbRatesRejectAnUnknownKind() {
    testing.assertThrows("badRateKind", "mutation");
}

func testPublishedRatesReachBothModelKinds() {
    def rates as map of string to float init aabb2003Rates("total");
    def spec as MutationSpec init withLocusRates(stepwiseModel(0.002, 0.1), $rates);
    testing.assertEqual(len($spec.perLocus), 18);
    near($spec.perLocus["TH01"], 100.0 / 779554.0);
    def km as kinship.MutationModel init kinshipModelWithRates(kinship.DEFAULT_MUTATION, $rates);
    near(kinship.mutationProbability($km, "TH01", "8", "8"), 1.0 - 100.0 / 779554.0);
    # a locus outside the table keeps the blanket rate
    near(kinship.mutationProbability($km, "D1S1656", "8", "8"), 0.998);
}

func badKinshipRate() {
    kinshipModelWithRates(kinship.DEFAULT_MUTATION, {"TH01": 2.0});
}

func testKinshipModelWithRatesValidates() {
    testing.assertThrows("badKinshipRate", "kinship");
}

# --- brute force ------------------------------------------------------------

# shows says whether ordered true alleles (x, y) display as typed genotype g:
# a heterozygote exactly, a single peak a as a/a or a/silent.
func shows(x as string, y as string, g as forensics.Genotype) {
    if ($g.a != $g.b) {
        return ($x == $g.a and $y == $g.b) or ($x == $g.b and $y == $g.a);
    }
    if ($x == $g.a) {
        return $y == $g.a or $y == SILENT_ALLELE;
    }
    return $y == $g.a and $x == SILENT_ALLELE;
}

# bruteTrio enumerates every ordered allele pair of mother, alleged father and
# child as Familias does (founders drawn from the frequencies, the child's
# paternal allele through the paternal matrix from the alleged father, or
# straight from the population for the random man), keeping the combinations
# each typed genotype can show. It shares only the matrices with the code.
func bruteTrio(
    la as LocusAlleles,
    fm as MutationMatrix,
    mm as MutationMatrix,
    gs as list of forensics.Genotype) {
    def n as int init len($la.alleles);
    def a as list of string init $la.alleles;
    def pi as list of float init $la.frequencies;
    def num as float init 0.0;
    def den as float init 0.0;
    for (def m1 in 0..$n) {
        for (def m2 in 0..$n) {
            if (shows($a[$m1], $a[$m2], $gs[0])) {
                for (def f1 in 0..$n) {
                    for (def f2 in 0..$n) {
                        if (shows($a[$f1], $a[$f2], $gs[2])) {
                            def w as float init $pi[$m1] * $pi[$m2] * $pi[$f1] * $pi[$f2];
                            for (def cp in 0..$n) {
                                for (def cm in 0..$n) {
                                    if (shows($a[$cp], $a[$cm], $gs[1])) {
                                        def tm as float init 0.5 * $mm.rows[$m1][$cm] +
                                            0.5 * $mm.rows[$m2][$cm];
                                        def tf as float init 0.5 * $fm.rows[$f1][$cp] +
                                            0.5 * $fm.rows[$f2][$cp];
                                        $num = $num + $w * $tm * $tf;
                                        $den = $den + $w * $tm * $pi[$cp];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return $num / $den;
}

func testTrioMatchesBruteForceEnumeration() {
    def db as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("b"),
        "L",
        10000,
        {"9": 0.15, "10": 0.1, "10.2": 0.05, "11": 0.2, "12": 0.3, "13": 0.2});
    def model as ParentageModel init withSilentFrequency(
        parentageModel(
            stationaryModel(extendedStepwiseModel(0.004, 0.4, 0.001)),
            proportionalModel(0.002)),
        "L",
        0.03);
    def la as LocusAlleles init locusAlleles($db, "L", [], 0.03);
    def fm as MutationMatrix init mutationMatrix($model.paternal, $la);
    def mm as MutationMatrix init mutationMatrix($model.maternal, $la);
    def trios as list of list of forensics.Genotype init [
        [gt("10", "12"), gt("10", "10"), gt("13", "13")],
        [gt("9", "11"), gt("11", "13"), gt("12", "12")],
        [gt("10.2", "13"), gt("13", "13"), gt("11", "12")],
        [gt("12", "12"), gt("11", "11"), gt("13", "13")]
    ];
    for (def t in $trios) {
        near(paternityIndexAt($db, $t[0], $t[1], $t[2], $model, 0.0), bruteTrio($la, $fm, $mm, $t));
    }
}
