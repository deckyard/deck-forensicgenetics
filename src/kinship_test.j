# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for kinship.j. Run with `jennifer test kinship_test.j`.
#
# The reference table is deliberately tiny - one locus "L" with alleles
# 10/11/12/13 at frequencies .1/.2/.3/.4 and N = 10000, so the 5/2N floor
# (0.00025) never bites and every expectation below is exactly hand-derivable
# from the textbook formulae.
use testing;
use math;

def const TOL as float init 0.000000001;

# near asserts two floats agree to within a relative tolerance, so a likelihood
# ratio in the hundreds is judged as strictly as a probability near zero.
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

# twoLocus repeats locus L as "M" so a product across loci can be checked.
func twoLocus() {
    def out as forensics.FrequencyDb init fixture();
    return forensics.withLocus($out, "M", 10000, {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

func gt(a as string, b as string) {
    return forensics.genotype("L", $a, $b);
}

# --- relationships ----------------------------------------------------------

func testStandardRelationshipCoefficients() {
    near(PARENT_CHILD.k1, 1.0);
    near(FULL_SIBLINGS.k0, 0.25);
    near(FULL_SIBLINGS.k1, 0.5);
    near(FULL_SIBLINGS.k2, 0.25);
    near(FIRST_COUSINS.k0, 0.75);
    near(UNRELATED.k0, 1.0);
    near(IDENTICAL.k2, 1.0);
}

func testHalfSibGrandparentAvuncularShareTheirCoefficients() {
    # Autosomal STRs cannot separate these three - the LR is identical, and
    # that is a property of the markers, not a shortcut in this code.
    near(HALF_SIBLINGS.k0, GRANDPARENT.k0);
    near(HALF_SIBLINGS.k1, AVUNCULAR.k1);
    def a as forensics.Profile init forensics.profile("A", [gt("10", "11")]);
    def b as forensics.Profile init forensics.profile("B", [gt("10", "12")]);
    near(kinshipLR($a, $b, HALF_SIBLINGS, fixture()), kinshipLR($a, $b, GRANDPARENT, fixture()));
}

func testCustomRelationship() {
    def r as Relationship init relationship("double first cousins", 0.5625, 0.375, 0.0625);
    near($r.k1, 0.375);
}

func coefficientsDoNotSumToOne() {
    relationship("bad", 0.5, 0.4, 0.2);
}

func negativeCoefficient() {
    relationship("bad", -0.1, 0.6, 0.5);
}

func testRelationshipCoefficientsAreValidated() {
    testing.assertThrows("coefficientsDoNotSumToOne", "kinship");
    testing.assertThrows("negativeCoefficient", "kinship");
}

# --- transmission -----------------------------------------------------------

func testTransmissionProbabilityIsMendelian() {
    near(transmissionProbability(gt("10", "10"), "10"), 1.0);
    near(transmissionProbability(gt("10", "11"), "10"), 0.5);
    near(transmissionProbability(gt("11", "12"), "10"), 0.0);
}

# --- joint genotype probabilities -------------------------------------------

func testJointProbabilityUnderParentChild() {
    # P(ga) * p11 = .01 * .2
    near(jointGenotypeProbability(fixture(), gt("10", "10"), gt("10", "11"), PARENT_CHILD), 0.002);
}

func testJointProbabilityUnderUnrelatedIsTheProduct() {
    near(
        jointGenotypeProbability(fixture(), gt("10", "10"), gt("10", "11"), UNRELATED),
        0.01 * 0.04);
}

func testJointProbabilityUnderIdentical() {
    near(jointGenotypeProbability(fixture(), gt("10", "11"), gt("10", "11"), IDENTICAL), 0.04);
    # Two different genotypes cannot both belong to the same person.
    near(jointGenotypeProbability(fixture(), gt("10", "11"), gt("10", "12"), IDENTICAL), 0.0);
}

func testJointProbabilitiesSumToTheMarginal() {
    # Summing over every genotype the second person could have must give back
    # P(first person's genotype) - the check that the IBD-1 term is a proper
    # conditional distribution and not merely plausible.
    def fdb as forensics.FrequencyDb init fixture();
    def labels as list of string init ["10", "11", "12", "13"];
    def total as float init 0.0;
    for (def i in 0..len($labels)) {
        for (def j in $i..len($labels)) {
            $total = $total +
                jointGenotypeProbability(
                $fdb,
                gt("10", "11"),
                gt($labels[$i], $labels[$j]),
                FULL_SIBLINGS);
        }
    }
    near($total, 0.04);
}

func mismatchedLoci() {
    jointGenotypeProbability(
        fixture(),
        gt("10", "11"),
        forensics.genotype("OTHER", "10", "11"),
        PARENT_CHILD);
}

func testJointProbabilityRejectsMismatchedLoci() {
    testing.assertThrows("mismatchedLoci", "kinship");
}

# --- kinship likelihood ratios ----------------------------------------------

func testParentChildLikelihoodRatio() {
    # Father 10/10, child 10/11: LR = 1 / (2 p10) = 5.
    near(kinshipLRAt(fixture(), gt("10", "10"), gt("10", "11"), PARENT_CHILD, UNRELATED), 5.0);
}

func testKinshipLrIsSymmetric() {
    def fdb as forensics.FrequencyDb init fixture();
    near(
        kinshipLRAt($fdb, gt("10", "10"), gt("10", "11"), PARENT_CHILD, UNRELATED),
        kinshipLRAt($fdb, gt("10", "11"), gt("10", "10"), PARENT_CHILD, UNRELATED));
    near(
        kinshipLRAt($fdb, gt("10", "11"), gt("12", "13"), FULL_SIBLINGS, UNRELATED),
        kinshipLRAt($fdb, gt("12", "13"), gt("10", "11"), FULL_SIBLINGS, UNRELATED));
}

func testFullSiblingLikelihoodRatio() {
    # Both 10/10: (1 + 1/p)^2 / 4 = (1 + 10)^2 / 4 = 30.25.
    near(kinshipLRAt(fixture(), gt("10", "10"), gt("10", "10"), FULL_SIBLINGS, UNRELATED), 30.25);
    # Both 10/11: 1/4 + (1 + p10 + p11) / (8 p10 p11) = .25 + 1.3/.16 = 8.375.
    near(kinshipLRAt(fixture(), gt("10", "11"), gt("10", "11"), FULL_SIBLINGS, UNRELATED), 8.375);
}

func testIdenticalLikelihoodRatioIsTheReciprocalGenotypeFrequency() {
    near(kinshipLRAt(fixture(), gt("10", "11"), gt("10", "11"), IDENTICAL, UNRELATED), 25.0);
}

func testSharingNoAlleleFavoursUnrelatedness() {
    # Two people with no allele in common cannot be parent and child at all.
    near(kinshipLRAt(fixture(), gt("10", "11"), gt("12", "13"), PARENT_CHILD, UNRELATED), 0.0);
    # Full siblings can still share nothing, but it is evidence against.
    testing.assertTrue(kinshipLRAt(
        fixture(),
        gt("10", "11"),
        gt("12", "13"),
        FULL_SIBLINGS,
        UNRELATED) < 1.0);
}

func testUnrelatedAgainstUnrelatedIsOne() {
    near(kinshipLRAt(fixture(), gt("10", "11"), gt("10", "12"), UNRELATED, UNRELATED), 1.0);
}

func testOneRelationshipAgainstAnother() {
    def a as forensics.Profile init forensics.profile("A", [gt("10", "10")]);
    def b as forensics.Profile init forensics.profile("B", [gt("10", "10")]);
    near(kinshipLRAgainst($a, $b, FULL_SIBLINGS, HALF_SIBLINGS, fixture()), 5.5);
}

func testKinshipLrMultipliesAcrossLoci() {
    def a as forensics.Profile init forensics.profile(
        "A",
        [gt("10", "10"), forensics.genotype("M", "10", "11")]);
    def b as forensics.Profile init forensics.profile(
        "B",
        [gt("10", "10"), forensics.genotype("M", "10", "11")]);
    # 30.25 * 8.375
    near(kinshipLRAgainst($a, $b, FULL_SIBLINGS, UNRELATED, twoLocus()), 253.34375);
}

func testPairLociIsTheIntersection() {
    def a as forensics.Profile init forensics.profile(
        "A",
        [gt("10", "10"), forensics.genotype("M", "10", "11")]);
    def b as forensics.Profile init forensics.profile("B", [gt("10", "10")]);
    testing.assertEqual(pairLoci(twoLocus(), $a, $b), ["L"]);
    # "M" is dropped because only A carries it; "L" alone drives the LR.
    near(kinshipLR($a, $b, FULL_SIBLINGS, twoLocus()), 30.25);
}

func noSharedLocus() {
    def a as forensics.Profile init forensics.profile("A", [gt("10", "10")]);
    def b as forensics.Profile init forensics.profile("B", [forensics.genotype("M", "10", "11")]);
    kinshipLR($a, $b, FULL_SIBLINGS, fixture());
}

func testNoSharedLocusThrows() {
    testing.assertThrows("noSharedLocus", "kinship");
}

# --- mutation model ---------------------------------------------------------

func testMutationProbabilityOfNoChange() {
    near(mutationProbability(DEFAULT_MUTATION, "L", "12", "12"), 0.998);
}

func testMutationProbabilityDecaysWithStepCount() {
    # rate * (1-decay) * decay^(k-1) / 2
    near(mutationProbability(DEFAULT_MUTATION, "L", "13", "12"), 0.0009);
    near(mutationProbability(DEFAULT_MUTATION, "L", "12", "10"), 0.00009);
    near(mutationProbability(DEFAULT_MUTATION, "L", "13", "10"), 0.000009);
}

func testMutationStepsSumToTheRate() {
    # Every non-zero change, both directions, must add up to exactly `rate`.
    def total as float init 0.0;
    for (def k in 1..40) {
        def target as string init convert.toString(12 + $k);
        $total = $total + 2.0 * mutationProbability(DEFAULT_MUTATION, "L", "12", $target);
    }
    near($total, DEFAULT_MUTATION.rate);
}

func testMicrovariantCountsAsOneStep() {
    near(mutationProbability(DEFAULT_MUTATION, "L", "15", "15.2"), 0.0009);
}

func testNoMutationModelIsStrictlyMendelian() {
    near(mutationProbability(NO_MUTATION, "L", "12", "12"), 1.0);
    near(mutationProbability(NO_MUTATION, "L", "13", "12"), 0.0);
}

func testPerLocusRateOverridesTheBlanketRate() {
    def m as MutationModel init withLocusRate(DEFAULT_MUTATION, "L", 0.01);
    near(mutationProbability($m, "L", "12", "12"), 0.99);
    # Another locus keeps the blanket rate.
    near(mutationProbability($m, "M", "12", "12"), 0.998);
}

func mutationBetweenNonNumericAlleles() {
    mutationProbability(DEFAULT_MUTATION, "AMEL", "X", "Y");
}

func badMutationRate() {
    mutationModel(1.5, 0.1);
}

func badMutationDecay() {
    mutationModel(0.002, 1.0);
}

func testMutationModelValidatesItsInput() {
    testing.assertThrows("mutationBetweenNonNumericAlleles", "kinship");
    testing.assertThrows("badMutationRate", "kinship");
    testing.assertThrows("badMutationDecay", "kinship");
}

# --- paternity --------------------------------------------------------------

func mother() {
    return forensics.profile("M", [gt("10", "11")]);
}

func child() {
    return forensics.profile("C", [gt("10", "12")]);
}

func father() {
    return forensics.profile("F", [gt("12", "13")]);
}

func testPaternityIndexAtOneLocus() {
    # X = T(M,10) * T(F,12) = .5 * .5 = .25
    # Y = T(M,10) * p12     = .5 * .3 = .15
    near(
        paternityIndexAt(fixture(), gt("10", "11"), gt("10", "12"), gt("12", "13"), NO_MUTATION),
        1.6666666666666667);
}

func testPaternityIndexUnderTheDefaultMutationModel() {
    # The alleged father's transmission is discounted by (1 - rate), so the
    # index drops very slightly against the strict Mendelian value.
    near(
        paternityIndexWith(mother(), child(), father(), fixture(), DEFAULT_MUTATION),
        1.6648333333333334);
    near(paternityIndex(mother(), child(), father(), fixture()), 1.6648333333333334);
}

func testStrictMendelianPaternityIndex() {
    near(
        paternityIndexWith(mother(), child(), father(), fixture(), NO_MUTATION),
        1.6666666666666667);
}

func excludedFather() {
    return forensics.profile("F2", [gt("13", "13")]);
}

func testASingleInconsistencyExcludesWithoutAMutationModel() {
    # The child needs a 12 from the father; a 13/13 man cannot transmit one.
    near(paternityIndexWith(mother(), child(), excludedFather(), fixture(), NO_MUTATION), 0.0);
    testing.assertEqual(inconsistentLoci(mother(), child(), excludedFather(), fixture()), ["L"]);
}

func testTheMutationModelSurvivesALoneInconsistency() {
    # A one-step mutation explains it, at the cost of a very small index.
    near(
        paternityIndexWith(mother(), child(), excludedFather(), fixture(), DEFAULT_MUTATION),
        0.003);
    testing.assertTrue(paternityIndexWith(
        mother(),
        child(),
        excludedFather(),
        fixture(),
        DEFAULT_MUTATION) > 0.0);
}

func testConsistentTrioHasNoInconsistentLoci() {
    testing.assertEqual(inconsistentLoci(mother(), child(), father(), fixture()), []);
}

func testPaternityIndexMultipliesAcrossLoci() {
    def m as forensics.Profile init forensics.profile(
        "M",
        [gt("10", "11"), forensics.genotype("M", "10", "11")]);
    def c as forensics.Profile init forensics.profile(
        "C",
        [gt("10", "12"), forensics.genotype("M", "10", "12")]);
    def f as forensics.Profile init forensics.profile(
        "F",
        [gt("12", "13"), forensics.genotype("M", "12", "13")]);
    near(
        paternityIndexWith($m, $c, $f, twoLocus(), NO_MUTATION),
        1.6666666666666667 * 1.6666666666666667);
    testing.assertEqual(trioLoci(twoLocus(), $m, $c, $f), ["L", "M"]);
}

func maternalExclusion() {
    def m as forensics.Profile init forensics.profile("M", [gt("13", "13")]);
    paternityIndexWith($m, child(), father(), fixture(), NO_MUTATION);
}

func testMaternalExclusionIsAnError() {
    # A child carrying neither of the mother's alleles is a sample or pedigree
    # problem, not a paternity index of zero.
    testing.assertThrows("maternalExclusion", "kinship");
}

func trioWithNoSharedLocus() {
    def f as forensics.Profile init forensics.profile("F", [forensics.genotype("M", "12", "13")]);
    paternityIndex(mother(), child(), $f, twoLocus());
}

func testTrioWithNoSharedLocusThrows() {
    testing.assertThrows("trioWithNoSharedLocus", "kinship");
}

# --- probability of paternity -----------------------------------------------

func testProbabilityOfPaternityUsesANeutralPrior() {
    near(probabilityOfPaternity(1000.0), 0.999000999000999);
    near(probabilityOfPaternity(1.0), 0.5);
    near(probabilityOfPaternity(0.0), 0.0);
}

func testProbabilityOfPaternityUnderAnExplicitPrior() {
    near(probabilityOfPaternityPrior(1000.0, 0.1), 0.9910802775024776);
    near(probabilityOfPaternityPrior(1000.0, 0.5), probabilityOfPaternity(1000.0));
}

func negativeCpi() {
    probabilityOfPaternity(-1.0);
}

func priorOfOne() {
    probabilityOfPaternityPrior(1000.0, 1.0);
}

func testProbabilityOfPaternityValidatesItsInput() {
    testing.assertThrows("negativeCpi", "kinship");
    testing.assertThrows("priorOfOne", "kinship");
}

# --- childProbability must be a probability, not twice one -------------------

func testChildProbabilityIsAProbability() {
    # Mother 10/11, child 10/10, and a father-side transmission of 0.5 for the
    # 10. The child's two alleles are the same allele, so there is one way to
    # assign them, not two: P = T(mother, 10) * 0.5 = 0.5 * 0.5.
    #
    # Dropping the homozygous branch and letting the two-term heterozygote form
    # handle this case returns 0.5 - twice a probability. That cancels in the
    # ratio `paternityIndexAt` takes, so no exported result moves and no other
    # test here would notice; this one would.
    near(childProbability(gt("10", "11"), gt("10", "10"), {"10": 0.5}), 0.25);
    # The heterozygous case, for contrast: both assignments are real events.
    near(childProbability(gt("10", "11"), gt("10", "12"), {"10": 0.5, "12": 0.5}), 0.25);
}

# --- documented edge cases with no assertion yet ----------------------------

func locusRateOfOne() {
    withLocusRate(DEFAULT_MUTATION, "L", 1.0);
}

func negativeLocusRate() {
    withLocusRate(DEFAULT_MUTATION, "L", -0.001);
}

func testWithLocusRateValidatesTheRate() {
    testing.assertThrows("locusRateOfOne", "kinship");
    testing.assertThrows("negativeLocusRate", "kinship");
}

func zeroDenominator() {
    # Parent and child must share an allele; these two share none.
    kinshipLRAt(fixture(), gt("10", "11"), gt("12", "13"), UNRELATED, PARENT_CHILD);
}

func testKinshipLrAtRefusesAZeroDenominator() {
    testing.assertThrows("zeroDenominator", "kinship");
}

func trioAtMixedLoci() {
    paternityIndexAt(
        fixture(),
        gt("10", "11"),
        gt("10", "12"),
        forensics.genotype("M", "12", "13"),
        NO_MUTATION);
}

func testPaternityIndexAtRejectsMixedLoci() {
    testing.assertThrows("trioAtMixedLoci", "kinship");
}

func priorOfZero() {
    probabilityOfPaternityPrior(1000.0, 0.0);
}

func testPriorMustBeStrictlyInsideTheUnitInterval() {
    testing.assertThrows("priorOfZero", "kinship");
}

# --- regressions from the validation audit ----------------------------------

# A locus the people are all typed at but the table does not cover must not be
# dropped: it may be the one locus that excludes. Here X is exactly that - a
# 13/13 man cannot give the child its 12 - and the table knows only L.
func trioAtL(ida as string, idb as string, idc as string) {
    def m as forensics.Profile init forensics.profile(
        $ida,
        [gt("10", "11"), forensics.genotype("X", "10", "11")]);
    def c as forensics.Profile init forensics.profile(
        $idb,
        [gt("10", "12"), forensics.genotype("X", "10", "12")]);
    def f as forensics.Profile init forensics.profile(
        $idc,
        [gt("12", "13"), forensics.genotype("X", "13", "13")]);
    return [$m, $c, $f];
}

func paternityOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    paternityIndexWith($t[0], $t[1], $t[2], fixture(), NO_MUTATION);
}

func inconsistenciesOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    inconsistentLoci($t[0], $t[1], $t[2], fixture());
}

func kinshipOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    kinshipLR($t[1], $t[2], PARENT_CHILD, fixture());
}

func testAnUncoveredSharedLocusThrowsInsteadOfHidingAnExclusion() {
    testing.assertThrows("paternityOverAnUncoveredLocus", "kinship");
    testing.assertThrows("inconsistenciesOverAnUncoveredLocus", "kinship");
    testing.assertThrows("kinshipOverAnUncoveredLocus", "kinship");
    # The query functions still report what the table can use, so a caller
    # can narrow the profiles deliberately and say so in the report.
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    testing.assertEqual(trioLoci(fixture(), $t[0], $t[1], $t[2]), ["L"]);
    testing.assertEqual(pairLoci(fixture(), $t[1], $t[2]), ["L"]);
    def keep as list of string init ["L"];
    near(
        paternityIndexWith(
            forensics.restrict($t[0], $keep),
            forensics.restrict($t[1], $keep),
            forensics.restrict($t[2], $keep),
            fixture(),
            NO_MUTATION),
        1.6666666666666667);
}

func testALocusOnlyOneSideCarriesIsStillSkipped() {
    # Nothing is hidden when one profile simply has no genotype there.
    def a as forensics.Profile init forensics.profile(
        "A",
        [gt("10", "10"), forensics.genotype("X", "10", "11")]);
    def b as forensics.Profile init forensics.profile("B", [gt("10", "10")]);
    near(kinshipLR($a, $b, FULL_SIBLINGS, fixture()), 30.25);
}

func testJointProbabilityIgnoresAlleleOrder() {
    # A Genotype built as a literal rather than through genotype() may carry
    # its alleles in either order; the k2 term must still recognise it.
    def raw as forensics.Genotype init forensics.Genotype{locus: "L", a: "11", b: "10"};
    near(jointGenotypeProbability(fixture(), gt("10", "11"), $raw, IDENTICAL), 0.04);
    near(
        jointGenotypeProbability(fixture(), gt("10", "11"), $raw, FULL_SIBLINGS),
        jointGenotypeProbability(fixture(), gt("10", "11"), gt("10", "11"), FULL_SIBLINGS));
}

func testZeroRateNeverCountsSteps() {
    # With no mutation a changed allele has probability 0 before any step
    # count is attempted, so non-numeric designations do not throw.
    near(mutationProbability(NO_MUTATION, "AMEL", "X", "Y"), 0.0);
    near(mutationProbability(NO_MUTATION, "AMEL", "X", "X"), 1.0);
}

# --- theta: the Balding-Nichols sampling formula ----------------------------
#
# Theta 0.1 throughout, so each sampling-formula factor is a short decimal.
# D(x | S) = (m theta + (1 - theta) p_x) / (1 + (n - 1) theta) for the next
# allele x after the n alleles S, m of which are x. With theta = 0.1:
#   (1 - theta) p10 = .09, (1 - theta) p11 = .18, (1 - theta) p12 = .27
# and the denominators for n = 2, 3, 4 already drawn are 1.1, 1.2, 1.3.

def const T as float init 0.1;

func testThetaZeroVariantsEqualTheOriginalsExactly() {
    def fdb as forensics.FrequencyDb init fixture();
    def gs as list of forensics.Genotype init [gt("10", "10"), gt("10", "11"), gt("12", "13")];
    def rels as list of Relationship init [UNRELATED, PARENT_CHILD, FULL_SIBLINGS, IDENTICAL];
    for (def i in 0..len($gs)) {
        for (def j in 0..len($gs)) {
            for (def k in 0..len($rels)) {
                testing.assertEqual(
                    jointGenotypeProbabilityTheta($fdb, $gs[$i], $gs[$j], $rels[$k], 0.0),
                    jointGenotypeProbability($fdb, $gs[$i], $gs[$j], $rels[$k]));
            }
            testing.assertEqual(
                kinshipLRAtTheta($fdb, $gs[$i], $gs[$j], FULL_SIBLINGS, UNRELATED, 0.0),
                kinshipLRAt($fdb, $gs[$i], $gs[$j], FULL_SIBLINGS, UNRELATED));
        }
    }
    def a as forensics.Profile init forensics.profile("A", [gt("10", "10")]);
    def b as forensics.Profile init forensics.profile("B", [gt("10", "11")]);
    testing.assertEqual(
        kinshipLRTheta($a, $b, PARENT_CHILD, $fdb, 0.0),
        kinshipLR($a, $b, PARENT_CHILD, $fdb));
    testing.assertEqual(
        kinshipLRAgainstTheta($a, $b, FULL_SIBLINGS, HALF_SIBLINGS, $fdb, 0.0),
        kinshipLRAgainst($a, $b, FULL_SIBLINGS, HALF_SIBLINGS, $fdb));
    testing.assertEqual(
        paternityIndexAtTheta(
            $fdb,
            gt("10", "11"),
            gt("10", "12"),
            gt("12", "13"),
            DEFAULT_MUTATION,
            0.0),
        paternityIndexAt($fdb, gt("10", "11"), gt("10", "12"), gt("12", "13"), DEFAULT_MUTATION));
    testing.assertEqual(
        paternityIndexTheta(mother(), child(), excludedFather(), $fdb, 0.0),
        paternityIndex(mother(), child(), excludedFather(), $fdb));
    testing.assertEqual(
        paternityIndexWithTheta(mother(), child(), father(), $fdb, NO_MUTATION, 0.0),
        paternityIndexWith(mother(), child(), father(), $fdb, NO_MUTATION));
}

func testTheSamplingFormulaReducesToTheProductRuleAtThetaZero() {
    # The theta-0 exports route through the original code; this checks that
    # the general form itself reduces, term by term, rather than relying on
    # that routing.
    def fdb as forensics.FrequencyDb init fixture();
    def gs as list of forensics.Genotype init [gt("10", "10"), gt("10", "11"), gt("10", "12")];
    def rels as list of Relationship init [UNRELATED, PARENT_CHILD, FULL_SIBLINGS, IDENTICAL];
    for (def i in 0..len($gs)) {
        for (def j in 0..len($gs)) {
            for (def k in 0..len($rels)) {
                near(
                    jointTheta($fdb, $gs[$i], $gs[$j], $rels[$k], 0.0),
                    jointGenotypeProbability($fdb, $gs[$i], $gs[$j], $rels[$k]));
            }
        }
    }
    near(
        trioIndexTheta($fdb, gt("10", "11"), gt("10", "12"), gt("12", "13"), DEFAULT_MUTATION, 0.0),
        paternityIndexAt($fdb, gt("10", "11"), gt("10", "12"), gt("12", "13"), DEFAULT_MUTATION));
}

func testDrawProbabilityIsTheSamplingFormula() {
    # First draw: (1 - theta) p / (1 - theta) = p.
    near(drawProbability(fixture(), "L", [], "10", T), 0.1);
    # After 10, 10: (2 * .1 + .09) / 1.1 = .29 / 1.1
    near(drawProbability(fixture(), "L", ["10", "10"], "10", T), 0.29 / 1.1);
    # After 10, 11, 10, 13: (2 * .1 + .09) / 1.3; an unseen 12: .27 / 1.3
    near(drawProbability(fixture(), "L", ["10", "11", "10", "13"], "10", T), 0.29 / 1.3);
    near(drawProbability(fixture(), "L", ["10", "11", "10", "13"], "12", T), 0.27 / 1.3);
}

func testParentChildWithThetaAtAHomozygousLocus() {
    # Both 10/10. IBD 1 draws A's 10, 10 and B's free 10; unrelated draws four.
    #   T1    = P(10,10) D(10 | 10,10)
    #   P0    = P(10,10) D(10 | 10,10) D(10 | 10,10,10)
    #   LR    = 1 / D(10 | 10,10,10) = (1 + 2t) / (3t + (1-t) p10)
    #         = 1.2 / (.3 + .09) = 1.2 / .39 = 3.0769...      (theta 0: 10)
    near(
        kinshipLRAtTheta(fixture(), gt("10", "10"), gt("10", "10"), PARENT_CHILD, UNRELATED, T),
        1.2 / 0.39);
}

func testParentChildWithThetaAtAHeterozygousLocus() {
    # Both 10/11. With u = t + (1-t) p10 = .19 and v = t + (1-t) p11 = .28:
    #   T1 = 2 p10 (1-t) p11 * [1/2 D(11 | 10,11) + 1/2 D(10 | 10,11)]
    #      = 2 p10 (1-t) p11 * (v + u) / (2 (1+t))
    #   P0 = 4 p10 (1-t) p11 * D(10 | 10,11) D(11 | 10,11,10)
    #      = 4 p10 (1-t) p11 * u v / ((1+t)(1+2t))
    #   LR = (1+2t)(u + v) / (4 u v) = 1.2 * .47 / (4 * .19 * .28)
    #      = .564 / .2128 = 2.6503...                          (theta 0: 3.75)
    near(
        kinshipLRAtTheta(fixture(), gt("10", "11"), gt("10", "11"), PARENT_CHILD, UNRELATED, T),
        0.564 / 0.2128);
    # 10/11 and 10/12 share only the 10:
    #   T1 = 2 p10 (1-t) p11 * 1/2 D(12 | 10,11) = 2 p10 (1-t) p11 * .27 / (2 * 1.1)
    #   P0 = 4 p10 (1-t) p11 * D(10 | 10,11) D(12 | 10,11,10) = ... * (.19/1.1)(.27/1.2)
    #   LR = (1+2t) / (4 u) = 1.2 / .76 = 1.5789...           (theta 0: 2.5)
    near(
        kinshipLRAtTheta(fixture(), gt("10", "11"), gt("10", "12"), PARENT_CHILD, UNRELATED, T),
        1.2 / 0.76);
}

func testFullSiblingsWithThetaAtAHomozygousLocus() {
    # LR(FS) = 1/4 + 1/2 LR(PC) + 1/4 LR(ID), each against unrelated at the same
    # theta, because the joint probability is linear in (k0, k1, k2).
    # Both 10/10: LR(PC) = 1.2 / .39 (above), and
    #   LR(ID) = P(10,10) / P0 = (1+t)(1+2t) / ((2t + .09)(3t + .09))
    #          = 1.32 / (.29 * .39) = 1.32 / .1131
    #   LR(FS) = .25 + .5 * 3.0769... + .25 * 11.6710... = 4.7062...  (theta 0: 30.25)
    near(
        kinshipLRAtTheta(fixture(), gt("10", "10"), gt("10", "10"), FULL_SIBLINGS, UNRELATED, T),
        0.25 + 0.5 * (1.2 / 0.39) + 0.25 * (1.32 / 0.1131));
}

func testFullSiblingsWithThetaAtAHeterozygousLocus() {
    # Both 10/11: LR(PC) = .564 / .2128 (above), and
    #   LR(ID) = 2 p10 (1-t) p11 / P0 = (1+t)(1+2t) / (2 u v)
    #          = 1.32 / (2 * .19 * .28) = 1.32 / .1064
    #   LR(FS) = .25 + .5 * 2.6503... + .25 * 12.4060... = 4.6766...  (theta 0: 8.375)
    near(
        kinshipLRAtTheta(fixture(), gt("10", "11"), gt("10", "11"), FULL_SIBLINGS, UNRELATED, T),
        0.25 + 0.5 * (0.564 / 0.2128) + 0.25 * (1.32 / 0.1064));
}

func testIdenticalWithThetaIsTheReciprocalNrcMatchProbability() {
    # NRC II 4.2 is P(second copy | first) by the same sampling formula, so the
    # IDENTICAL LR must be its reciprocal - two independent code paths.
    def fdb as forensics.FrequencyDb init fixture();
    for (def g in [gt("10", "10"), gt("10", "11"), gt("12", "13")]) {
        near(
            kinshipLRAtTheta($fdb, $g, $g, IDENTICAL, UNRELATED, T),
            1.0 / forensics.genotypeFrequency($fdb, $g, T));
    }
}

func testJointProbabilitiesWithThetaSumToTheMarginal() {
    # Summed over B's genotypes, every IBD term gives back P(A) under theta:
    # P(10,11) = 2 p10 (1-t) p11 = 2 * .1 * .18 = .036.
    def fdb as forensics.FrequencyDb init fixture();
    def labels as list of string init ["10", "11", "12", "13"];
    for (def rel in [UNRELATED, PARENT_CHILD, FULL_SIBLINGS, HALF_SIBLINGS]) {
        def total as float init 0.0;
        for (def i in 0..len($labels)) {
            for (def j in $i..len($labels)) {
                $total = $total +
                    jointGenotypeProbabilityTheta(
                    $fdb,
                    gt("10", "11"),
                    gt($labels[$i], $labels[$j]),
                    $rel,
                    T);
            }
        }
        near($total, 0.036);
    }
}

func testKinshipLrWithThetaIsSymmetric() {
    def fdb as forensics.FrequencyDb init fixture();
    near(
        kinshipLRAtTheta($fdb, gt("10", "10"), gt("10", "11"), FULL_SIBLINGS, UNRELATED, T),
        kinshipLRAtTheta($fdb, gt("10", "11"), gt("10", "10"), FULL_SIBLINGS, UNRELATED, T));
    near(
        kinshipLRAtTheta($fdb, gt("10", "11"), gt("11", "12"), HALF_SIBLINGS, UNRELATED, T),
        kinshipLRAtTheta($fdb, gt("11", "12"), gt("10", "11"), HALF_SIBLINGS, UNRELATED, T));
}

func testThetaLowersTheLrForASharedRareAllele() {
    # A sanity check alongside the exact values above, not instead of them.
    def fdb as forensics.FrequencyDb init fixture();
    def prev as float init kinshipLRAtTheta(
        $fdb,
        gt("10", "10"),
        gt("10", "11"),
        PARENT_CHILD,
        UNRELATED,
        0.0);
    for (def t in [0.01, 0.03, 0.1]) {
        def lr as float init kinshipLRAtTheta(
            $fdb,
            gt("10", "10"),
            gt("10", "11"),
            PARENT_CHILD,
            UNRELATED,
            $t);
        testing.assertTrue($lr < $prev);
        $prev = $lr;
    }
}

func testKinshipLrWithThetaAcrossProfiles() {
    def a as forensics.Profile init forensics.profile(
        "A",
        [gt("10", "10"), forensics.genotype("M", "10", "11")]);
    def b as forensics.Profile init forensics.profile(
        "B",
        [gt("10", "10"), forensics.genotype("M", "10", "11")]);
    def lrL as float init 0.25 + 0.5 * (1.2 / 0.39) + 0.25 * (1.32 / 0.1131);
    def lrM as float init 0.25 + 0.5 * (0.564 / 0.2128) + 0.25 * (1.32 / 0.1064);
    near(kinshipLRTheta($a, $b, FULL_SIBLINGS, twoLocus(), T), $lrL * $lrM);
    # Against half siblings: LR(FS) / LR(HS), LR(HS) = 1/2 + 1/2 LR(PC).
    near(
        kinshipLRAgainstTheta($a, $b, FULL_SIBLINGS, HALF_SIBLINGS, twoLocus(), T),
        $lrL / (0.5 + 0.5 * (1.2 / 0.39)) * $lrM / (0.5 + 0.5 * (0.564 / 0.2128)));
}

func testTrioPaternityIndexWithTheta() {
    # Mother 10/12, child 10/11, father 11/13, NO_MUTATION.
    #   X = T(M,10) T(F,11) = .5 * .5 = .25
    #   Y = T(M,10) D(11 | 10,12,11,13) = .5 * (.1 + .18) / 1.3 = .5 * .28 / 1.3
    #   PI = .5 * 1.3 / .28 = .65 / .28 = 2.3214...           (theta 0: 2.5)
    near(
        paternityIndexAtTheta(
            fixture(),
            gt("10", "12"),
            gt("10", "11"),
            gt("11", "13"),
            NO_MUTATION,
            T),
        0.65 / 0.28);
    # All three 10/10: PI = (1 + 3t) / (4t + (1-t) p10) = 1.3 / .49 = 2.6530...
    # (Buckleton et al. 2016, ch. 9, the mother aa / child aa / father aa row.)
    near(
        paternityIndexAtTheta(
            fixture(),
            gt("10", "10"),
            gt("10", "10"),
            gt("10", "10"),
            NO_MUTATION,
            T),
        1.3 / 0.49);
}

func testTrioPaternityIndexWithThetaAndMutation() {
    # Mother 10/10, child 10/12, father 11/11: the 12 needs a one-step change.
    #   X = 1 * (.002 * .9 / 2) = .0009
    #   Y = 1 * D(12 | 10,10,11,11) = .27 / 1.3
    #   PI = .0009 * 1.3 / .27 = .0043333...                  (theta 0: .003)
    def m as forensics.Profile init forensics.profile("M", [gt("10", "10")]);
    def c as forensics.Profile init forensics.profile("C", [gt("10", "12")]);
    def f as forensics.Profile init forensics.profile("F", [gt("11", "11")]);
    near(paternityIndexWithTheta($m, $c, $f, fixture(), DEFAULT_MUTATION, T), 0.0009 * 1.3 / 0.27);
    near(paternityIndexTheta($m, $c, $f, fixture(), T), 0.0009 * 1.3 / 0.27);
    near(paternityIndexWithTheta($m, $c, $f, fixture(), NO_MUTATION, T), 0.0);
}

func testDuoPaternityIndexWithTheta() {
    # Without mutation the duo index is the parent/child kinship LR.
    def fdb as forensics.FrequencyDb init fixture();
    near(duoPaternityIndexAt($fdb, gt("10", "10"), gt("10", "10"), NO_MUTATION, T), 1.2 / 0.39);
    near(duoPaternityIndexAt($fdb, gt("10", "11"), gt("10", "11"), NO_MUTATION, T), 0.564 / 0.2128);
    near(
        duoPaternityIndexAt($fdb, gt("10", "12"), gt("10", "11"), NO_MUTATION, 0.0),
        kinshipLRAt($fdb, gt("10", "11"), gt("10", "12"), PARENT_CHILD, UNRELATED));
    # Child 10/12, father 11/11, DEFAULT_MUTATION: F(10) = F(12) = .0009.
    #   X = F(10) D(12 | 11,11) + F(12) D(10 | 11,11) = .0009 * (.27 + .09) / 1.1
    #   Y = 2 D(10 | 11,11) D(12 | 11,11,10) = 2 * (.09 / 1.1) * (.27 / 1.2)
    #   PI = .0009 * .36 * 1.2 / (2 * .09 * .27) = .0003888 / .0486 = .008
    #   (theta 0: .0009 * .4 / (2 * .1 * .3) = .006)
    near(duoPaternityIndexAt($fdb, gt("10", "12"), gt("11", "11"), DEFAULT_MUTATION, T), 0.008);
    near(duoPaternityIndexAt($fdb, gt("10", "12"), gt("11", "11"), DEFAULT_MUTATION, 0.0), 0.006);
    def c as forensics.Profile init forensics.profile(
        "C",
        [gt("10", "12"), forensics.genotype("M", "10", "10")]);
    def f as forensics.Profile init forensics.profile(
        "F",
        [gt("11", "11"), forensics.genotype("M", "10", "10")]);
    near(
        duoPaternityIndex($c, $f, twoLocus(), DEFAULT_MUTATION, T),
        0.008 *
            duoPaternityIndexAt(
            twoLocus(),
            forensics.genotype("M", "10", "10"),
            forensics.genotype("M", "10", "10"),
            DEFAULT_MUTATION,
            T));
}

func thetaOfOneJoint() {
    jointGenotypeProbabilityTheta(fixture(), gt("10", "10"), gt("10", "10"), PARENT_CHILD, 1.0);
}

func negativeThetaKinship() {
    def a as forensics.Profile init forensics.profile("A", [gt("10", "10")]);
    kinshipLRTheta($a, $a, PARENT_CHILD, fixture(), -0.01);
}

func negativeThetaKinshipAt() {
    kinshipLRAtTheta(fixture(), gt("10", "10"), gt("10", "10"), PARENT_CHILD, UNRELATED, -0.01);
}

func thetaOfOneTrio() {
    paternityIndexTheta(mother(), child(), father(), fixture(), 1.0);
}

func thetaOfOneTrioAt() {
    paternityIndexAtTheta(
        fixture(),
        gt("10", "11"),
        gt("10", "12"),
        gt("12", "13"),
        NO_MUTATION,
        1.0);
}

func thetaOfOneDuo() {
    duoPaternityIndex(child(), father(), fixture(), NO_MUTATION, 1.0);
}

func mismatchedLociTheta() {
    jointGenotypeProbabilityTheta(
        fixture(),
        gt("10", "11"),
        forensics.genotype("M", "10", "11"),
        PARENT_CHILD,
        T);
}

func duoAtMixedLoci() {
    duoPaternityIndexAt(
        twoLocus(),
        gt("10", "11"),
        forensics.genotype("M", "10", "11"),
        NO_MUTATION,
        T);
}

func trioAtMixedLociTheta() {
    paternityIndexAtTheta(
        fixture(),
        gt("10", "11"),
        gt("10", "12"),
        forensics.genotype("M", "12", "13"),
        NO_MUTATION,
        T);
}

func maternalExclusionTheta() {
    def m as forensics.Profile init forensics.profile("M", [gt("13", "13")]);
    paternityIndexWithTheta($m, child(), father(), fixture(), NO_MUTATION, T);
}

func zeroDenominatorTheta() {
    kinshipLRAtTheta(fixture(), gt("10", "11"), gt("12", "13"), UNRELATED, PARENT_CHILD, T);
}

func duoWithNoSharedLocus() {
    def f as forensics.Profile init forensics.profile("F", [forensics.genotype("M", "12", "13")]);
    duoPaternityIndex(child(), $f, twoLocus(), NO_MUTATION, T);
}

func kinshipThetaWithNoSharedLocus() {
    def a as forensics.Profile init forensics.profile("A", [gt("10", "10")]);
    def b as forensics.Profile init forensics.profile("B", [forensics.genotype("M", "10", "11")]);
    kinshipLRAgainstTheta($a, $b, FULL_SIBLINGS, UNRELATED, twoLocus(), T);
}

func trioThetaWithNoSharedLocus() {
    def f as forensics.Profile init forensics.profile("F", [forensics.genotype("M", "12", "13")]);
    paternityIndexWithTheta(mother(), child(), $f, twoLocus(), NO_MUTATION, T);
}

func paternityThetaOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    paternityIndexWithTheta($t[0], $t[1], $t[2], fixture(), DEFAULT_MUTATION, T);
}

func kinshipThetaOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    kinshipLRTheta($t[1], $t[2], PARENT_CHILD, fixture(), T);
}

func duoOverAnUncoveredLocus() {
    def t as list of forensics.Profile init trioAtL("M", "C", "F");
    duoPaternityIndex($t[1], $t[2], fixture(), DEFAULT_MUTATION, T);
}

func testThetaVariantsValidateTheirInput() {
    testing.assertThrows("thetaOfOneJoint", "kinship");
    testing.assertThrows("negativeThetaKinship", "kinship");
    testing.assertThrows("negativeThetaKinshipAt", "kinship");
    testing.assertThrows("thetaOfOneTrio", "kinship");
    testing.assertThrows("thetaOfOneTrioAt", "kinship");
    testing.assertThrows("thetaOfOneDuo", "kinship");
    testing.assertThrows("mismatchedLociTheta", "kinship");
    testing.assertThrows("duoAtMixedLoci", "kinship");
    testing.assertThrows("trioAtMixedLociTheta", "kinship");
    testing.assertThrows("maternalExclusionTheta", "kinship");
    testing.assertThrows("zeroDenominatorTheta", "kinship");
    testing.assertThrows("duoWithNoSharedLocus", "kinship");
    testing.assertThrows("kinshipThetaWithNoSharedLocus", "kinship");
    testing.assertThrows("trioThetaWithNoSharedLocus", "kinship");
}

func testThetaVariantsNeverDropAnUncoveredLocus() {
    testing.assertThrows("paternityThetaOverAnUncoveredLocus", "kinship");
    testing.assertThrows("kinshipThetaOverAnUncoveredLocus", "kinship");
    testing.assertThrows("duoOverAnUncoveredLocus", "kinship");
}
