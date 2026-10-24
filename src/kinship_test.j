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

func testZeroRateNeverCountsSteps() {
    # With no mutation a changed allele has probability 0 before any step
    # count is attempted, so non-numeric designations do not throw.
    near(mutationProbability(NO_MUTATION, "AMEL", "X", "Y"), 0.0);
    near(mutationProbability(NO_MUTATION, "AMEL", "X", "X"), 1.0);
}
