# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for forensicgenetics.j. Run with `jennifer test forensicgenetics_test.j`.
#
# Expected values are hand-computed from the published NRC II formulae rather
# than captured from this implementation, so a formula regression fails here
# instead of being blessed.
use testing;
use math;

def const TOL as float init 0.000000001;

# near asserts two floats agree to within TOL.
func near(actual as float, expected as float) {
    testing.assertTrue(math.abs($actual - $expected) < TOL);
}

# fixture builds a two-locus reference table:
#   D3S1358, N = 100 (floor 5/200 = 0.025): 15 -> .25, 16 -> .20, 17 -> .30, 18 -> .25
#   TH01,    N =  50 (floor 5/100 = 0.050): 9.3 -> .30, 6 -> .20, 7 -> .50
func fixture() {
    def out as FrequencyDb init db("test");
    $out = withLocus($out, "D3S1358", 100, {"15": 0.25, "16": 0.20, "17": 0.30, "18": 0.25});
    $out = withLocus($out, "TH01", 50, {"9.3": 0.30, "6": 0.20, "7": 0.50});
    return $out;
}

# suspect is the reference profile the probability tests score.
func suspect() {
    return profile("S1", [genotype("D3S1358", "15", "16"), genotype("TH01", "9.3", "6")]);
}

# --- allele labels ----------------------------------------------------------

func testAlleleValueReadsMicrovariants() {
    near(alleleValue("15"), 15.0);
    near(alleleValue("9.3"), 9.3);
    near(alleleValue("15.2"), 15.2);
}

func testIsNumericAlleleRejectsDesignations() {
    testing.assertTrue(isNumericAllele("15"));
    testing.assertTrue(isNumericAllele("9.3"));
    testing.assertFalse(isNumericAllele("X"));
    testing.assertFalse(isNumericAllele("OL"));
    testing.assertFalse(isNumericAllele(""));
    testing.assertFalse(isNumericAllele("1.2.3"));
    testing.assertFalse(isNumericAllele("."));
}

func nonNumericAlleleValue() {
    alleleValue("X");
}

func testAlleleValueThrowsOnDesignation() {
    testing.assertThrows("nonNumericAlleleValue", "forensicgenetics");
}

# --- genotypes --------------------------------------------------------------

func testGenotypeOrdersItsPairNumerically() {
    # "9.3" must sort after "6" by repeat count, not by string order.
    def g as Genotype init genotype("TH01", "9.3", "6");
    testing.assertEqual($g.a, "6");
    testing.assertEqual($g.b, "9.3");
    # Either writing order yields the same value, so genotypes compare equal.
    testing.assertEqual(genotype("TH01", "6", "9.3"), $g);
}

func testGenotypeOrdersNonNumericLexicographically() {
    def g as Genotype init genotype("AMEL", "Y", "X");
    testing.assertEqual($g.a, "X");
    testing.assertEqual($g.b, "Y");
}

func testHomozygosityAndAlleles() {
    testing.assertTrue(isHomozygous(genotype("TH01", "6", "6")));
    testing.assertFalse(isHomozygous(genotype("TH01", "6", "7")));
    testing.assertEqual(alleles(genotype("TH01", "6", "6")), ["6"]);
    testing.assertEqual(alleles(genotype("TH01", "6", "7")), ["6", "7"]);
}

func emptyAlleleLabel() {
    genotype("TH01", "6", "");
}

func testGenotypeRejectsEmptyLabel() {
    testing.assertThrows("emptyAlleleLabel", "forensicgenetics");
}

# --- profiles ---------------------------------------------------------------

func testProfileBuildAndAccess() {
    def p as Profile init suspect();
    testing.assertEqual(profileLoci($p), ["D3S1358", "TH01"]);
    testing.assertTrue(hasGenotype($p, "TH01"));
    testing.assertFalse(hasGenotype($p, "FGA"));
    testing.assertEqual(genotypeAt($p, "TH01"), genotype("TH01", "6", "9.3"));
}

func testWithGenotypeReplacesAndCopies() {
    def p as Profile init suspect();
    def q as Profile init withGenotype($p, genotype("TH01", "7", "7"));
    testing.assertEqual(genotypeAt($q, "TH01"), genotype("TH01", "7", "7"));
    # Value semantics: the original is untouched.
    testing.assertEqual(genotypeAt($p, "TH01"), genotype("TH01", "6", "9.3"));
}

func duplicateLocusProfile() {
    profile("dup", [genotype("TH01", "6", "7"), genotype("TH01", "8", "9")]);
}

func testProfileRejectsDuplicateLocus() {
    testing.assertThrows("duplicateLocusProfile", "forensicgenetics");
}

func testCommonLociAndRestrict() {
    def fdb as FrequencyDb init fixture();
    def p as Profile init withGenotype(suspect(), genotype("FGA", "22", "23"));
    testing.assertEqual(commonLoci($fdb, $p), ["D3S1358", "TH01"]);
    def narrowed as Profile init restrict($p, commonLoci($fdb, $p));
    testing.assertEqual(profileLoci($narrowed), ["D3S1358", "TH01"]);
}

# --- frequencies and the 5/2N floor -----------------------------------------

func testMinFrequencyIsFiveOverTwoN() {
    def fdb as FrequencyDb init fixture();
    near(minFrequency($fdb, "D3S1358"), 0.025);
    near(minFrequency($fdb, "TH01"), 0.05);
}

func testUnobservedAlleleTakesTheFloor() {
    def fdb as FrequencyDb init fixture();
    near(alleleFrequency($fdb, "D3S1358", "19"), 0.025);
    # An observed allele rarer than the floor is raised to it.
    def rare as FrequencyDb init withLocus(db("rare"), "TPOX", 100, {"8": 0.001, "9": 0.999});
    near(alleleFrequency($rare, "TPOX", "8"), 0.025);
    near(alleleFrequency($rare, "TPOX", "9"), 0.999);
}

func missingLocusFrequency() {
    alleleFrequency(fixture(), "FGA", "22");
}

func testUnknownLocusThrows() {
    testing.assertThrows("missingLocusFrequency", "forensicgenetics");
}

func testSampleSizeAndLocusListing() {
    def fdb as FrequencyDb init fixture();
    testing.assertEqual(sampleSize($fdb, "TH01"), 50);
    testing.assertEqual(loci($fdb), ["D3S1358", "TH01"]);
    testing.assertTrue(hasLocus($fdb, "D3S1358"));
    testing.assertFalse(hasLocus($fdb, "FGA"));
}

# --- Hardy-Weinberg and the theta correction --------------------------------

func testHardyWeinbergAtThetaZero() {
    def fdb as FrequencyDb init fixture();
    # heterozygote 2 pa pb = 2 * .25 * .20
    near(genotypeFrequency($fdb, genotype("D3S1358", "15", "16"), THETA_NONE), 0.1);
    # homozygote p^2 = .25^2
    near(genotypeFrequency($fdb, genotype("D3S1358", "15", "15"), THETA_NONE), 0.0625);
}

func testNrc2ThetaCorrection() {
    def fdb as FrequencyDb init fixture();
    # 2 [t+(1-t)pa][t+(1-t)pb] / [(1+t)(1+2t)], t = .01, pa = .25, pb = .20
    near(
        genotypeFrequency($fdb, genotype("D3S1358", "15", "16"), THETA_GENERAL),
        0.10397980974568045);
    # [2t+(1-t)p][3t+(1-t)p] / [(1+t)(1+2t)], t = .01, p = .25
    near(
        genotypeFrequency($fdb, genotype("D3S1358", "15", "15"), THETA_GENERAL),
        0.07205518345952241);
}

func testThetaMakesTheEstimateMoreConservative() {
    def fdb as FrequencyDb init fixture();
    def g as Genotype init genotype("D3S1358", "15", "15");
    def plain as float init genotypeFrequency($fdb, $g, THETA_NONE);
    def general as float init genotypeFrequency($fdb, $g, THETA_GENERAL);
    def isolated as float init genotypeFrequency($fdb, $g, THETA_ISOLATED);
    testing.assertTrue($plain < $general);
    testing.assertTrue($general < $isolated);
}

func negativeTheta() {
    genotypeFrequency(fixture(), genotype("D3S1358", "15", "16"), -0.1);
}

func thetaOne() {
    genotypeFrequency(fixture(), genotype("D3S1358", "15", "16"), 1.0);
}

func testThetaRangeIsChecked() {
    testing.assertThrows("negativeTheta", "forensicgenetics");
    testing.assertThrows("thetaOne", "forensicgenetics");
}

# --- RMP and the single-source LR -------------------------------------------

func testRandomMatchProbabilityAtThetaZero() {
    # (2 * .25 * .20) * (2 * .30 * .20) = .1 * .12
    near(randomMatchProbability(fixture(), suspect(), THETA_NONE), 0.012);
}

func testRandomMatchProbabilityWithTheta() {
    near(randomMatchProbability(fixture(), suspect(), THETA_GENERAL), 0.012890185849582939);
}

func testLog10RmpAgreesWithTheProduct() {
    def fdb as FrequencyDb init fixture();
    near(log10RandomMatchProbability($fdb, suspect(), THETA_GENERAL), -1.8897408209812383);
    near(
        math.pow(10.0, log10RandomMatchProbability($fdb, suspect(), THETA_GENERAL)),
        randomMatchProbability($fdb, suspect(), THETA_GENERAL));
}

func testMatchLikelihoodRatioIsTheReciprocal() {
    near(matchLikelihoodRatio(fixture(), suspect(), THETA_GENERAL), 77.57840047219761);
}

func rmpWithUncoveredLocus() {
    randomMatchProbability(fixture(), profile("X", [genotype("FGA", "22", "23")]), THETA_NONE);
}

func rmpOfEmptyProfile() {
    randomMatchProbability(fixture(), profile("X", []), THETA_NONE);
}

func testRmpRefusesToDropLoci() {
    # Silently skipping an uncovered locus would inflate the RMP without saying so.
    testing.assertThrows("rmpWithUncoveredLocus", "forensicgenetics");
    testing.assertThrows("rmpOfEmptyProfile", "forensicgenetics");
}

# --- mixtures: CPI / CPE ----------------------------------------------------

func stain() {
    return mixture("M1", {"D3S1358": ["15", "16"], "TH01": ["9.3", "6"]});
}

func testInclusionAndExclusionAtOneLocus() {
    def fdb as FrequencyDb init fixture();
    # (.25 + .20)^2
    near(inclusionProbability($fdb, "D3S1358", ["15", "16"]), 0.2025);
    near(exclusionProbability($fdb, "D3S1358", ["15", "16"]), 0.7975);
}

func testInclusionIgnoresRepeatedAlleles() {
    def fdb as FrequencyDb init fixture();
    near(inclusionProbability($fdb, "D3S1358", ["15", "16", "15"]), 0.2025);
}

func testMixtureLociListing() {
    testing.assertEqual(mixtureLoci(stain()), ["D3S1358", "TH01"]);
}

func testCombinedInclusionAndExclusion() {
    def fdb as FrequencyDb init fixture();
    # .2025 * .25
    near(cpi($fdb, stain()), 0.050625);
    near(cpe($fdb, stain()), 0.949375);
}

func testInclusionProbabilityIsCappedAtOne() {
    # N = 10 puts the 5/2N floor at 0.25, so five floored alleles sum past 1.0.
    # PI must still not exceed 1.0.
    def fdb as FrequencyDb init withLocus(
        db("small"),
        "W",
        10,
        {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2});
    near(minFrequency($fdb, "W"), 0.25);
    near(inclusionProbability($fdb, "W", ["1", "2", "3", "4", "5"]), 1.0);
}

func testIsIncludedMatchesTheCpiPremise() {
    testing.assertTrue(isIncluded(suspect(), stain()));
    # An allele the stain never showed excludes the profile.
    def other as Profile init profile(
        "S2",
        [genotype("D3S1358", "17", "18"), genotype("TH01", "9.3", "6")]);
    testing.assertFalse(isIncluded($other, stain()));
    # Loci the mixture does not cover are simply skipped.
    def extra as Profile init withGenotype(suspect(), genotype("FGA", "22", "23"));
    testing.assertTrue(isIncluded($extra, stain()));
}

func emptyMixtureLocus() {
    mixture("bad", {"TH01": []});
}

func testMixtureRejectsEmptyLocus() {
    testing.assertThrows("emptyMixtureLocus", "forensicgenetics");
}

# --- regressions from the validation audit ----------------------------------

func tableSummingWellPastOne() {
    withLocus(db("broken"), "X", 100, {"1": 0.9, "2": 0.9, "3": 0.9});
}

func testAMisreadTableIsRejected() {
    # Left alone, a table like this yields genotype "probabilities" above 1.
    testing.assertThrows("tableSummingWellPastOne", "forensicgenetics");
    # A partial table summing below 1 stays legal - it is a real thing to have.
    def partial as FrequencyDb init withLocus(db("partial"), "X", 100, {"1": 0.4, "2": 0.4});
    near(alleleFrequency($partial, "X", "1"), 0.4);
}

func testGenotypeFrequencyNeverExceedsOne() {
    # N = 1 puts the floor at 2.5, so both alleles floor to 1.0 and the
    # heterozygote formula would otherwise return 2.0.
    def tiny as FrequencyDb init withLocus(db("tiny"), "Y", 1, {"1": 0.5, "2": 0.5});
    near(alleleFrequency($tiny, "Y", "1"), 1.0);
    near(genotypeFrequency($tiny, genotype("Y", "1", "2"), THETA_NONE), 1.0);
    near(randomMatchProbability($tiny, profile("P", [genotype("Y", "1", "2")]), THETA_NONE), 1.0);
}

func inclusionWithNoSharedLocus() {
    isIncluded(profile("P", [genotype("FGA", "22", "24")]), stain());
}

func testInclusionRefusesAVacuousAnswer() {
    # No shared locus is not "not excluded" - there is nothing to compare.
    testing.assertThrows("inclusionWithNoSharedLocus", "forensicgenetics");
}

# --- pooling several populations --------------------------------------------

# Two populations that disagree sharply at one locus. Equal sample sizes, so
# the pooled table sits exactly halfway between them.
func popA() {
    return withLocus(db("A"), "TH01", 1000, {"6": 0.8, "9.3": 0.2});
}

func popB() {
    return withLocus(db("B"), "TH01", 1000, {"6": 0.2, "9.3": 0.8});
}

func testPoolingAveragesEqualSizedSamples() {
    def pooled as FrequencyDb init combineFrequencies([popA(), popB()], "A+B");
    testing.assertEqual($pooled.name, "A+B");
    # 0.8*2000 + 0.2*2000 = 2000 of 4000 alleles.
    near(alleleFrequency($pooled, "TH01", "6"), 0.5);
    near(alleleFrequency($pooled, "TH01", "9.3"), 0.5);
    # The pooled N is the sum, so the 5/2N floor drops accordingly.
    testing.assertEqual(sampleSize($pooled, "TH01"), 2000);
    near(minFrequency($pooled, "TH01"), 0.00125);
}

func testPoolingWeightsBySampleSize() {
    def small as FrequencyDb init withLocus(db("B"), "TH01", 500, {"6": 0.2, "9.3": 0.8});
    def pooled as FrequencyDb init combineFrequencies([popA(), $small], "A+B");
    # 0.8*2000 + 0.2*1000 = 1800 of 3000 alleles: the larger sample dominates.
    near(alleleFrequency($pooled, "TH01", "6"), 0.6);
    testing.assertEqual(sampleSize($pooled, "TH01"), 1500);
}

func testPoolingKeepsOnlySharedLoci() {
    def a as FrequencyDb init withLocus(popA(), "FGA", 1000, {"21": 1.0});
    def b as FrequencyDb init withLocus(popB(), "VWA", 1000, {"17": 1.0});
    testing.assertEqual(loci(combineFrequencies([$a, $b], "shared")), ["TH01"]);
}

func noSharedLocus() {
    combineFrequencies(
        [
            withLocus(db("A"), "FGA", 100, {"21": 1.0}),
            withLocus(db("B"), "VWA", 100, {"17": 1.0})
        ],
        "none");
}

func testPoolingWithNoSharedLocusThrows() {
    testing.assertThrows("noSharedLocus", "forensicgenetics");
}

func noDatabasesToPool() {
    combineFrequencies([], "empty");
}

func testPoolingNothingThrows() {
    testing.assertThrows("noDatabasesToPool", "forensicgenetics");
}

# --- an uncertain population ------------------------------------------------

func homozygote() {
    return profile("S1", [genotype("TH01", "6", "6")]);
}

func testAcrossIsTheWeightedMeanOfTheSeparateAnswers() {
    def dbs as list of FrequencyDb init [popA(), popB()];
    # 0.5*0.8^2 + 0.5*0.2^2 = 0.32 + 0.02.
    near(randomMatchProbabilityAcross($dbs, [0.5, 0.5], homozygote(), THETA_NONE), 0.34);
    near(matchLikelihoodRatioAcross($dbs, [0.5, 0.5], homozygote(), THETA_NONE), 1.0 / 0.34);
}

func testAsymmetricWeightsShiftTheAnswer() {
    def dbs as list of FrequencyDb init [popA(), popB()];
    # 0.9*0.64 + 0.1*0.04 = 0.576 + 0.004.
    near(randomMatchProbabilityAcross($dbs, [0.9, 0.1], homozygote(), THETA_NONE), 0.58);
}

func testPoolingUnderstatesAnUncertainPopulation() {
    # The Wahlund effect, in one assertion. The pooled table says 0.5^2 = 0.25;
    # the law of total probability over the same two populations says 0.34.
    # Pooling would overstate this match by a third.
    def pooled as FrequencyDb init combineFrequencies([popA(), popB()], "A+B");
    near(randomMatchProbability($pooled, homozygote(), THETA_NONE), 0.25);
    near(
        randomMatchProbabilityAcross([popA(), popB()], [0.5, 0.5], homozygote(), THETA_NONE),
        0.34);
}

func weightsDoNotSumToOne() {
    randomMatchProbabilityAcross([popA(), popB()], [0.5, 0.9], homozygote(), THETA_NONE);
}

func testWeightsMustSumToOne() {
    testing.assertThrows("weightsDoNotSumToOne", "forensicgenetics");
}

func weightCountMismatch() {
    randomMatchProbabilityAcross([popA(), popB()], [1.0], homozygote(), THETA_NONE);
}

func testWeightCountMustMatchDatabaseCount() {
    testing.assertThrows("weightCountMismatch", "forensicgenetics");
}
