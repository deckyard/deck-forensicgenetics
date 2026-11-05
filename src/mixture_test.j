# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for mixture.j. Run with `jennifer test mixture_test.j`.
#
# The table is the kinship / pedigree fixture: loci "L" and "M", alleles
# 10/11/12/13 at .1/.2/.3/.4, N = 10000, so the 5/2N floor (0.00025) never
# bites. Every expectation is either a hand calculation written out beside it
# or an identity with an independent function of the core
# (`matchLikelihoodRatio`, NRC II 4.2). None is captured from this code.
#
# The d = 0, c = 0 closed forms used below: with no drop-out and no drop-in the
# evidence must be exactly the union of the contributors' alleles, so for u
# unknowns and no knowns, by inclusion-exclusion over the 2u allele draws,
#
#   P(union = {a})       = pa^(2u)
#   P(union = {a, b})    = (pa+pb)^(2u) - pa^(2u) - pb^(2u)
#   P(union = {a, b, c}) = (pa+pb+pc)^(2u) - Σ pairs (pi+pj)^(2u) + Σ pi^(2u)
use testing;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

# nearRel compares probabilities far below 1 relative to their own size.
func nearRel(actual as float, expected as float) {
    testing.assertTrue(math.abs($actual - $expected) <= TOL * math.abs($expected));
}

func fixture() {
    def freqs as map of string to float init {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4};
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("test"),
        "L",
        10000,
        $freqs);
    return forensics.withLocus($out, "M", 10000, $freqs);
}

func person(id as string, a as string, b as string) {
    return forensics.profile($id, [forensics.genotype("L", $a, $b)]);
}

func stain(observed as list of string) {
    return forensics.mixture("E", {"L": $observed});
}

func unknownsOnly(n as int, d as float) {
    return hypothesis("Hd", [], [], $n, $d);
}

# --- hypotheses -------------------------------------------------------------

func testHypothesisKeepsItsFields() {
    def s as forensics.Profile init person("S", "10", "11");
    def h as Hypothesis init hypothesis("Hp", [$s], [0.1], 2, 0.2);
    testing.assertEqual(len($h.known), 1);
    testing.assertEqual($h.unknowns, 2);
    near($h.knownDropout[0], 0.1);
    near($h.unknownDropout, 0.2);
}

func dropoutLengthMismatch() {
    hypothesis("H", [person("S", "10", "11")], [], 1, 0.0);
}

func tooManyUnknowns() {
    unknownsOnly(MAX_UNKNOWNS + 1, 0.0);
}

func negativeUnknowns() {
    unknownsOnly(-1, 0.0);
}

func noContributor() {
    unknownsOnly(0, 0.0);
}

func dropoutAboveOne() {
    unknownsOnly(1, 1.5);
}

func knownDropoutNegative() {
    hypothesis("H", [person("S", "10", "11")], [-0.1], 0, 0.0);
}

func samePersonTwice() {
    def s as forensics.Profile init person("S", "10", "11");
    hypothesis("H", [$s, $s], [0.0, 0.0], 0, 0.0);
}

func testHypothesisIsValidated() {
    testing.assertThrows("dropoutLengthMismatch", "mixture");
    testing.assertThrows("tooManyUnknowns", "mixture");
    testing.assertThrows("negativeUnknowns", "mixture");
    testing.assertThrows("noContributor", "mixture");
    testing.assertThrows("dropoutAboveOne", "mixture");
    testing.assertThrows("knownDropoutNegative", "mixture");
    testing.assertThrows("samePersonTwice", "mixture");
}

func testConditioningProfilesUnionBothHypotheses() {
    def s as forensics.Profile init person("S", "10", "11");
    def v as forensics.Profile init person("V", "12", "12");
    def hp as Hypothesis init hypothesis("Hp", [$s, $v], [0.0, 0.0], 0, 0.0);
    def hd as Hypothesis init hypothesis("Hd", [$v], [0.0], 1, 0.0);
    def typed as list of forensics.Profile init conditioningProfiles($hp, $hd);
    testing.assertEqual(len($typed), 2);
    testing.assertEqual($typed[0].id, "S");
    testing.assertEqual($typed[1].id, "V");
    # a person named only under Hd is typed too, after Hp's
    def alt as forensics.Profile init person("A", "13", "13");
    def hd2 as Hypothesis init hypothesis("Hd", [$alt], [0.0], 1, 0.0);
    def both as list of forensics.Profile init conditioningProfiles($hp, $hd2);
    testing.assertEqual(len($both), 3);
    testing.assertEqual($both[2].id, "A");
}

func conflictingProfiles() {
    def hp as Hypothesis init hypothesis("Hp", [person("S", "10", "11")], [0.0], 0, 0.0);
    def hd as Hypothesis init hypothesis("Hd", [person("S", "10", "12")], [0.0], 1, 0.0);
    conditioningProfiles($hp, $hd);
}

func testConditioningProfilesRejectsAConflict() {
    testing.assertThrows("conflictingProfiles", "mixture");
}

# --- mixture helpers ----------------------------------------------------------

func testMinimumContributorsIsHalfTheMostAllelesRoundedUp() {
    # L shows 3 alleles, M shows 5 (a repeat collapses): ceil(5 / 2) = 3.
    def m as forensics.Mixture init forensics.mixture(
        "E",
        {"L": ["10", "11", "12"], "M": ["10", "11", "12", "13", "14", "10"]});
    testing.assertEqual(minimumContributors($m), 3);
    # 4 alleles: ceil(4 / 2) = 2; 1 allele: 1.
    testing.assertEqual(minimumContributors(stain(["10", "11", "12", "13"])), 2);
    testing.assertEqual(minimumContributors(stain(["10"])), 1);
}

func minimumOfNothing() {
    minimumContributors(forensics.Mixture{id: "E", loci: {}});
}

func testMinimumContributorsRejectsAnEmptyMixture() {
    testing.assertThrows("minimumOfNothing", "mixture");
}

func testRestrictMixtureKeepsTheNamedLoci() {
    def m as forensics.Mixture init forensics.mixture("E", {"L": ["10"], "M": ["11"]});
    def r as forensics.Mixture init restrictMixture($m, ["M"]);
    testing.assertEqual(forensics.mixtureLoci($r), ["M"]);
    testing.assertEqual($r.loci["M"], ["11"]);
}

func restrictToMissingLocus() {
    restrictMixture(stain(["10"]), ["Q"]);
}

func testRestrictMixtureRejectsAMissingLocus() {
    testing.assertThrows("restrictToMissingLocus", "mixture");
}

# --- limiting cases: d = 0, c = 0 ----------------------------------------------

# singleSourceLR is the mixture LR of "the suspect alone" against "one unknown"
# with no drop-out and no drop-in.
func singleSourceLR(s as forensics.Profile, m as forensics.Mixture, theta as float) {
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.0], 0, 0.0);
    def hd as Hypothesis init unknownsOnly(1, 0.0);
    return mixtureLR(fixture(), $m, $hp, $hd, 0.0, $theta);
}

func testSingleSourceEqualsMatchLikelihoodRatio() {
    # Two loci, a heterozygote at L and a homozygote at M; the stain shows
    # exactly the suspect's alleles. The mixture LR must be 1/RMP, computed by
    # the core's NRC II 4.2 genotypeFrequency, at every theta.
    def s as forensics.Profile init forensics.profile(
        "S",
        [forensics.genotype("L", "10", "12"), forensics.genotype("M", "13", "13")]);
    def m as forensics.Mixture init forensics.mixture("E", {"L": ["10", "12"], "M": ["13"]});
    for (def theta in [0.0, 0.01, 0.03]) {
        nearRel(
            singleSourceLR($s, $m, $theta),
            forensics.matchLikelihoodRatio(fixture(), $s, $theta));
    }
    # theta = 0 by hand: 1 / (2 · .1 · .3 · .4²) = 1 / .0096 = 104.1666...
    nearRel(singleSourceLR($s, $m, 0.0), 1.0 / 0.0096);
}

func testTwoUnknownsShowingTwoAlleles() {
    # (pa+pb)^4 - pa^4 - pb^4 = .3^4 - .1^4 - .2^4 = .0081 - .0001 - .0016 = .0064
    def e as float init evidenceProbability(
        fixture(),
        stain(["10", "11"]),
        "L",
        unknownsOnly(2, 0.0),
        [],
        0.0,
        0.0);
    nearRel($e, 0.0064);
}

func testTwoUnknownsShowingThreeAlleles() {
    # (.1+.2+.3)^4                     = .6^4             = .1296
    # - (.3^4 + .4^4 + .5^4)           = .0081+.0256+.0625 = .0962
    # + (.1^4 + .2^4 + .3^4)           = .0001+.0016+.0081 = .0098
    # = .1296 - .0962 + .0098 = .0432
    def e as float init evidenceProbability(
        fixture(),
        stain(["10", "11", "12"]),
        "L",
        unknownsOnly(2, 0.0),
        [],
        0.0,
        0.0);
    nearRel($e, 0.0432);
}

func testThreeUnknowns() {
    # One allele: .1^6 = 1e-6. Two alleles: .3^6 - .1^6 - .2^6
    # = .000729 - .000001 - .000064 = .000664.
    def fdb as forensics.FrequencyDb init fixture();
    def h as Hypothesis init unknownsOnly(3, 0.0);
    nearRel(evidenceProbability($fdb, stain(["10"]), "L", $h, [], 0.0, 0.0), 0.000001);
    nearRel(evidenceProbability($fdb, stain(["10", "11"]), "L", $h, [], 0.0, 0.0), 0.000664);
}

func testSuspectPlusUnknownAgainstTwoUnknowns() {
    # Hp: S (10/11) and one unknown, who may carry only 10 and 11:
    #   (p10 + p11)^2 = .3^2 = .09
    # Hd: two unknowns, union exactly {10, 11}: .0064 (above).
    # LR = .09 / .0064 = 14.0625
    def s as forensics.Profile init person("S", "10", "11");
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.0], 1, 0.0);
    def hd as Hypothesis init unknownsOnly(2, 0.0);
    near(locusLR(fixture(), stain(["10", "11"]), "L", $hp, $hd, 0.0, 0.0), 14.0625);
}

func testVictimConditionedOnBothSides() {
    # Hp: V (12/12) and S (10/11); Hd: V and one unknown. Evidence {10, 11, 12}.
    # Hp: everything explained, 1. Hd: the unknown must carry 10 and 11
    # (V covers 12, so 12 is optional): genotype 10/11 only, 2 · .1 · .2 = .04.
    # LR = 1 / .04 = 25.
    def s as forensics.Profile init person("S", "10", "11");
    def v as forensics.Profile init person("V", "12", "12");
    def hp as Hypothesis init hypothesis("Hp", [$v, $s], [0.0, 0.0], 0, 0.0);
    def hd as Hypothesis init hypothesis("Hd", [$v], [0.0], 1, 0.0);
    near(locusLR(fixture(), stain(["10", "11", "12"]), "L", $hp, $hd, 0.0, 0.0), 25.0);
}

func testThetaSamplingFormulaForTwoUnknowns() {
    # Two unknowns, nobody typed, evidence {10}: both 10/10, four draws of 10.
    # Sampling formula (n_a θ + (1-θ) p) / (1 + (n-1) θ) after n alleles:
    #   n=0: (1-θ)p / (1-θ) = p
    #   n=1: (θ + (1-θ)p) / 1
    #   n=2: (2θ + (1-θ)p) / (1+θ)
    #   n=3: (3θ + (1-θ)p) / (1+2θ)
    # θ = .01, p = .1: .1 · .109 · (.119/1.01) · (.129/1.02)
    def t as float init 0.01;
    def p as float init 0.1;
    def expected as float init $p *
        ($t + (1.0 - $t) * $p) *
        ((2.0 * $t + (1.0 - $t) * $p) / (1.0 + $t)) *
        ((3.0 * $t + (1.0 - $t) * $p) / (1.0 + 2.0 * $t));
    nearRel($expected, 0.1 * 0.109 * (0.119 / 1.01) * (0.129 / 1.02));
    def e as float init evidenceProbability(
        fixture(),
        stain(["10"]),
        "L",
        unknownsOnly(2, 0.0),
        [],
        0.0,
        $t);
    nearRel($e, $expected);
}

# --- drop-out and drop-in -------------------------------------------------------

func testDropoutOfTheSuspectsSecondAllele() {
    # Evidence {10}; S = 10/11; d = .2 for everyone; c = 0.
    # Hp: 10 seen (1-d), 11 dropped (d): .8 · .2 = .16
    # Hd: one unknown, who must carry 10 (no drop-in):
    #   10/10: p² (1 - d²)          = .01 · .96        = .0096
    #   10/x : 2 p (1-p) (1-d) d    = 2·.1·.9·.8·.2    = .0288
    #   total .0384  (x summed over 11, 12, 13: checks the lumped class)
    # LR = .16 / .0384 = 4.1666...
    def s as forensics.Profile init person("S", "10", "11");
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.2], 0, 0.0);
    def hd as Hypothesis init unknownsOnly(1, 0.2);
    def typed as list of forensics.Profile init conditioningProfiles($hp, $hd);
    def fdb as forensics.FrequencyDb init fixture();
    def m as forensics.Mixture init stain(["10"]);
    nearRel(evidenceProbability($fdb, $m, "L", $hp, $typed, 0.0, 0.0), 0.16);
    nearRel(evidenceProbability($fdb, $m, "L", $hd, $typed, 0.0, 0.0), 0.0384);
    near(locusLR($fdb, $m, "L", $hp, $hd, 0.0, 0.0), 0.16 / 0.0384);
}

func testDropIn() {
    # As above with c = .05.
    # Hp: .16 · (1 - c) = .16 · .95 = .152
    # Hd: genotypes carrying 10 as before, times (1 - c): .0384 · .95 = .03648
    #     genotypes without 10 (prob (1-p10)² = .81): every allele drops (d² for
    #     a heterozygote d·d, for a homozygote d²) and 10 drops in (c · p10):
    #     .81 · .04 · .05 · .1 = .000162
    #     total .036642
    # LR = .152 / .036642
    def s as forensics.Profile init person("S", "10", "11");
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.2], 0, 0.0);
    def hd as Hypothesis init unknownsOnly(1, 0.2);
    def typed as list of forensics.Profile init conditioningProfiles($hp, $hd);
    def fdb as forensics.FrequencyDb init fixture();
    def m as forensics.Mixture init stain(["10"]);
    nearRel(evidenceProbability($fdb, $m, "L", $hp, $typed, 0.05, 0.0), 0.152);
    nearRel(evidenceProbability($fdb, $m, "L", $hd, $typed, 0.05, 0.0), 0.036642);
    near(locusLR($fdb, $m, "L", $hp, $hd, 0.05, 0.0), 0.152 / 0.036642);
}

func testHomozygoteDropsOutWithDSquared() {
    # S = 10/10 with d = .3, evidence empty (a struct literal: the builder
    # rejects an empty locus, the model does not): d² · (1 - c) = .09.
    def s as forensics.Profile init person("S", "10", "10");
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.3], 0, 0.0);
    def m as forensics.Mixture init forensics.Mixture{id: "E", loci: {"L": []}};
    nearRel(evidenceProbability(fixture(), $m, "L", $hp, [$s], 0.0, 0.0), 0.09);
}

func testCompleteDropoutIsUninformative() {
    # Nothing seen, d = .5 for one contributor either way: every genotype,
    # homozygote or not, drops out with d² = .25, so P = .25 under both and
    # the LR is 1.
    def s as forensics.Profile init person("S", "10", "11");
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.5], 0, 0.0);
    def hd as Hypothesis init unknownsOnly(1, 0.5);
    def m as forensics.Mixture init forensics.Mixture{id: "E", loci: {"L": []}};
    near(locusLR(fixture(), $m, "L", $hp, $hd, 0.0, 0.0), 1.0);
}

func testSharedAlleleDropsOnlyIfEveryCopyDoes() {
    # Two knowns, A = 10/11 (d = .5) and B = 10/12 (d = .4), no unknowns,
    # evidence {11, 12}: 10 must drop from both (.5 · .4 = .2), 11 stays (.5),
    # 12 stays (.6), no drop-in (1 - c = 1): .2 · .5 · .6 = .06
    def a as forensics.Profile init person("A", "10", "11");
    def b as forensics.Profile init person("B", "10", "12");
    def h as Hypothesis init hypothesis("H", [$a, $b], [0.5, 0.4], 0, 0.0);
    nearRel(evidenceProbability(fixture(), stain(["11", "12"]), "L", $h, [$a, $b], 0.0, 0.0), 0.06);
}

# --- across loci ----------------------------------------------------------------

func testLocusLRsMultiplyAndLogsAdd() {
    # L as testSuspectPlusUnknownAgainstTwoUnknowns (14.0625); M: S = 12/13,
    # stain {12, 13}: Hp (.3+.4)² = .49, Hd .7^4 - .3^4 - .4^4
    # = .2401 - .0081 - .0256 = .2064, LR .49 / .2064.
    def s as forensics.Profile init forensics.profile(
        "S",
        [forensics.genotype("L", "10", "11"), forensics.genotype("M", "12", "13")]);
    def m as forensics.Mixture init forensics.mixture("E", {"L": ["10", "11"], "M": ["12", "13"]});
    def hp as Hypothesis init hypothesis("Hp", [$s], [0.0], 1, 0.0);
    def hd as Hypothesis init unknownsOnly(2, 0.0);
    def fdb as forensics.FrequencyDb init fixture();
    def per as map of string to float init locusLRs($fdb, $m, $hp, $hd, 0.0, 0.0);
    near($per["L"], 14.0625);
    near($per["M"], 0.49 / 0.2064);
    near(mixtureLR($fdb, $m, $hp, $hd, 0.0, 0.0), 14.0625 * 0.49 / 0.2064);
    near(log10MixtureLR($fdb, $m, $hp, $hd, 0.0, 0.0), math.log10(14.0625 * 0.49 / 0.2064));
}

func testExclusionGivesZero() {
    # S = 12/13 cannot explain {10, 11} with no drop-out or drop-in.
    def hp as Hypothesis init hypothesis("Hp", [person("S", "12", "13")], [0.0], 0, 0.0);
    def hd as Hypothesis init unknownsOnly(1, 0.0);
    near(mixtureLR(fixture(), stain(["10", "11"]), $hp, $hd, 0.0, 0.0), 0.0);
}

func logOfExclusion() {
    def hp as Hypothesis init hypothesis("Hp", [person("S", "12", "13")], [0.0], 0, 0.0);
    log10MixtureLR(fixture(), stain(["10", "11"]), $hp, unknownsOnly(1, 0.0), 0.0, 0.0);
}

func impossibleDefence() {
    # one unknown, no drop-in, three alleles: P(E | Hd) = 0
    def hp as Hypothesis init hypothesis("Hp", [person("S", "10", "11")], [0.0], 1, 0.0);
    locusLR(fixture(), stain(["10", "11", "12"]), "L", $hp, unknownsOnly(1, 0.0), 0.0, 0.0);
}

func testZeroProbabilitiesThrowWhereTheRatioIsUndefined() {
    testing.assertThrows("logOfExclusion", "mixture");
    testing.assertThrows("impossibleDefence", "mixture");
}

# --- invalid input ----------------------------------------------------------------

func badTheta() {
    evidenceProbability(fixture(), stain(["10"]), "L", unknownsOnly(1, 0.0), [], 0.0, 1.0);
}

func badDropIn() {
    evidenceProbability(fixture(), stain(["10"]), "L", unknownsOnly(1, 0.0), [], 1.0, 0.0);
}

func locusNotInMixture() {
    evidenceProbability(fixture(), stain(["10"]), "M", unknownsOnly(1, 0.0), [], 0.0, 0.0);
}

func locusNotInTable() {
    def m as forensics.Mixture init forensics.mixture("E", {"Q": ["10"]});
    evidenceProbability(fixture(), $m, "Q", unknownsOnly(1, 0.0), [], 0.0, 0.0);
}

func knownUntypedAtLocus() {
    def m as forensics.Mixture init forensics.mixture("E", {"M": ["10"]});
    def hp as Hypothesis init hypothesis("Hp", [person("S", "10", "10")], [0.0], 0, 0.0);
    mixtureLR(fixture(), $m, $hp, unknownsOnly(1, 0.0), 0.0, 0.0);
}

func knownNotTyped() {
    def hp as Hypothesis init hypothesis("Hp", [person("S", "10", "10")], [0.0], 0, 0.0);
    evidenceProbability(fixture(), stain(["10"]), "L", $hp, [], 0.0, 0.0);
}

func tooManyAlleles() {
    def freqs as map of string to float init {};
    def seen as list of string init [];
    for (def i in 0..MAX_LOCUS_ALLELES + 1) {
        $freqs["{$i + 1}"] = 0.05;
        $seen[] = "{$i + 1}";
    }
    def fdb as forensics.FrequencyDb init forensics.withLocus(forensics.db("t"), "L", 1000, $freqs);
    evidenceProbability($fdb, stain($seen), "L", unknownsOnly(1, 0.0), [], 0.0, 0.0);
}

func flooredPastOne() {
    # N = 2: the 5/2N floor is 1.25, capped at 1 per allele; two alleles sum to 2.
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("t"),
        "L",
        2,
        {"10": 0.5, "11": 0.5});
    evidenceProbability($fdb, stain(["10", "11"]), "L", unknownsOnly(1, 0.0), [], 0.0, 0.0);
}

func emptyMixture() {
    mixtureLR(
        fixture(),
        forensics.Mixture{id: "E", loci: {}},
        unknownsOnly(1, 0.0),
        unknownsOnly(1, 0.0),
        0.0,
        0.0);
}

func literalWithBadHypothesis() {
    def h as Hypothesis init Hypothesis{
        name: "H",
        known: [],
        knownDropout: [],
        unknowns: 9,
        unknownDropout: 0.0
    };
    evidenceProbability(fixture(), stain(["10"]), "L", $h, [], 0.0, 0.0);
}

func testInvalidInputThrowsMixtureErrors() {
    testing.assertThrows("badTheta", "mixture");
    testing.assertThrows("badDropIn", "mixture");
    testing.assertThrows("locusNotInMixture", "mixture");
    testing.assertThrows("locusNotInTable", "mixture");
    testing.assertThrows("knownUntypedAtLocus", "mixture");
    testing.assertThrows("knownNotTyped", "mixture");
    testing.assertThrows("tooManyAlleles", "mixture");
    testing.assertThrows("flooredPastOne", "mixture");
    testing.assertThrows("emptyMixture", "mixture");
    testing.assertThrows("literalWithBadHypothesis", "mixture");
}

# --- internals ------------------------------------------------------------------

func testOrderingsCountsDistinctPermutations() {
    testing.assertEqual(orderings([]), 1);
    testing.assertEqual(orderings([4]), 1);
    testing.assertEqual(orderings([1, 1, 1]), 1);
    testing.assertEqual(orderings([1, 1, 2]), 3);
    testing.assertEqual(orderings([1, 2, 3]), 6);
}

func testNextComboVisitsEveryMultiset() {
    # multisets of size 2 from 3 genotypes: C(3 + 2 - 1, 2) = 6
    def combo as list of int init [0, 0];
    def visited as int init 0;
    while (len($combo) > 0) {
        $visited = $visited + 1;
        $combo = nextCombo($combo, 3);
    }
    testing.assertEqual($visited, 6);
}
