# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for xstr.j. Run with `jennifer test xstr_test.j`.
#
# Loci "L" and "M" both have alleles 10/11/12/13 at .1/.2/.3/.4 and
# N = 10000 (floor 0.00025, never reached), as in kinship_test.j. Every
# expectation is written out by hand from the X-IBD formulae, checked against
# an identity with existing code (kinship.j, the independent model), or taken
# from Gomes et al. (2020), Front. Genet. 11:926, Table 5.
use testing;
use math;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

func fixture() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("test"),
        "L",
        10000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    return forensics.withLocus($out, "M", 10000, {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

func model() {
    return independentModel(fixture());
}

# linked is a synthetic two-locus group {L, M} out of linkage equilibrium:
# 10-10 .3, 11-11 .3, 10-11 .1, 11-10 .1, 12-12 .2 (sum 1), floor .001.
func linked() {
    def g as LinkageGroup init linkageGroup(
        "LM",
        ["L", "M"],
        0.001,
        {"10-10": 0.3, "11-11": 0.3, "10-11": 0.1, "11-10": 0.1, "12-12": 0.2});
    return withLinkageGroup(model(), $g);
}

# equilibrium is the group {L, M} with every haplotype at the product of its
# allele frequencies: linkage equilibrium.
func equilibrium() {
    def freqs as map of string to float init {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4};
    def table as map of string to float init {};
    for (def x in $freqs) {
        for (def y in $freqs) {
            $table[haplotypeKey([$x, $y])] = $freqs[$x] * $freqs[$y];
        }
    }
    return withLinkageGroup(model(), linkageGroup("LM", ["L", "M"], 0.0, $table));
}

func male(id as string, l as string) {
    return maleProfile($id, {"L": $l});
}

func female(id as string, a as string, b as string) {
    return femaleProfile($id, [forensics.genotype("L", $a, $b)]);
}

func male2(id as string, l as string, m as string) {
    return maleProfile($id, {"L": $l, "M": $m});
}

func female2(id as string, la as string, lb as string, ma as string, mb as string) {
    return femaleProfile(
        $id,
        [forensics.genotype("L", $la, $lb), forensics.genotype("M", $ma, $mb)]);
}

# --- relationships ----------------------------------------------------------

func testCoancestryMatchesGomes2020Table5() {
    # Gomes et al. (2020), Table 5, X-chr rows: female-female, male-male,
    # female-male, male-female.
    near(coancestry(MOTHER_DAUGHTER), 0.25); # parenthood F-F
    near(coancestry(FATHER_SON), 0.0); # parenthood M-M
    near(coancestry(MOTHER_SON), 0.5); # parenthood F-M
    near(coancestry(FATHER_DAUGHTER), 0.5); # parenthood M-F
    near(coancestry(FULL_SISTERS), 0.375); # full sibship F-F, 3/8
    near(coancestry(FULL_BROTHERS), 0.5); # full sibship M-M
    near(coancestry(BROTHER_SISTER), 0.25); # full sibship M-F
    near(coancestry(PATERNAL_HALF_SISTERS), 0.25); # paternal half-sibship F-F
    near(coancestry(MATERNAL_HALF_SISTERS), 0.125); # maternal half-sibship F-F
    near(coancestry(MATERNAL_HALF_BROTHERS), 0.5); # maternal half-sibship M-M
    near(coancestry(PATERNAL_GRANDMOTHER_GRANDDAUGHTER), 0.25);
    near(coancestry(MATERNAL_GRANDMOTHER_GRANDDAUGHTER), 0.125);
    near(coancestry(MATERNAL_GRANDFATHER_GRANDDAUGHTER), 0.25);
    near(coancestry(PATERNAL_AUNT_NIECE), 0.125); # paternal avuncular F-F
    near(coancestry(PATERNAL_UNCLE_NIECE), 0.25); # paternal avuncular M-F
    near(coancestry(MATERNAL_AUNT_NIECE), 0.1875); # maternal avuncular F-F, 3/16
}

func testCoefficientsSumToOne() {
    for (def r in [
        FATHER_DAUGHTER,
        FATHER_SON,
        MOTHER_SON,
        MOTHER_DAUGHTER,
        FULL_SISTERS,
        FULL_BROTHERS,
        BROTHER_SISTER,
        PATERNAL_HALF_SISTERS,
        MATERNAL_HALF_SISTERS,
        MATERNAL_HALF_BROTHERS,
        PATERNAL_GRANDMOTHER_GRANDDAUGHTER,
        MATERNAL_GRANDMOTHER_GRANDDAUGHTER,
        MATERNAL_GRANDFATHER_GRANDDAUGHTER,
        PATERNAL_AUNT_NIECE,
        PATERNAL_UNCLE_NIECE,
        MATERNAL_AUNT_NIECE
    ]) {
        # xRelationship validates exactly what a constant must satisfy.
        def same as XRelationship init xRelationship(
            $r.name,
            $r.sexA,
            $r.sexB,
            $r.k0,
            $r.k1,
            $r.k2);
        near($same.k0 + $same.k1 + $same.k2, 1.0);
    }
}

func testUnrelated() {
    def u as XRelationship init unrelated(MALE, FEMALE);
    near($u.k0, 1.0);
    testing.assertEqual($u.sexB, "female");
    near(coancestry($u), 0.0);
}

func maleK2() {
    xRelationship("bad", MALE, FEMALE, 0.0, 0.5, 0.5);
}

func badSum() {
    xRelationship("bad", FEMALE, FEMALE, 0.5, 0.4, 0.2);
}

func badRange() {
    xRelationship("bad", FEMALE, FEMALE, -0.5, 1.5, 0.0);
}

func badSex() {
    xRelationship("bad", "m", FEMALE, 1.0, 0.0, 0.0);
}

func badSexCoancestry() {
    coancestry(XRelationship{name: "x", sexA: "?", sexB: MALE, k0: 1.0, k1: 0.0, k2: 0.0});
}

func testRelationshipValidation() {
    testing.assertThrows("maleK2", "xstr");
    testing.assertThrows("badSum", "xstr");
    testing.assertThrows("badRange", "xstr");
    testing.assertThrows("badSex", "xstr");
    testing.assertThrows("badSexCoancestry", "xstr");
}

# --- profiles ---------------------------------------------------------------

func testProfiles() {
    def m as XProfile init male("AF", "10");
    testing.assertEqual($m.sex, MALE);
    testing.assertEqual($m.loci["L"].a, "10");
    testing.assertEqual($m.loci["L"].b, "10");
    def f as XProfile init female("C", "11", "10");
    testing.assertEqual($f.sex, FEMALE);
    testing.assertEqual($f.loci["L"].a, "10");
}

func femaleTwice() {
    femaleProfile("F", [forensics.genotype("L", "10", "11"), forensics.genotype("L", "12", "12")]);
}

func maleEmptyAllele() {
    maleProfile("M", {"L": ""});
}

func maleWithTwoAlleles() {
    def p as XProfile init XProfile{
        id: "bad",
        sex: MALE,
        loci: {"L": forensics.genotype("L", "10", "11")}
    };
    matchProbability(model(), $p);
}

func testProfileValidation() {
    testing.assertThrows("femaleTwice", "xstr");
    testing.assertThrows("maleEmptyAllele", "xstr");
    testing.assertThrows("maleWithTwoAlleles", "xstr");
}

# --- linkage groups ---------------------------------------------------------

func testHaplotypeKeyAndFrequency() {
    testing.assertEqual(haplotypeKey(["17", "24.1", "12"]), "17-24.1-12");
    def g as LinkageGroup init linked().groups[0];
    near(haplotypeFrequency($g, ["10", "11"]), 0.1);
    # unlisted: scored at the floor, .001
    near(haplotypeFrequency($g, ["13", "13"]), 0.001);
}

func keyWithSeparator() {
    haplotypeKey(["1-2"]);
}

func emptyKey() {
    haplotypeKey([]);
}

func wrongArity() {
    haplotypeFrequency(linked().groups[0], ["10"]);
}

func unlistedNoFloor() {
    def g as LinkageGroup init linkageGroup("G", ["L", "M"], 0.0, {"10-10": 1.0});
    haplotypeFrequency($g, ["11", "11"]);
}

func groupNoLoci() {
    linkageGroup("G", [], 0.0, {});
}

func groupRepeatedLocus() {
    linkageGroup("G", ["L", "L"], 0.0, {});
}

func groupBadFloor() {
    linkageGroup("G", ["L", "M"], 1.0, {});
}

func groupBadKey() {
    linkageGroup("G", ["L", "M"], 0.0, {"10": 0.5});
}

func groupEmptyAllele() {
    linkageGroup("G", ["L", "M"], 0.0, {"10-": 0.5});
}

func groupBadFrequency() {
    linkageGroup("G", ["L", "M"], 0.0, {"10-10": 1.5});
}

func groupBadSum() {
    linkageGroup("G", ["L", "M"], 0.0, {"10-10": 0.6, "11-11": 0.6});
}

func groupTwiceByName() {
    withLinkageGroup(linked(), linkageGroup("LM", ["X1", "X2"], 0.0, {}));
}

func groupOverlap() {
    withLinkageGroup(linked(), linkageGroup("other", ["M", "N"], 0.0, {}));
}

func testLinkageGroupValidation() {
    testing.assertThrows("keyWithSeparator", "xstr");
    testing.assertThrows("emptyKey", "xstr");
    testing.assertThrows("wrongArity", "xstr");
    testing.assertThrows("unlistedNoFloor", "xstr");
    testing.assertThrows("groupNoLoci", "xstr");
    testing.assertThrows("groupRepeatedLocus", "xstr");
    testing.assertThrows("groupBadFloor", "xstr");
    testing.assertThrows("groupBadKey", "xstr");
    testing.assertThrows("groupEmptyAllele", "xstr");
    testing.assertThrows("groupBadFrequency", "xstr");
    testing.assertThrows("groupBadSum", "xstr");
    testing.assertThrows("groupTwiceByName", "xstr");
    testing.assertThrows("groupOverlap", "xstr");
}

func testTwoDisjointGroups() {
    def xf as XFrequencies init withLinkageGroup(linked(), linkageGroup("NO", ["N", "O"], 0.0, {}));
    testing.assertEqual(len($xf.groups), 2);
    testing.assertEqual($xf.groups[1].name, "NO");
}

func testArgusGroupsAsInTheHandbook() {
    # QIAGEN Investigator Argus X-12 QS Handbook (07/2022), Table 1: four groups
    # of three, twelve distinct loci.
    def groups as map of string to list of string init argusX12Groups();
    testing.assertEqual(len($groups), 4);
    def seen as map of string to bool init {};
    for (def name in $groups) {
        testing.assertEqual(len($groups[$name]), 3);
        for (def locus in $groups[$name]) {
            $seen[$locus] = true;
        }
    }
    testing.assertEqual(len($seen), 12);
    testing.assertEqual($groups["LG1"], ["DXS8378", "DXS10135", "DXS10148"]);
    testing.assertEqual($groups["LG2"], ["DXS7132", "DXS10074", "DXS10079"]);
    testing.assertEqual($groups["LG3"], ["HPRTB", "DXS10101", "DXS10103"]);
    testing.assertEqual($groups["LG4"], ["DXS7423", "DXS10134", "DXS10146"]);
}

# --- match probability ------------------------------------------------------

func testMatchProbabilityIndependent() {
    # male L 10, M 12: .1 * .3 = .03
    near(matchProbability(model(), male2("m", "10", "12")), 0.03);
    # female L 10/11, M 12/12: (2 * .1 * .2) * .3^2 = .04 * .09 = .0036
    near(matchProbability(model(), female2("f", "10", "11", "12", "12")), 0.0036);
    near(log10MatchProbability(model(), female2("f", "10", "11", "12", "12")), math.log10(0.0036));
}

func testMatchProbabilityLinked() {
    # male haplotype 10-11: listed .1
    near(matchProbability(linked(), male2("m", "10", "11")), 0.1);
    # female 10/11 at both: phases {10-10, 11-11} and {10-11, 11-10}:
    # 2 * .3 * .3 + 2 * .1 * .1 = .18 + .02 = .2
    near(matchProbability(linked(), female2("f", "10", "11", "10", "11")), 0.2);
    # female 10/10 at L, 10/11 at M: one phase {10-10, 10-11}: 2 * .3 * .1 = .06
    near(matchProbability(linked(), female2("f", "10", "10", "10", "11")), 0.06);
    near(log10MatchProbability(linked(), male2("m", "10", "11")), -1.0);
}

func testLinkageEquilibriumGroupEqualsIndependent() {
    # Double heterozygote, haplotypes at product frequencies:
    # 2 (.1 * .1)(.2 * .2) + 2 (.1 * .2)(.2 * .1) = .0016 = (2 * .1 * .2)^2
    near(matchProbability(equilibrium(), female2("f", "10", "11", "10", "11")), 0.0016);
    near(matchProbability(equilibrium(), male2("m", "12", "13")), 0.12);
}

func partialGroup() {
    matchProbability(linked(), male("m", "10"));
}

func uncoveredLocus() {
    matchProbability(model(), maleProfile("m", {"Q": "10"}));
}

func emptyProfile() {
    matchProbability(model(), maleProfile("m", {}));
}

func testMatchProbabilityErrors() {
    testing.assertThrows("partialGroup", "xstr");
    testing.assertThrows("uncoveredLocus", "xstr");
    testing.assertThrows("emptyProfile", "xstr");
}

# --- pairwise kinship, independent loci -------------------------------------

func testFatherDaughter() {
    # k1 = 1, M-F: px Q(G, x) / (px P(G)) = Q(G, x) / P(G).
    # 10 vs 10/11: .2 / (2 * .1 * .2) = 5
    near(kinshipLR(male("AF", "10"), female("C", "10", "11"), FATHER_DAUGHTER, model()), 5.0);
    # 10 vs 10/10: .1 / .01 = 10
    near(kinshipLR(male("AF", "10"), female("C", "10", "10"), FATHER_DAUGHTER, model()), 10.0);
    # 10 vs 11/12: the daughter lacks his allele - exclusion
    near(kinshipLR(male("AF", "10"), female("C", "11", "12"), FATHER_DAUGHTER, model()), 0.0);
}

func testMotherSon() {
    # F-M, k1 = 1: Q(G, y) / P(G) = (.2) / (2 * .1 * .2)... per son allele 10:
    # Q(10/11, 10) = p11 = .2; P(10/11) = .04; ratio 5 = T(M, 10) / p10 = .5 / .1
    near(kinshipLR(female("M", "10", "11"), male("S", "10"), MOTHER_SON, model()), 5.0);
}

func testFatherSonIsUninformative() {
    near(kinshipLR(male("F", "10"), male("S", "11"), FATHER_SON, model()), 1.0);
}

func testFullBrothers() {
    # M-M: (k0 p^2 + k1 p) / p^2 = .5 + .5 / .1 = 5.5
    near(kinshipLR(male("A", "10"), male("B", "10"), FULL_BROTHERS, model()), 5.5);
    # different alleles: k0 only, .5
    near(kinshipLR(male("A", "10"), male("B", "11"), FULL_BROTHERS, model()), 0.5);
}

func testBrotherSister() {
    # k0 + k1 Q(G, x) / P(G) = .5 + .5 * .2 / .04 = 3
    near(kinshipLR(male("B", "10"), female("S", "10", "11"), BROTHER_SISTER, model()), 3.0);
}

func testFullSisters() {
    # 10/11 vs 10/11: k1 term (.5 * .2 + .5 * .1) / (2 * .1 * .2) = .15 / .04 = 3.75;
    # k2 term 1 / .04 = 25; LR = .5 * 3.75 + .5 * 25 = 14.375
    def a as XProfile init female("A", "10", "11");
    def b as XProfile init female("B", "10", "11");
    near(kinshipLR($a, $b, FULL_SISTERS, model()), 14.375);
    # full sisters always share the paternal allele: 10/11 vs 12/13 is excluded
    near(kinshipLR($a, female("C", "12", "13"), FULL_SISTERS, model()), 0.0);
    # paternal half sisters (k1 = 1) against full sisters: 3.75 / 14.375
    near(kinshipLRAgainst($a, $b, PATERNAL_HALF_SISTERS, FULL_SISTERS, model()), 3.75 / 14.375);
}

func testProductAcrossLoci() {
    # L: 10 vs 10/11 -> 5; M: 12 vs 12/13 -> .4 / (2 * .3 * .4) = 1 / .6
    def lr as float init kinshipLR(
        male2("AF", "10", "12"),
        female2("C", "10", "11", "12", "13"),
        FATHER_DAUGHTER,
        model());
    near($lr, 5.0 / 0.6);
    testing.assertEqual(
        pairUnits(model(), male2("AF", "10", "12"), female("C", "10", "11")),
        ["L"]);
}

func testFemalePairsEqualAutosomalKinship() {
    # Two females: the X formula is the autosomal one with the X coefficients.
    def gs as list of forensics.Genotype init [
        forensics.genotype("L", "10", "11"),
        forensics.genotype("L", "10", "10"),
        forensics.genotype("L", "11", "12"),
        forensics.genotype("L", "12", "13")
    ];
    for (def r in [FULL_SISTERS, MATERNAL_HALF_SISTERS, MATERNAL_AUNT_NIECE, MOTHER_DAUGHTER]) {
        def k as kinship.Relationship init kinship.relationship($r.name, $r.k0, $r.k1, $r.k2);
        for (def ga in $gs) {
            for (def gb in $gs) {
                def x as float init kinshipLR(
                    femaleProfile("A", [$ga]),
                    femaleProfile("B", [$gb]),
                    $r,
                    model());
                near($x, kinship.kinshipLRAt(fixture(), $ga, $gb, $k, kinship.UNRELATED));
            }
        }
    }
}

func testJointProbabilitySumsToMarginal() {
    # Summing the joint over every genotype of B gives P(A): male B for a
    # female A, female B for a male A.
    def u as Unit init Unit{name: "L", loci: ["L"], group: -1};
    def alleles as list of string init ["10", "11", "12", "13"];
    def a as XProfile init female("A", "10", "11");
    def total as float init 0.0;
    for (def x in $alleles) {
        $total = $total + unitJoint(model(), $u, $a, male("B", $x), MOTHER_SON);
    }
    near($total, 0.04);
    def m as XProfile init male("A", "12");
    $total = 0.0;
    for (def i in 0..4) {
        for (def j in $i..4) {
            def b as XProfile init female("B", $alleles[$i], $alleles[$j]);
            $total = $total + unitJoint(model(), $u, $m, $b, BROTHER_SISTER);
        }
    }
    near($total, 0.3);
}

func testSymmetry() {
    def a as XProfile init female("A", "10", "11");
    def b as XProfile init female("B", "11", "12");
    near(
        kinshipLR($a, $b, MATERNAL_AUNT_NIECE, model()),
        kinshipLR($b, $a, MATERNAL_AUNT_NIECE, model()));
}

func testLog10KinshipLR() {
    def lr as float init log10KinshipLRAgainst(
        male("AF", "10"),
        female("C", "10", "10"),
        FATHER_DAUGHTER,
        unrelated(MALE, FEMALE),
        model());
    near($lr, 1.0);
}

func log10Exclusion() {
    log10KinshipLRAgainst(
        male("AF", "10"),
        female("C", "11", "12"),
        FATHER_DAUGHTER,
        unrelated(MALE, FEMALE),
        model());
}

func sexMismatch() {
    kinshipLR(female("A", "10", "11"), female("B", "10", "11"), FATHER_DAUGHTER, model());
}

func altZero() {
    kinshipLRAgainst(
        female("A", "10", "11"),
        female("B", "12", "13"),
        MATERNAL_HALF_SISTERS,
        FULL_SISTERS,
        model());
}

func noSharedLocus() {
    kinshipLR(maleProfile("A", {"M": "10"}), female("B", "10", "11"), FATHER_DAUGHTER, model());
}

func pairUncovered() {
    kinshipLR(
        maleProfile("A", {"Q": "10"}),
        femaleProfile("B", [forensics.genotype("Q", "10", "11")]),
        FATHER_DAUGHTER,
        model());
}

func testKinshipErrors() {
    testing.assertThrows("log10Exclusion", "xstr");
    testing.assertThrows("sexMismatch", "xstr");
    testing.assertThrows("altZero", "xstr");
    testing.assertThrows("noSharedLocus", "xstr");
    testing.assertThrows("pairUncovered", "xstr");
}

# --- pairwise kinship, linkage groups ---------------------------------------

func testLinkedBrothers() {
    # Haplotype 10-10 at .1 * .1 = .01 (equilibrium table), shared as one
    # block: .5 + .5 / .01 = 50.5; the independent model multiplies two loci:
    # (.5 + .5 / .1)^2 = 30.25.
    def a as XProfile init male2("A", "10", "10");
    def b as XProfile init male2("B", "10", "10");
    near(kinshipLR($a, $b, FULL_BROTHERS, equilibrium()), 50.5);
    near(kinshipLR($a, $b, FULL_BROTHERS, model()), 30.25);
    testing.assertEqual(pairUnits(equilibrium(), $a, $b), ["LM"]);
}

func testSingleMeiosisUnderEquilibriumMatchesIndependent() {
    # Mother-daughter: one meiosis, the mother's phase unknown; with haplotypes
    # in equilibrium the block and the per-locus product agree.
    def m as XProfile init female2("M", "10", "11", "10", "12");
    def d as XProfile init female2("D", "10", "13", "12", "12");
    near(
        kinshipLR($m, $d, MOTHER_DAUGHTER, equilibrium()),
        kinshipLR($m, $d, MOTHER_DAUGHTER, model()));
}

func testLinkedMotherDaughterByHand() {
    # Mother 10/11, 10/11: phases {10-10, 11-11} (.18), {10-11, 11-10} (.02).
    # Daughter 10/12, 10/12: phases {10-10, 12-12} (2 * .3 * .2 = .12) and
    # {10-12, 12-10} (both unlisted, floor .001: 2e-6).
    # k1 = 1: only the mother's 10-10 meets the daughter's first phase:
    #   .18 * (1/2) * Q({10-10, 12-12}, 10-10) = .18 * .5 * .2 = .018
    # unrelated: P(M) P(D) = .2 * (.12 + .000002) = .0240004
    def m as XProfile init female2("M", "10", "11", "10", "11");
    def d as XProfile init female2("D", "10", "12", "10", "12");
    near(kinshipLR($m, $d, MOTHER_DAUGHTER, linked()), 0.018 / 0.0240004);
}

func testOneLocusGroupEqualsIndependent() {
    def g as LinkageGroup init linkageGroup(
        "L only",
        ["L"],
        0.0,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    def xf as XFrequencies init withLinkageGroup(model(), $g);
    def a as XProfile init female("A", "10", "11");
    def b as XProfile init female("B", "10", "12");
    near(kinshipLR($a, $b, FULL_SISTERS, $xf), kinshipLR($a, $b, FULL_SISTERS, model()));
}

func pairPartialGroup() {
    kinshipLR(male2("A", "10", "10"), female("B", "10", "11"), FATHER_DAUGHTER, linked());
}

func testPartialGroupThrowsForEitherPerson() {
    testing.assertThrows("pairPartialGroup", "xstr");
}

# --- paternity of a daughter ------------------------------------------------

func testTrioIndex() {
    def m as XProfile init female("M", "10", "11");
    # daughter 10/12, father 12: X = T(10) * 1 = .5, Y = .5 * .3 -> 1 / .3
    near(
        daughterTrioIndex(
            $m,
            female("C", "10", "12"),
            male("AF", "12"),
            model(),
            kinship.NO_MUTATION),
        1.0 / 0.3);
    # daughter 10/11, father 11: X = .5 * 1 + .5 * 0 = .5; Y = .5 * .2 + .5 * .1 = .15
    def c as XProfile init female("C", "10", "11");
    near(daughterTrioIndex($m, $c, male("AF", "11"), model(), kinship.NO_MUTATION), 0.5 / 0.15);
    # father 12: no allele to give - 0 without mutation
    near(daughterTrioIndex($m, $c, male("AF", "12"), model(), kinship.NO_MUTATION), 0.0);
    # rate .002, decay .1: F(11 | 12) = .002 * .9 / 2 = .0009,
    # F(10 | 12) = .002 * .9 * .1 / 2 = .00009;
    # X = .5 * .0009 + .5 * .00009 = .000495, PI = .000495 / .15 = .0033
    def mm as kinship.MutationModel init kinship.mutationModel(0.002, 0.1);
    near(daughterTrioIndex($m, $c, male("AF", "12"), model(), $mm), 0.0033);
}

func testDuoIndex() {
    # daughter 10/11, father 10: X = F(10) p11 = .2, Y = .04 -> 5
    def c as XProfile init female("C", "10", "11");
    near(daughterDuoIndex($c, male("AF", "10"), model(), kinship.NO_MUTATION), 5.0);
    # father 12, rate .002, decay .1: X = .00009 * .2 + .0009 * .1 = .000108;
    # Y = .04 -> .0027
    def mm as kinship.MutationModel init kinship.mutationModel(0.002, 0.1);
    near(daughterDuoIndex($c, male("AF", "12"), model(), $mm), 0.0027);
}

func testDuoEqualsFatherDaughterLR() {
    def alleles as list of string init ["10", "11", "12", "13"];
    for (def f in $alleles) {
        for (def i in 0..4) {
            for (def j in $i..4) {
                def c as XProfile init female("C", $alleles[$i], $alleles[$j]);
                near(
                    daughterDuoIndex($c, male("AF", $f), model(), kinship.NO_MUTATION),
                    kinshipLR(male("AF", $f), $c, FATHER_DAUGHTER, model()));
            }
        }
    }
}

func testLinkedTrioAndDuo() {
    # Mother 10/11, 10/11: phase weights .18 / .2 = .9 and .02 / .2 = .1.
    # Daughter 10/12, 10/12: phases {10-10, 12-12} and {10-12, 12-10};
    # father 12-12.
    # T_M(10-10) = .9 * .5 = .45; T_M(12-12) = T_M(10-12) = T_M(12-10) = 0.
    # X = .45 * 1 = .45; Y = .45 * p(12-12) = .45 * .2 = .09; PI = 5.
    def m as XProfile init female2("M", "10", "11", "10", "11");
    def d as XProfile init female2("D", "10", "12", "10", "12");
    near(daughterTrioIndex($m, $d, male2("AF", "12", "12"), linked(), kinship.NO_MUTATION), 5.0);
    # Father 12-13 under rate .002, decay .1: F(12-12 | 12-13) =
    # (1 - .002) * (.002 * .9 / 2) = .998 * .0009 = .0008982; PI = 5 * .0008982
    def mm as kinship.MutationModel init kinship.mutationModel(0.002, 0.1);
    near(daughterTrioIndex($m, $d, male2("AF", "12", "13"), linked(), $mm), 5.0 * 0.0008982);
    # Duo, father 12-12: X = F(12-12) p(10-10) = .3 (second phase: 0);
    # Y = .12 + .000002
    near(
        daughterDuoIndex($d, male2("AF", "12", "12"), linked(), kinship.NO_MUTATION),
        0.3 / 0.120002);
}

func maternalExclusion() {
    daughterTrioIndex(
        female("M", "12", "13"),
        female("C", "10", "11"),
        male("AF", "10"),
        model(),
        kinship.NO_MUTATION);
}

func trioWrongSex() {
    daughterTrioIndex(
        female("M", "10", "11"),
        female("C", "10", "11"),
        female("AF", "10", "10"),
        model(),
        kinship.NO_MUTATION);
}

func duoWrongSex() {
    daughterDuoIndex(male("C", "10"), male("AF", "10"), model(), kinship.NO_MUTATION);
}

func testPaternityErrors() {
    testing.assertThrows("maternalExclusion", "xstr");
    testing.assertThrows("trioWrongSex", "xstr");
    testing.assertThrows("duoWrongSex", "xstr");
}
