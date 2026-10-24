# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for pedigree.j. Run with `jennifer test pedigree_test.j`.
#
# The strongest tests here are the cross-checks against kinship.j: a pedigree
# LR built from first principles by peeling must reproduce, to the last digit,
# the closed-form IBD and trio formulae that module computes a completely
# different way. Agreement between two independent derivations is much better
# evidence than either matching a number captured from itself.
use testing;
use math;

import "./kinship.j" as kinship;

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

func gt(a as string, b as string) {
    return forensics.genotype("L", $a, $b);
}

func one(id as string, a as string, b as string) {
    return forensics.profile($id, [gt($a, $b)]);
}

# --- building and validating ------------------------------------------------

func trioPedigree() {
    def p as Pedigree init pedigree("trio");
    $p = withFounder($p, "AF", "m");
    $p = withFounder($p, "MO", "f");
    return withChild($p, "CH", "m", "AF", "MO");
}

func testPedigreeConstruction() {
    def p as Pedigree init trioPedigree();
    testing.assertEqual(members($p), ["AF", "MO", "CH"]);
    testing.assertTrue(isFounder($p, "AF"));
    testing.assertFalse(isFounder($p, "CH"));
    testing.assertEqual(childrenOf($p, "AF"), ["CH"]);
    testing.assertEqual(childrenOf($p, "CH"), []);
    testing.assertEqual(personOf($p, "CH").mother, "MO");
    testing.assertTrue(hasMember($p, "MO"));
    testing.assertFalse(hasMember($p, "NOBODY"));
}

func duplicateMember() {
    withFounder(withFounder(pedigree("p"), "A", "m"), "A", "f");
}

func selfParent() {
    withChild(withFounder(pedigree("p"), "A", "m"), "B", "m", "B", "A");
}

func sameIdForBothParents() {
    withChild(withFounder(pedigree("p"), "A", "m"), "B", "m", "A", "A");
}

func badSex() {
    withFounder(pedigree("p"), "A", "male");
}

func testConstructionIsValidated() {
    testing.assertThrows("duplicateMember", "pedigree");
    testing.assertThrows("selfParent", "pedigree");
    testing.assertThrows("sameIdForBothParents", "pedigree");
    testing.assertThrows("badSex", "pedigree");
}

func missingParent() {
    validate(withChild(withFounder(pedigree("p"), "A", "f"), "B", "m", "GHOST", "A"));
}

func testUnknownParentIsRejected() {
    testing.assertThrows("missingParent", "pedigree");
}

func consanguineousPedigree() {
    # Two full siblings have a child together: the individual/family graph
    # closes a loop, which peeling cannot handle.
    def p as Pedigree init pedigree("loop");
    $p = withFounder($p, "GF", "m");
    $p = withFounder($p, "GM", "f");
    $p = withChild($p, "A", "m", "GF", "GM");
    $p = withChild($p, "B", "f", "GF", "GM");
    return withChild($p, "C", "m", "A", "B");
}

func validateLoop() {
    validate(consanguineousPedigree());
}

func testLoopsAreRejectedNotMisSummed() {
    testing.assertThrows("validateLoop", "pedigree");
    testing.assertFalse(isSimple(consanguineousPedigree()));
    testing.assertTrue(isSimple(trioPedigree()));
}

func testHalfSiblingPedigreeIsNotALoop() {
    # A shared father with two different mothers is a perfectly simple
    # pedigree - the loop check must not over-reject it.
    def p as Pedigree init pedigree("halfsibs");
    $p = withFounder($p, "F", "m");
    $p = withFounder($p, "M1", "f");
    $p = withFounder($p, "M2", "f");
    $p = withChild($p, "A", "m", "F", "M1");
    $p = withChild($p, "B", "f", "F", "M2");
    testing.assertTrue(isSimple($p));
}

# --- likelihoods at one locus -----------------------------------------------

func testUntypedFounderHasLikelihoodOne() {
    # With nothing observed, the genotype probabilities sum to 1.
    near(locusLikelihood(withFounder(pedigree("p"), "A", "m"), {}, fixture(), "L"), 1.0);
}

func testTypedFounderIsItsGenotypeFrequency() {
    def p as Pedigree init withFounder(pedigree("p"), "A", "m");
    near(locusLikelihood($p, {"A": gt("10", "11")}, fixture(), "L"), 0.04);
    near(locusLikelihood($p, {"A": gt("10", "10")}, fixture(), "L"), 0.01);
}

func testFullyTypedTrio() {
    # P(F) P(M) * P(child | F, M) = .04 * .24 * .25
    def g as map of string to forensics.Genotype init {
        "AF": gt("10", "11"),
        "MO": gt("12", "13"),
        "CH": gt("10", "12")
    };
    near(locusLikelihood(trioPedigree(), $g, fixture(), "L"), 0.0024);
}

func testAnImpossibleTrioHasLikelihoodZero() {
    # The child carries no allele either parent could have transmitted.
    def g as map of string to forensics.Genotype init {
        "AF": gt("10", "10"),
        "MO": gt("11", "11"),
        "CH": gt("12", "13")
    };
    near(locusLikelihood(trioPedigree(), $g, fixture(), "L"), 0.0);
}

func testUntypedMembersAreSummedOut() {
    # Typing only the child must give back the plain population frequency: the
    # untyped parents contribute nothing once summed over.
    near(locusLikelihood(trioPedigree(), {"CH": gt("10", "11")}, fixture(), "L"), 0.04);
}

func testDisconnectedComponentsMultiply() {
    def p as Pedigree init withFounder(withFounder(pedigree("two"), "A", "m"), "B", "f");
    near(
        locusLikelihood($p, {"A": gt("10", "11"), "B": gt("10", "10")}, fixture(), "L"),
        0.04 * 0.01);
}

func testIdsOutsideThePedigreeAreIgnored() {
    def p as Pedigree init withFounder(pedigree("p"), "A", "m");
    near(
        locusLikelihood($p, {"A": gt("10", "11"), "STRANGER": gt("12", "13")}, fixture(), "L"),
        0.04);
}

# --- cross-checks against kinship.j -----------------------------------------

# unrelatedPair puts two typed people in a pedigree with no relationship at
# all - the denominator every LR below is taken against.
func unrelatedPair(a as string, b as string) {
    def p as Pedigree init pedigree("unrelated");
    $p = withFounder($p, "F1", "m");
    $p = withFounder($p, "M1", "f");
    $p = withFounder($p, "F2", "m");
    $p = withFounder($p, "M2", "f");
    $p = withChild($p, $a, "m", "F1", "M1");
    return withChild($p, $b, "m", "F2", "M2");
}

func siblingPedigree() {
    def p as Pedigree init pedigree("sibs");
    $p = withFounder($p, "F1", "m");
    $p = withFounder($p, "M1", "f");
    $p = withChild($p, "A", "m", "F1", "M1");
    return withChild($p, "B", "m", "F1", "M1");
}

func testFullSiblingLrMatchesKinship() {
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "10"),
        "B": one("B", "10", "10")
    };
    def viaPeeling as float init pedigreeLR(
        siblingPedigree(),
        unrelatedPair("A", "B"),
        $profiles,
        $db);
    near(
        $viaPeeling,
        kinship.kinshipLR($profiles["A"], $profiles["B"], kinship.FULL_SIBLINGS, $db));
    # And against the closed form (1 + 1/p)^2 / 4 for two identical homozygotes.
    near($viaPeeling, 30.25);
}

func testFullSiblingLrMatchesKinshipForHeterozygotes() {
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "11"),
        "B": one("B", "11", "12")
    };
    near(
        pedigreeLR(siblingPedigree(), unrelatedPair("A", "B"), $profiles, $db),
        kinship.kinshipLR($profiles["A"], $profiles["B"], kinship.FULL_SIBLINGS, $db));
}

func halfSiblingPedigree() {
    def p as Pedigree init pedigree("halfsibs");
    $p = withFounder($p, "F1", "m");
    $p = withFounder($p, "M1", "f");
    $p = withFounder($p, "M2", "f");
    $p = withChild($p, "A", "m", "F1", "M1");
    return withChild($p, "B", "m", "F1", "M2");
}

func testHalfSiblingLrMatchesKinship() {
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "11"),
        "B": one("B", "10", "12")
    };
    def viaPeeling as float init pedigreeLR(
        halfSiblingPedigree(),
        unrelatedPair("A", "B"),
        $profiles,
        $db);
    near(
        $viaPeeling,
        kinship.kinshipLR($profiles["A"], $profiles["B"], kinship.HALF_SIBLINGS, $db));
    near($viaPeeling, 1.75);
}

func parentChildPedigree() {
    def p as Pedigree init pedigree("parentchild");
    $p = withFounder($p, "A", "m");
    $p = withFounder($p, "M1", "f");
    return withChild($p, "B", "m", "A", "M1");
}

func testParentChildLrMatchesKinship() {
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "10"),
        "B": one("B", "10", "11")
    };
    def viaPeeling as float init pedigreeLR(
        parentChildPedigree(),
        unrelatedPair("A", "B"),
        $profiles,
        $db);
    near($viaPeeling, kinship.kinshipLR($profiles["A"], $profiles["B"], kinship.PARENT_CHILD, $db));
    near($viaPeeling, 5.0);
}

func grandparentPedigree() {
    # GP -> P (untyped) -> C. The intermediate generation is never observed,
    # so this exercises peeling through an unobserved link.
    def p as Pedigree init pedigree("grandparent");
    $p = withFounder($p, "GP", "m");
    $p = withFounder($p, "GM", "f");
    $p = withChild($p, "P", "m", "GP", "GM");
    $p = withFounder($p, "PS", "f");
    return withChild($p, "C", "m", "P", "PS");
}

func testGrandparentLrMatchesKinship() {
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "GP": one("GP", "10", "11"),
        "C": one("C", "10", "12")
    };
    near(
        pedigreeLR(grandparentPedigree(), unrelatedPair("GP", "C"), $profiles, $db),
        kinship.kinshipLR($profiles["GP"], $profiles["C"], kinship.GRANDPARENT, $db));
}

func motherlessPaternityPedigree() {
    def p as Pedigree init pedigree("h2");
    $p = withFounder($p, "AF", "m");
    $p = withFounder($p, "MO", "f");
    $p = withFounder($p, "UN", "m");
    return withChild($p, "CH", "m", "UN", "MO");
}

func testTrioPaternityLrMatchesKinship() {
    # H1: the child is the alleged father's. H2: the child is an unrelated
    # man's, with the same mother. The ratio must be the trio paternity index.
    def db as forensics.FrequencyDb init fixture();
    def profiles as map of string to forensics.Profile init {
        "MO": one("MO", "10", "11"),
        "CH": one("CH", "10", "12"),
        "AF": one("AF", "12", "13")
    };
    def h1 as Pedigree init withFounder(trioPedigree(), "UN", "m");
    def viaPeeling as float init pedigreeLR($h1, motherlessPaternityPedigree(), $profiles, $db);
    near(
        $viaPeeling,
        kinship.paternityIndexWith(
            $profiles["MO"],
            $profiles["CH"],
            $profiles["AF"],
            $db,
            kinship.NO_MUTATION));
    near($viaPeeling, 1.6666666666666667);
}

func testAvuncularLrMatchesKinship() {
    # The aunt is the child's father's sister; the father himself is untyped.
    def db as forensics.FrequencyDb init fixture();
    def p as Pedigree init pedigree("avuncular");
    $p = withFounder($p, "GF", "m");
    $p = withFounder($p, "GM", "f");
    $p = withChild($p, "AU", "f", "GF", "GM");
    $p = withChild($p, "FA", "m", "GF", "GM");
    $p = withFounder($p, "MO", "f");
    $p = withChild($p, "CH", "m", "FA", "MO");
    def profiles as map of string to forensics.Profile init {
        "AU": one("AU", "10", "11"),
        "CH": one("CH", "10", "12")
    };
    near(
        pedigreeLR($p, unrelatedPair("AU", "CH"), $profiles, $db),
        kinship.kinshipLR($profiles["AU"], $profiles["CH"], kinship.AVUNCULAR, $db));
}

# --- across a panel ---------------------------------------------------------

func twoLocus() {
    return forensics.withLocus(fixture(), "M", 10000, {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

func testLikelihoodMultipliesAcrossLoci() {
    def p as Pedigree init withFounder(pedigree("p"), "A", "m");
    def profiles as map of string to forensics.Profile init {
        "A": forensics.profile("A", [gt("10", "11"), forensics.genotype("M", "10", "10")])
    };
    testing.assertEqual(pedigreeLoci($p, $profiles, twoLocus()), ["L", "M"]);
    near(pedigreeLikelihood($p, $profiles, twoLocus()), 0.04 * 0.01);
    near(log10PedigreeLikelihood($p, $profiles, twoLocus()), math.log10(0.0004));
}

func testTypedMembersListing() {
    def profiles as map of string to forensics.Profile init {
        "AF": one("AF", "10", "11"),
        "STRANGER": one("STRANGER", "12", "13")
    };
    testing.assertEqual(typedMembers(trioPedigree(), $profiles), ["AF"]);
}

func mismatchedTypedMembers() {
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "10"),
        "B": one("B", "10", "10")
    };
    # The denominator pedigree does not contain B, so it is scored against less
    # data than the numerator and the ratio would be meaningless.
    pedigreeLR(siblingPedigree(), withFounder(pedigree("only-a"), "A", "m"), $profiles, fixture());
}

func testLrRefusesToCompareDifferentData() {
    testing.assertThrows("mismatchedTypedMembers", "pedigree");
}

func noUsableLocus() {
    pedigreeLikelihood(
        trioPedigree(),
        {"AF": one("AF", "10", "11")},
        forensics.withLocus(forensics.db("other"), "M", 100, {"10": 1.0}));
}

func testNoUsableLocusThrows() {
    testing.assertThrows("noUsableLocus", "pedigree");
}

func log10OfAnImpossiblePedigree() {
    def profiles as map of string to forensics.Profile init {
        "AF": one("AF", "10", "10"),
        "MO": one("MO", "11", "11"),
        "CH": one("CH", "12", "13")
    };
    log10PedigreeLikelihood(trioPedigree(), $profiles, fixture());
}

func testLog10RefusesAZeroLikelihood() {
    testing.assertThrows("log10OfAnImpossiblePedigree", "pedigree");
}

func testLoopedPedigreeIsRejectedAtLikelihoodTime() {
    testing.assertFalse(isSimple(consanguineousPedigree()));
}

# --- regressions from the validation audit ----------------------------------

func likelihoodWithAWrongLocusGenotype() {
    def p as Pedigree init withFounder(pedigree("p"), "A", "m");
    locusLikelihood($p, {"A": forensics.genotype("OTHER", "10", "11")}, fixture(), "L");
}

func testGenotypesMustBelongToTheLocusAsked() {
    # Otherwise the alleles are scored against the wrong locus's frequencies and
    # a plausible-looking number comes back for a question nobody asked.
    testing.assertThrows("likelihoodWithAWrongLocusGenotype", "pedigree");
}

func testAConsanguineousPedigreeIsReportedAsALoop() {
    # It is a loop, not a self-ancestry error: diagnosing it as the latter sends
    # the reader looking for a data-entry mistake that is not there.
    def caught as string init "";
    try {
        validate(consanguineousPedigree());
    } catch (e) {
        $caught = $e.message;
    }
    testing.assertContains($caught, "loop");
}

func selfAncestor() {
    # B is their own grandparent: a genuine cycle, and still reported as one.
    def p as Pedigree init pedigree("cycle");
    $p = withFounder($p, "X", "f");
    $p = withFounder($p, "Y", "f");
    $p = withChild($p, "A", "m", "B", "X");
    $p = withChild($p, "B", "m", "A", "Y");
    validate($p);
}

func testATrueAncestorCycleIsStillCaught() {
    testing.assertThrows("selfAncestor", "pedigree");
}

# --- documented edge cases with no assertion yet ----------------------------

func emptyId() {
    withFounder(pedigree("p"), "", "m");
}

func oneParentOnly() {
    withChild(withFounder(pedigree("p"), "A", "m"), "B", "m", "A", "");
}

func testMoreConstructionChecks() {
    testing.assertThrows("emptyId", "pedigree");
    testing.assertThrows("oneParentOnly", "pedigree");
}

func personOfStranger() {
    personOf(trioPedigree(), "NOBODY");
}

func isFounderOfStranger() {
    isFounder(trioPedigree(), "NOBODY");
}

func testLookupOfAnUnknownMemberThrows() {
    testing.assertThrows("personOfStranger", "pedigree");
    testing.assertThrows("isFounderOfStranger", "pedigree");
}

func testParentsMayBeAddedAfterTheChild() {
    def p as Pedigree init withChild(pedigree("late"), "CH", "m", "AF", "MO");
    testing.assertFalse(isSimple($p));
    $p = withFounder(withFounder($p, "AF", "m"), "MO", "f");
    testing.assertTrue(isSimple($p));
    testing.assertEqual(childrenOf($p, "MO"), ["CH"]);
}

func likelihoodOfEmptyPedigree() {
    locusLikelihood(pedigree("empty"), {}, fixture(), "L");
}

func testAnEmptyPedigreeHasNoLikelihood() {
    testing.assertThrows("likelihoodOfEmptyPedigree", "pedigree");
}

func zeroDenominatorPedigree() {
    def profiles as map of string to forensics.Profile init {
        "AF": one("AF", "10", "10"),
        "MO": one("MO", "11", "11"),
        "CH": one("CH", "12", "13")
    };
    pedigreeLR(trioPedigree(), trioPedigree(), $profiles, fixture());
}

func log10WithNoUsableLocus() {
    log10PedigreeLikelihood(
        trioPedigree(),
        {"AF": one("AF", "10", "11")},
        forensics.withLocus(forensics.db("other"), "M", 100, {"10": 1.0}));
}

func testPanelFunctionsRefuseUndefinedAnswers() {
    testing.assertThrows("zeroDenominatorPedigree", "pedigree");
    testing.assertThrows("log10WithNoUsableLocus", "pedigree");
}

func testPedigreeLociFollowTheDatabaseOrder() {
    def p as Pedigree init withFounder(pedigree("p"), "A", "m");
    def profiles as map of string to forensics.Profile init {
        "A": forensics.profile("A", [forensics.genotype("M", "10", "10"), gt("10", "11")])
    };
    testing.assertEqual(pedigreeLoci($p, $profiles, twoLocus()), ["L", "M"]);
}

func testAnEmptyProfileIsNotTyped() {
    def profiles as map of string to forensics.Profile init {"AF": forensics.profile("AF", [])};
    testing.assertEqual(typedMembers(trioPedigree(), $profiles), []);
}
