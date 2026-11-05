# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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

# tinyTable has N = 20, so the floor is 5/40 = 0.125 and allele 10, listed at
# 0.05, floors up to it. The floored frequencies then sum to 1.075.
func tinyTable() {
    return forensics.withLocus(
        forensics.db("tiny"),
        "L",
        20,
        {"10": 0.05, "11": 0.2, "12": 0.3, "13": 0.45});
}

func testAFlooredAlleleStillAgreesWithKinship() {
    # Siblings 10/11 and 10/12 against unrelated, with p10 at its floor:
    #   LR = k0 + k1 (p12 / 2) / (2 p10 p12) = 1/4 + 1/(8 p10) = 0.25 + 1 = 1.25
    # Only 10, 11 and 12 are observed, so the floored sum is 0.625 and the
    # lumped class carries the rest; the floor is in force and both modules
    # must still agree.
    def db as forensics.FrequencyDb init tinyTable();
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "11"),
        "B": one("B", "10", "12")
    };
    def viaPeeling as float init pedigreeLR(
        siblingPedigree(),
        unrelatedPair("A", "B"),
        $profiles,
        $db);
    near($viaPeeling, 1.25);
    near(
        $viaPeeling,
        kinship.kinshipLR($profiles["A"], $profiles["B"], kinship.FULL_SIBLINGS, $db));
}

func flooredSumPastOne() {
    # A typed bystander C shows the 13 too, so every allele is observed and
    # the floored frequencies sum to 1.075. Rescaling them to 1 would put p10
    # at 0.125/1.075, under the NRC II floor, and give 1.325 where kinship.j
    # gives 1.25; keeping them makes an untyped founder's prior sum to
    # 1.075^2 instead of 1, so the ratio would depend on how many untyped
    # people each hypothesis names. The locus has no honest answer.
    def profiles as map of string to forensics.Profile init {
        "A": one("A", "10", "11"),
        "B": one("B", "10", "12"),
        "C": one("C", "13", "13")
    };
    pedigreeLR(
        withFounder(siblingPedigree(), "C", "m"),
        withFounder(unrelatedPair("A", "B"), "C", "m"),
        $profiles,
        tinyTable());
}

func testAFlooredSumPastOneIsRefusedNotRescaled() {
    testing.assertThrows("flooredSumPastOne", "pedigree");
    def caught as string init "";
    try {
        flooredSumPastOne();
    } catch (e) {
        $caught = $e.message;
    }
    testing.assertContains($caught, "sum to 1.075");
}

func likelihoodOverAnUncoveredLocus() {
    def profiles as map of string to forensics.Profile init {
        "A": forensics.profile("A", [gt("10", "11"), forensics.genotype("X", "10", "10")])
    };
    pedigreeLikelihood(withFounder(pedigree("p"), "A", "m"), $profiles, fixture());
}

func lrOverAnUncoveredLocus() {
    def profiles as map of string to forensics.Profile init {
        "A": forensics.profile("A", [gt("10", "10"), forensics.genotype("X", "10", "10")]),
        "B": forensics.profile("B", [gt("10", "10"), forensics.genotype("X", "12", "13")])
    };
    pedigreeLR(siblingPedigree(), unrelatedPair("A", "B"), $profiles, fixture());
}

func testAnUncoveredTypedLocusThrowsInsteadOfBeingDropped() {
    # At X the two could not be siblings of a 10/10 - that evidence must not
    # vanish because the table happens to lack the locus.
    testing.assertThrows("likelihoodOverAnUncoveredLocus", "pedigree");
    testing.assertThrows("lrOverAnUncoveredLocus", "pedigree");
    # pedigreeLoci remains the query that says what the table can use.
    def profiles as map of string to forensics.Profile init {
        "A": forensics.profile("A", [gt("10", "11"), forensics.genotype("X", "10", "10")])
    };
    testing.assertEqual(
        pedigreeLoci(withFounder(pedigree("p"), "A", "m"), $profiles, fixture()),
        ["L"]);
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

# --- review fixes -------------------------------------------------------------

func testAFlooredSumOfExactlyOneIsAccepted() {
    # .05 + .55 + .3 + .1 is 1 but adds up to 1.0000000000000002 in floats,
    # which was refused as past 1. Two unrelated typed founders:
    #   2 * .05 * .55 * 2 * .3 * .1 = .055 * .06 = .0033
    def db as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("t"),
        "L",
        10000,
        {"10": 0.05, "11": 0.55, "12": 0.3, "13": 0.1});
    def p as Pedigree init withFounder(withFounder(pedigree("p"), "A", "m"), "B", "f");
    near(locusLikelihood($p, {"A": gt("10", "11"), "B": gt("12", "13")}, $db, "L"), 0.0033);
}

# chain is `n` generations: C1 is the child of the founders F0 and M0, and
# each C(i+1) is the child of C(i) and a founder spouse S(i+1).
func chain(n as int) {
    def m as map of string to Person init {};
    $m["F0"] = Person{id: "F0", sex: "m", father: "", mother: ""};
    $m["M0"] = Person{id: "M0", sex: "f", father: "", mother: ""};
    def prevF as string init "F0";
    def prevM as string init "M0";
    for (def i in 1..($n + 1)) {
        $m["C{$i}"] = Person{id: "C{$i}", sex: "m", father: $prevF, mother: $prevM};
        $m["S{$i}"] = Person{id: "S{$i}", sex: "f", father: "", mother: ""};
        $prevF = "C{$i}";
        $prevM = "S{$i}";
    }
    return Pedigree{name: "chain", members: $m};
}

func testADeepPedigreeNeedsNoCallStack() {
    # 1200 generations (recursive peeling ran out of call stack near 1000).
    # F0 10/11 and C1200 12/13 share an allele IBD with probability about
    # 2^-1200, so the likelihood is the unrelated product to float precision:
    #   (2 * .1 * .2) * (2 * .3 * .4) = .04 * .24 = .0096
    def g as map of string to forensics.Genotype init {
        "F0": gt("10", "11"),
        "C1200": gt("12", "13")
    };
    near(locusLikelihood(chain(1200), $g, fixture(), "L"), 0.0096);
}

# farPanel builds `n` loci listing only allele 11, so an observed 10 scores
# the floor 5/20000 = .00025.
func farPanel(n as int) {
    def out as forensics.FrequencyDb init forensics.db("t");
    for (def i in 0..$n) {
        $out = forensics.withLocus($out, "L{$i}", 10000, {"11": 0.5});
    }
    return $out;
}

func farProfile(id as string, n as int, a as string, b as string) {
    def gs as list of forensics.Genotype init [];
    for (def i in 0..$n) {
        $gs[] = forensics.genotype("L{$i}", $a, $b);
    }
    return forensics.profile($id, $gs);
}

# fatherOf: AF is CH's father, the mother an untyped founder. Against
# unrelatedPair("AF", "CH") with AF and CH both 10/10 the LR per locus is
#   P(CH | AF) / P(CH) = p10 / p10^2 = 1 / .00025 = 4000,
# so 100 loci give 4000^100 = 10^360.2.
func fatherOf() {
    return withChild(
        withFounder(withFounder(pedigree("duo"), "AF", "m"), "MO", "f"),
        "CH",
        "m",
        "AF",
        "MO");
}

func pedigreeLrPastTheFloatRange() {
    def profiles as map of string to forensics.Profile init {
        "AF": farProfile("AF", 100, "10", "10"),
        "CH": farProfile("CH", 100, "10", "10")
    };
    pedigreeLR(fatherOf(), unrelatedPair("AF", "CH"), $profiles, farPanel(100));
}

func testPedigreeLrPastTheFloatRangeHasALog10Form() {
    # Before the fix: a runtime float-overflow error.
    testing.assertThrows("pedigreeLrPastTheFloatRange", "pedigree");
    def profiles as map of string to forensics.Profile init {
        "AF": farProfile("AF", 100, "10", "10"),
        "CH": farProfile("CH", 100, "10", "10")
    };
    near(
        log10PedigreeLR(fatherOf(), unrelatedPair("AF", "CH"), $profiles, farPanel(100)),
        100.0 * math.log10(4000.0));
}

func log10PedigreeLrOfAnExclusion() {
    # AF 10/10 cannot father CH 11/11: the numerator is 0 at L.
    def profiles as map of string to forensics.Profile init {
        "AF": one("AF", "10", "10"),
        "CH": one("CH", "11", "11")
    };
    log10PedigreeLR(fatherOf(), unrelatedPair("AF", "CH"), $profiles, fixture());
}

func testLog10PedigreeLrRefusesAnExclusion() {
    testing.assertThrows("log10PedigreeLrOfAnExclusion", "pedigree");
}

# --- round 2 review fixes ---------------------------------------------------

# literalTrio is the trio as a Pedigree literal, with CH's parents as given.
func literalTrio(father as string, mother as string) {
    return Pedigree{
        name: "lit",
        members: {
            "AF": Person{id: "AF", sex: "m", father: "", mother: ""},
            "MO": Person{id: "MO", sex: "f", father: "", mother: ""},
            "CH": Person{id: "CH", sex: "m", father: $father, mother: $mother}
        }
    };
}

func motherOnlyLikelihood() {
    # MO 10/10, CH 11/11: a maternal exclusion, L = 0. Treated as a founder,
    # CH gave .01 * .04 = .0004.
    locusLikelihood(
        literalTrio("", "MO"),
        {"MO": gt("10", "10"), "CH": gt("11", "11")},
        fixture(),
        "L");
}

func fatherOnlyLiteral() {
    validate(literalTrio("AF", ""));
}

func sameParentLiteral() {
    validate(literalTrio("MO", "MO"));
}

func messageOf(name as string) {
    def caught as string init "";
    try {
        if ($name == "motherOnlyLikelihood") {
            motherOnlyLikelihood();
        } elseif ($name == "fatherOnlyLiteral") {
            fatherOnlyLiteral();
        } elseif ($name == "sameParentLiteral") {
            sameParentLiteral();
        } else {
            tooManyAlleles();
        }
    } catch (e) {
        $caught = $e.message;
    }
    return $caught;
}

func testALiteralWithOneParentIsRejected() {
    testing.assertThrows("motherOnlyLikelihood", "pedigree");
    testing.assertThrows("fatherOnlyLiteral", "pedigree");
    testing.assertContains(messageOf("motherOnlyLikelihood"), "must name both parents or neither");
    testing.assertContains(messageOf("fatherOnlyLiteral"), "must name both parents or neither");
    testing.assertFalse(isFounder(literalTrio("", "MO"), "CH"));
    testing.assertFalse(isSimple(literalTrio("", "MO")));
}

func testALiteralWithOneIdForBothParentsSaysSo() {
    testing.assertThrows("sameParentLiteral", "pedigree");
    testing.assertContains(messageOf("sameParentLiteral"), "names 'MO' as both parents");
}

# untypedCouple is founders F and M with `k` children C0..C(k-1).
func untypedCouple(k as int) {
    def p as Pedigree init withFounder(withFounder(pedigree("d"), "F", "m"), "M", "f");
    for (def i in 0..$k) {
        $p = withChild($p, "C{$i}", "", "F", "M");
    }
    return $p;
}

func tooManyAlleles() {
    # Eleven typed children with pairwise-distinct alleles show 22 > 20.
    def g as map of string to forensics.Genotype init {};
    for (def i in 0..11) {
        $g["C{$i}"] = forensics.genotype("L", "{100 + 2 * $i}", "{101 + 2 * $i}");
    }
    locusLikelihood(
        untypedCouple(11),
        $g,
        forensics.withLocus(forensics.db("t"), "L", 100000, {"1": 0.5}),
        "L");
}

func testTooManyAllelesAtALocusThrows() {
    testing.assertThrows("tooManyAlleles", "pedigree");
    testing.assertContains(messageOf("tooManyAlleles"), "MAX_PEDIGREE_ALLELES (20)");
    # At the limit: ten unrelated founders show 20 unlisted alleles, each at
    # the floor f = 5 / 200000, so L = (2 f^2)^10.
    def p as Pedigree init pedigree("f");
    def g as map of string to forensics.Genotype init {};
    for (def i in 0..10) {
        $p = withFounder($p, "P{$i}", "");
        $g["P{$i}"] = forensics.genotype("L", "{100 + 2 * $i}", "{101 + 2 * $i}");
    }
    def f as float init 5.0 / 200000.0;
    def l as float init locusLikelihood(
        $p,
        $g,
        forensics.withLocus(forensics.db("t"), "L", 100000, {"1": 0.5}),
        "L");
    near($l / math.pow(2.0 * $f * $f, 10.0), 1.0);
}

func testTypedParentsSkipNothingTheySum() {
    # Typed parents make childWeights fill one cell, a typed spouse one
    # column. AF 10/11, MO 12/13 typed; children C1 untyped with typed child
    # G1 (C1 x founder W 10/10), C2 10/12. G1 10/12:
    #   P(C2 10/12 | AF, MO) = 1/4
    #   C1 in {10/12, 10/13, 11/12, 11/13}, 1/4 each; G1 gets 10 from W (1)
    #   and 12 from C1: 1/2 for 10/12 and 11/12, 0 otherwise -> 1/4.
    #   L = (2 .1 .2)(2 .3 .4)(.1^2) * 1/4 * 1/4 = .04 * .24 * .01 / 16
    def p as Pedigree init withFounder(
        withFounder(withFounder(pedigree("t"), "AF", "m"), "MO", "f"),
        "W",
        "f");
    $p = withChild(withChild($p, "C1", "m", "AF", "MO"), "C2", "f", "AF", "MO");
    $p = withChild($p, "G1", "m", "C1", "W");
    def g as map of string to forensics.Genotype init {
        "AF": gt("10", "11"),
        "MO": gt("12", "13"),
        "W": gt("10", "10"),
        "C2": gt("10", "12"),
        "G1": gt("10", "12")
    };
    near(locusLikelihood($p, $g, fixture(), "L"), 0.04 * 0.24 * 0.01 / 16.0);
}

func pedigreeOfMotherOnly() {
    pedigreeOf(
        "p",
        [
            Person{id: "MO", sex: "f", father: "", mother: ""},
            Person{id: "CH", sex: "m", father: "", mother: "MO"}
        ]);
}

func pedigreeOfDuplicate() {
    pedigreeOf(
        "p",
        [
            Person{id: "A", sex: "f", father: "", mother: ""},
            Person{id: "A", sex: "m", father: "", mother: ""}
        ]);
}

func pedigreeOfBadSex() {
    pedigreeOf("p", [Person{id: "A", sex: "x", father: "", mother: ""}]);
}

func testPedigreeOfMatchesTheFold() {
    def built as Pedigree init pedigreeOf(
        "trio",
        [
            Person{id: "AF", sex: "m", father: "", mother: ""},
            Person{id: "MO", sex: "f", father: "", mother: ""},
            Person{id: "CH", sex: "m", father: "AF", mother: "MO"}
        ]);
    testing.assertEqual($built, trioPedigree());
    testing.assertThrows("pedigreeOfMotherOnly", "pedigree");
    testing.assertThrows("pedigreeOfDuplicate", "pedigree");
    testing.assertThrows("pedigreeOfBadSex", "pedigree");
}
