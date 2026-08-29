# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for lineage.j. Run with `jennifer test lineage_test.j`.
#
# The Clopper-Pearson expectations are computed independently by bisecting the
# exact binomial CDF, not captured from this implementation.
use testing;
use math;
use strings;
use convert;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    testing.assertTrue(math.abs($actual - $expected) < TOL);
}

# --- counting estimators ----------------------------------------------------

func testHaplotypeFrequencyIsThePointEstimate() {
    near(haplotypeFrequency(2, 5000), 0.0004);
    near(haplotypeFrequency(0, 100), 0.0);
    near(haplotypeFrequency(7, 7), 1.0);
}

func testClopperPearsonUpperBoundAtZeroSightings() {
    # k = 0 has the closed form 1 - alpha^(1/N); at 95% that is the rule of three.
    near(haplotypeUpperBound(0, 100, 0.95), 0.029513049607039987);
    near(haplotypeUpperBound(0, 100, 0.95), 1.0 - math.pow(0.05, 1.0 / 100.0));
}

func testClopperPearsonUpperBound() {
    near(haplotypeUpperBound(2, 5000, 0.95), 0.0012586179249663432);
    near(haplotypeUpperBound(1, 1000, 0.95), 0.004734993575499724);
    near(haplotypeUpperBound(5, 200, 0.95), 0.051843339120290544);
}

func testUpperBoundIsAlwaysAboveThePointEstimate() {
    # The whole point of the conservative estimate: never quote k/N directly.
    testing.assertTrue(haplotypeUpperBound(2, 5000, 0.95) > haplotypeFrequency(2, 5000));
    testing.assertTrue(haplotypeUpperBound(0, 100, 0.95) > haplotypeFrequency(0, 100));
}

func testUpperBoundIsOneWhenEveryProfileMatches() {
    near(haplotypeUpperBound(50, 50, 0.95), 1.0);
}

func testHigherConfidenceIsMoreConservative() {
    testing.assertTrue(haplotypeUpperBound(2, 5000, 0.99) > haplotypeUpperBound(2, 5000, 0.95));
}

func testLineageMatchProbabilityAndLr() {
    near(lineageMatchProbability(2, 5000), 0.0012586179249663432);
    near(lineageLikelihoodRatio(2, 5000), 794.5222931945301);
}

func zeroDatabase() {
    haplotypeFrequency(1, 0);
}

func countAboveDatabaseSize() {
    haplotypeFrequency(11, 10);
}

func negativeCount() {
    haplotypeFrequency(-1, 10);
}

func confidenceBelowHalf() {
    haplotypeUpperBound(1, 10, 0.4);
}

func testCountAndConfidenceAreValidated() {
    testing.assertThrows("zeroDatabase", "lineage");
    testing.assertThrows("countAboveDatabaseSize", "lineage");
    testing.assertThrows("negativeCount", "lineage");
    testing.assertThrows("confidenceBelowHalf", "lineage");
}

func testLrRefusesAWorthlessBound() {
    # k == N gives a bound of 1.0, so the LR is 1 - no weight to report.
    near(lineageLikelihoodRatio(10, 10), 1.0);
}

# --- rCRS difference nomenclature -------------------------------------------

func testParseSubstitution() {
    def d as Difference init parseDifference("16189C");
    testing.assertEqual($d.position, 16189);
    testing.assertEqual($d.insertion, 0);
    testing.assertEqual($d.base, "C");
    testing.assertEqual($d.kind, "sub");
}

func testParseDropsTheRedundantReferenceBase() {
    # "A263G" and "263G" name the same difference: the position fixes the
    # reference base, so it carries no information.
    testing.assertEqual(parseDifference("A263G"), parseDifference("263G"));
}

func testParseInsertion() {
    def d as Difference init parseDifference("315.1C");
    testing.assertEqual($d.position, 315);
    testing.assertEqual($d.insertion, 1);
    testing.assertEqual($d.base, "C");
    testing.assertEqual($d.kind, "ins");
}

func testParseDeletionForms() {
    def d as Difference init parseDifference("249d");
    testing.assertEqual($d.position, 249);
    testing.assertEqual($d.kind, "del");
    testing.assertEqual($d.base, "");
    testing.assertEqual(parseDifference("249D"), $d);
    testing.assertEqual(parseDifference("249DEL"), $d);
    testing.assertEqual(parseDifference("249-"), $d);
}

func testParseAcceptsIupacAmbiguity() {
    # A heteroplasmic position is reported with an ambiguity code.
    testing.assertEqual(parseDifference("16093Y").base, "Y");
}

func testFormatRoundTrips() {
    for (def s in ["16189C", "315.1C", "249d", "73G"]) {
        testing.assertEqual(formatDifference(parseDifference($s)), $s);
    }
}

func badNoPosition() {
    parseDifference("CG");
}

func badNoBase() {
    parseDifference("263");
}

func badDoubleBase() {
    parseDifference("263GG");
}

func badInsertionIndex() {
    parseDifference("315.C");
}

func badDeletedInsertion() {
    parseDifference("315.1d");
}

func testMalformedDesignationsThrow() {
    testing.assertThrows("badNoPosition", "lineage");
    testing.assertThrows("badNoBase", "lineage");
    testing.assertThrows("badDoubleBase", "lineage");
    testing.assertThrows("badInsertionIndex", "lineage");
    testing.assertThrows("badDeletedInsertion", "lineage");
}

# --- mtDNA haplotypes -------------------------------------------------------

func testMtHaplotypeSortsByPositionThenInsertion() {
    def h as MtHaplotype init mtHaplotype(
        "Q1",
        ["16024-16365", "73-340"],
        "16189C 315.1C 263G 315.2C");
    testing.assertEqual(mtDesignations($h), ["263G", "315.1C", "315.2C", "16189C"]);
}

func testMtHaplotypeSplitsOnAnyWhitespace() {
    def h as MtHaplotype init mtHaplotype("Q1", [], "263G\t315.1C\n16189C");
    testing.assertEqual(len($h.differences), 3);
}

func testMtMatchesAndDistance() {
    def a as MtHaplotype init mtHaplotype("A", ["16024-16365"], "16189C 16223T 16278T");
    def b as MtHaplotype init mtHaplotype("B", ["16024-16365"], "16223T 16189C 16278T");
    testing.assertTrue(mtMatches($a, $b));
    testing.assertEqual(mtDistance($a, $b), 0);
    # One extra difference: distance 1.
    def c as MtHaplotype init mtHaplotype("C", ["16024-16365"], "16189C 16223T 16278T 16311C");
    testing.assertEqual(mtDistance($a, $c), 1);
    testing.assertFalse(mtMatches($a, $c));
    # A different base at a shared position counts twice - each side carries a
    # designation the other lacks.
    def d as MtHaplotype init mtHaplotype("D", ["16024-16365"], "16189C 16223C 16278T");
    testing.assertEqual(mtDistance($a, $d), 2);
}

func testMtRangesAgree() {
    def a as MtHaplotype init mtHaplotype("A", ["16024-16365", "73-340"], "263G");
    def b as MtHaplotype init mtHaplotype("B", ["73-340", "16024-16365"], "263G");
    def c as MtHaplotype init mtHaplotype("C", ["16024-16365"], "263G");
    testing.assertTrue(mtRangesAgree($a, $b));
    testing.assertFalse(mtRangesAgree($a, $c));
    # Same differences, different coverage: matching is not like-for-like.
    testing.assertTrue(mtMatches($a, $c));
    testing.assertFalse(mtRangesAgree($a, $c));
}

# --- alignment --------------------------------------------------------------

func designations(ds as list of Difference) {
    def out as list of string init [];
    for (def d in $ds) {
        $out[] = formatDifference($d);
    }
    return $out;
}

func testAlignFindsSubstitutions() {
    # positions 100..103 = A C G T; the sample reads A C T T.
    testing.assertEqual(designations(align("ACGT", "ACTT", 100)), ["102T"]);
    testing.assertEqual(designations(align("ACGT", "ACGT", 100)), []);
}

func testAlignFindsADeletion() {
    # positions 100..104 = A C G T A; G is missing from the sample.
    testing.assertEqual(designations(align("ACGTA", "ACTA", 100)), ["102d"]);
}

func testAlignPlacesAnInsertionAtTheThreePrimeEndOfAHomopolymer() {
    # The C-stretch case: an extra C in a run of three starting at 311 must be
    # numbered from the run's 3' end (313.1C), not its 5' end (310.1C).
    testing.assertEqual(designations(align("ACCCG", "ACCCCG", 310)), ["313.1C"]);
}

func testAlignHandlesLeadingAndTrailingChange() {
    testing.assertEqual(designations(align("ACGT", "TCGT", 1)), ["1T"]);
    testing.assertEqual(designations(align("ACGT", "ACGA", 1)), ["4A"]);
}

func testAlignIgnoresWhitespaceInSequences() {
    # A FASTA block pasted with line breaks aligns as one sequence.
    testing.assertEqual(designations(align("AC GT\n", "ACTT", 100)), ["102T"]);
}

func testMtHaplotypeFromRecordsItsRange() {
    def h as MtHaplotype init mtHaplotypeFrom("S1", "ACGT", "ACTT", 100);
    testing.assertEqual($h.ranges, ["100-103"]);
    testing.assertEqual(mtDesignations($h), ["102T"]);
}

func alignEmptySample() {
    align("ACGT", "", 100);
}

func alignBadCharacter() {
    align("ACGT", "ACZT", 100);
}

func alignZeroStart() {
    align("ACGT", "ACGT", 0);
}

func testAlignValidatesItsInput() {
    testing.assertThrows("alignEmptySample", "lineage");
    testing.assertThrows("alignBadCharacter", "lineage");
    testing.assertThrows("alignZeroStart", "lineage");
}

# --- Y-STR haplotypes -------------------------------------------------------

func yA() {
    return yHaplotype("A", {"DYS19": "14", "DYS385": "11,14", "DYS389I": "13", "DYS390": "24"});
}

func testYMatchesIsOrderInsensitiveAtMultiCopyLoci() {
    def b as YHaplotype init yHaplotype(
        "B",
        {"DYS19": "14", "DYS385": "14,11", "DYS389I": "13", "DYS390": "24"});
    testing.assertTrue(yMatches(yA(), $b));
    near(yStepDistance(yA(), $b), 0.0);
}

func testYStepDistanceCountsRepeatUnits() {
    def b as YHaplotype init yHaplotype(
        "B",
        {"DYS19": "15", "DYS385": "11,14", "DYS389I": "13", "DYS390": "24"});
    near(yStepDistance(yA(), $b), 1.0);
    testing.assertFalse(yMatches(yA(), $b));
    testing.assertTrue(yMatchesWithin(yA(), $b, 1.0));
    testing.assertFalse(yMatchesWithin(yA(), $b, 0.0));
}

func testYStepDistanceSumsAcrossLoci() {
    def b as YHaplotype init yHaplotype(
        "B",
        {"DYS19": "15", "DYS385": "11,16", "DYS389I": "13", "DYS390": "24"});
    near(yStepDistance(yA(), $b), 3.0);
    testing.assertFalse(yMatchesWithin(yA(), $b, 1.0));
}

func testYComparisonUsesOnlySharedLoci() {
    def b as YHaplotype init yHaplotype("B", {"DYS19": "14", "DYS390": "24", "DYS391": "10"});
    testing.assertEqual(ySharedLoci(yA(), $b), ["DYS19", "DYS390"]);
    testing.assertTrue(yMatches(yA(), $b));
    near(yStepDistance(yA(), $b), 0.0);
}

func testYLociListing() {
    testing.assertEqual(yLoci(yA()), ["DYS19", "DYS385", "DYS389I", "DYS390"]);
}

func yNoSharedLocus() {
    yMatches(yA(), yHaplotype("B", {"DYS437": "15"}));
}

func yCopyCountMismatch() {
    yStepDistance(yA(), yHaplotype("B", {"DYS385": "11"}));
}

func yNegativeSteps() {
    yMatchesWithin(yA(), yA(), -1.0);
}

func testYComparisonValidatesItsInput() {
    testing.assertThrows("yNoSharedLocus", "lineage");
    testing.assertThrows("yCopyCountMismatch", "lineage");
    testing.assertThrows("yNegativeSteps", "lineage");
}

# --- regressions from the validation audit ----------------------------------

func positionPastTheSortWidth() {
    parseDifference("1234567G");
}

func positionZero() {
    parseDifference("0G");
}

func testPositionIsBounded() {
    # Differences sort through a fixed-width key, so a position wider than that
    # key would order wrongly. The mitochondrial genome is 16569 bases, so the
    # bound is far above anything real.
    testing.assertThrows("positionPastTheSortWidth", "lineage");
    testing.assertThrows("positionZero", "lineage");
    # Six digits still parse, and still sort correctly.
    def h as MtHaplotype init mtHaplotype("H", [], "999999G 16189C 263G");
    testing.assertEqual(mtDesignations($h), ["263G", "16189C", "999999G"]);
}

func stepDistanceOnANonNumericAllele() {
    yStepDistance(yA(), yHaplotype("b", {"DYS19": "OL"}));
}

func buildWithANonNumericAllele() {
    yHaplotype("bad", {"DYS19": "OL"});
}

func testYAlleleValidationRaisesLineageErrors() {
    # Not a bare toFloat failure from `convert` - the caller catches "lineage".
    testing.assertThrows("stepDistanceOnANonNumericAllele", "lineage");
    testing.assertThrows("buildWithANonNumericAllele", "lineage");
}

func alignBeyondTheLengthCap() {
    def long as string init strings.repeat("ACGT", (MAX_ALIGN_LENGTH // 4) + 2);
    align($long, $long, 1);
}

func testAlignRefusesAnOversizedRegion() {
    # O(nm): a whole-genome region would appear to hang. Refuse it instead.
    testing.assertThrows("alignBeyondTheLengthCap", "lineage");
}

# --- banded alignment: differential against a full-width pass ---------------

# The aligner computes only a band around the diagonal and widens it if the
# traceback touches the edge. That is a correctness-critical shortcut in a
# module whose whole job is producing a haplotype designation, so it is checked
# the only way worth checking: against a full-width pass, on randomised input,
# demanding byte-identical output.

func randomSeq(n as int) {
    def bases as list of string init ["A", "C", "G", "T"];
    def out as string init "";
    for (def i in 0..$n) {
        $out = $out + $bases[math.randInt(0, 3)];
    }
    return $out;
}

func mutate(seq as string, subs as int, dels as int, ins as int) {
    def ch as list of string init strings.chars($seq);
    def bases as list of string init ["A", "C", "G", "T"];
    for (def k in 0..$subs) {
        $ch[math.randInt(0, len($ch) - 1)] = $bases[math.randInt(0, 3)];
    }
    for (def k in 0..$dels) {
        $ch = dropAt($ch, math.randInt(0, len($ch) - 1));
    }
    for (def k in 0..$ins) {
        $ch = insertAt($ch, math.randInt(0, len($ch)), $bases[math.randInt(0, 3)]);
    }
    return strings.join($ch, "");
}

func dropAt(ch as list of string, at as int) {
    def out as list of string init [];
    for (def i in 0..len($ch)) {
        if ($i != $at) {
            $out[] = $ch[$i];
        }
    }
    return $out;
}

func insertAt(ch as list of string, at as int, base as string) {
    def out as list of string init [];
    for (def i in 0..len($ch)) {
        if ($i == $at) {
            $out[] = $base;
        }
        $out[] = $ch[$i];
    }
    if ($at == len($ch)) {
        $out[] = $base;
    }
    return $out;
}

# assertSameAsFullWidth aligns one pair both ways and compares. It returns
# whether the initial band was too narrow, so a test can show the widening path
# is actually being exercised rather than silently never taken.
func assertSameAsFullWidth(ref as string, smp as string) {
    def rb as bytes init convert.bytesFromString($ref, "utf-8");
    def sb as bytes init convert.bytesFromString($smp, "utf-8");
    def r as list of string init strings.chars($ref);
    def s as list of string init strings.chars($smp);
    def full as Band init alignBand($rb, $sb, $r, $s, math.max(len($r), len($s)));
    def banded as list of string init alignPair($ref, $smp);
    testing.assertEqual($banded[0], $full.top);
    testing.assertEqual($banded[1], $full.bottom);
    return not alignBand($rb, $sb, $r, $s, INITIAL_BAND).exact;
}

func testBandedAlignmentMatchesFullWidth() {
    math.randSeed(20260828);
    for (def trial in 0..60) {
        def ref as string init randomSeq(math.randInt(20, 90));
        assertSameAsFullWidth(
            $ref,
            mutate($ref, math.randInt(0, 6), math.randInt(0, 4), math.randInt(0, 4)));
    }
}

func testBandedAlignmentAgreesOnDivergentSequences() {
    # Unrelated sequences: the optimal path wanders far from the diagonal, so
    # the band has to widen. Every one of these must still match full width.
    math.randSeed(99);
    def widened as int init 0;
    for (def trial in 0..12) {
        def a as string init randomSeq(math.randInt(30, 60));
        def b as string init randomSeq(math.randInt(30, 60));
        if (assertSameAsFullWidth($a, $b)) {
            $widened = $widened + 1;
        }
    }
    testing.assertTrue($widened > 0);
}

func testALongIndelForcesTheBandWider() {
    # A 40-base deletion pushes the path further off the diagonal than the
    # initial band half-width, so this exercises the widening loop directly.
    def ref as string init strings.repeat("ACGT", 40);
    def smp as string init strings.substring($ref, 0, 20) + strings.substring($ref, 120, 160);
    testing.assertTrue(assertSameAsFullWidth($ref, $smp));
}
