# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# White-box test overlay for lineage.j. Run with `jennifer test lineage_test.j`.
#
# The Clopper-Pearson expectations are computed independently by bisecting the
# exact binomial CDF, not captured from this implementation.
use testing;
use math;
use strings;
use convert;
use time;

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

func testMtDistanceOnLiteralsKeepsSetSemantics() {
    def a as MtHaplotype init mtHaplotype("A", [], "16189C 16223T 263G");
    # out of order and repeated: the merge does not apply; as sets, equal
    def lit as MtHaplotype init MtHaplotype{
        id: "L",
        ranges: [],
        differences: [
            parseDifference("16223T"),
            parseDifference("263G"),
            parseDifference("16189C"),
            parseDifference("263G")
        ]
    };
    testing.assertEqual(mtDistance($a, $lit), 0);
    testing.assertEqual(mtDistance($lit, $a), 0);
    # a deletion literal carrying a base still renders "249d"
    def del as MtHaplotype init MtHaplotype{
        id: "D",
        ranges: [],
        differences: [Difference{position: 249, insertion: 0, base: "X", kind: "del"}]
    };
    testing.assertEqual(mtDistance($del, mtHaplotype("E", [], "249d")), 0);
    # 249d against 249A: one designation each side
    testing.assertEqual(mtDistance($del, mtHaplotype("F", [], "249A")), 2);
}

func testMtDistanceMergeEqualsTheSetDifference() {
    # random designation sets over a few positions, with substitutions,
    # insertions and deletions, against the designation-set count
    def kinds as list of string init ["C", "T", "d", ".1C", ".2A", "A"];
    math.randSeed(7);
    for (def round in 0..200) {
        def ta as list of string init [];
        def tb as list of string init [];
        for (def pos in [73, 150, 263, 315, 16189, 16223]) {
            if (math.randInt(0, 2) > 0) {
                $ta[] = "{$pos}" + $kinds[math.randInt(0, 5)];
            }
            if (math.randInt(0, 2) > 0) {
                $tb[] = "{$pos}" + $kinds[math.randInt(0, 5)];
            }
        }
        def a as MtHaplotype init mtHaplotype("A", [], strings.join($ta, " "));
        def b as MtHaplotype init mtHaplotype("B", [], strings.join($tb, " "));
        testing.assertEqual(mergeDistance($a, $b), setDistance($a, $b));
    }
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

func alignPastMaxPosition() {
    # the last base would be 999999 + 3 = 1000002; before the check this
    # emitted 1000000G, which parseDifference rejects
    align("ACGT", "AGGT", 999999);
}

func fromAtMaxInt() {
    # before the check: a runtime integer overflow computing the range end
    mtHaplotypeFrom("S", "ACGT", "ACGT", 9223372036854775807);
}

func testAlignValidatesItsInput() {
    testing.assertThrows("alignEmptySample", "lineage");
    testing.assertThrows("alignBadCharacter", "lineage");
    testing.assertThrows("alignZeroStart", "lineage");
    testing.assertThrows("alignPastMaxPosition", "lineage");
    testing.assertThrows("fromAtMaxInt", "lineage");
}

func testAlignAtTheLastPositionRoundTrips() {
    # refStart 999996: the four bases sit at 999996..999999, the C→G at 999997
    def ds as list of Difference init align("ACGT", "AGGT", 999996);
    testing.assertEqual(designations($ds), ["999997G"]);
    testing.assertEqual(parseDifference("999997G"), $ds[0]);
    def h as MtHaplotype init mtHaplotypeFrom("S", "ACGT", "AGGT", 999996);
    testing.assertEqual($h.ranges, ["999996-999999"]);
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

func buildWithAHugeAllele() {
    # 400 digits: before the cap, a runtime toFloat failure
    yHaplotype("bad", {"DYS19": strings.repeat("9", 400)});
}

func stepDistanceOnASevenDigitAllele() {
    yStepDistance(yA(), yHaplotype("b", {"DYS19": "1234567"}));
}

func matchOnALongFraction() {
    yMatches(yA(), YHaplotype{id: "b", loci: {"DYS19": "14.1234567"}});
}

func testYAlleleValidationRaisesLineageErrors() {
    # Not a bare toFloat failure from `convert` - the caller catches "lineage".
    testing.assertThrows("stepDistanceOnANonNumericAllele", "lineage");
    testing.assertThrows("buildWithANonNumericAllele", "lineage");
    testing.assertThrows("buildWithAHugeAllele", "lineage");
    testing.assertThrows("stepDistanceOnASevenDigitAllele", "lineage");
    testing.assertThrows("matchOnALongFraction", "lineage");
    # six digits on each side still parse
    testing.assertEqual(len(yLoci(yHaplotype("ok", {"DYS19": "123456.123456"}))), 1);
}

func testYSingleCopyStepsAreExact() {
    # 18.2 - 17.2 = 1 exactly (the float difference is 0.99999...);
    # " 15 " against "14": trimmed, 1; "11" against "11.0": 0 steps.
    def a as YHaplotype init yHaplotype("A", {"X": "17.2", "Y": " 15 ", "Z": "11"});
    def b as YHaplotype init yHaplotype("B", {"X": "18.2", "Y": "14", "Z": "11.0"});
    testing.assertEqual(yStepDistance($a, $b), 2.0);
    testing.assertTrue(yMatchesWithin($a, $b, 2.0));
    testing.assertFalse(yMatches($a, $b));
    testing.assertTrue(yMatches(yHaplotype("C", {"Z": "11"}), yHaplotype("D", {"Z": "11.0"})));
    # 17.2 against 18 differs by a fraction: no step count
    testing.assertFalse(yMatchesWithin($a, yHaplotype("E", {"X": "18"}), 5.0));
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

# --- documented edge cases with no assertion yet ----------------------------

func confidenceOfOne() {
    haplotypeUpperBound(1, 10, 1.0);
}

func testConfidenceMustBeBelowOne() {
    testing.assertThrows("confidenceOfOne", "lineage");
}

func testZeroSightingBoundAtAnotherConfidence() {
    # k = 0 closed form at 99%: 1 - 0.01^(1/N).
    near(haplotypeUpperBound(0, 300, 0.99), 1.0 - math.pow(0.01, 1.0 / 300.0));
    near(lineageLikelihoodRatio(0, 100), 1.0 / (1.0 - math.pow(0.05, 1.0 / 100.0)));
}

func testParseTrimsAndUpperCases() {
    testing.assertEqual(parseDifference(" 16189c "), parseDifference("16189C"));
}

func badEmpty() {
    parseDifference("   ");
}

func badInsertionZero() {
    parseDifference("315.0C");
}

func testMoreMalformedDesignationsThrow() {
    testing.assertThrows("badEmpty", "lineage");
    testing.assertThrows("badInsertionZero", "lineage");
}

func formatUnknownKind() {
    formatDifference(Difference{position: 1, insertion: 0, base: "A", kind: "inv"});
}

func testFormatRejectsAnUnknownKind() {
    testing.assertThrows("formatUnknownKind", "lineage");
}

func testMtHaplotypeFromNumbersInReferenceCoordinates() {
    # Four reference bases from 16024 cover 16024-16027; the change is the third.
    def h as MtHaplotype init mtHaplotypeFrom("S1", "ACGT", "ACAT", 16024);
    testing.assertEqual($h.ranges, ["16024-16027"]);
    testing.assertEqual(mtDesignations($h), ["16026A"]);
}

func yEmptyLocusName() {
    yHaplotype("bad", {"": "14"});
}

func yNoAlleles() {
    yHaplotype("bad", {"DYS19": " , "});
}

func testYHaplotypeValidatesItsInput() {
    testing.assertThrows("yEmptyLocusName", "lineage");
    testing.assertThrows("yNoAlleles", "lineage");
}

func testYMultiCopyAcceptsASlash() {
    def b as YHaplotype init yHaplotype("B", {"DYS385": "14/11"});
    testing.assertTrue(yMatches(yA(), $b));
}

# --- regressions from the casework review -----------------------------------

# parkMiller is a Park-Miller generator private to these tests, so the
# sequences below do not depend on math.rand's algorithm.
func parkMiller(state as int) {
    return ($state * 16807) % 2147483647;
}

# seededSeq returns n bases drawn from `seed`, and the next state, as
# [sequence, state].
func seededSeq(n as int, seed as int) {
    def bases as list of string init ["A", "C", "G", "T"];
    def out as list of string init [];
    def state as int init $seed;
    for (def i in 0..$n) {
        $state = parkMiller($state);
        $out[] = $bases[$state % 4];
    }
    return [strings.join($out, ""), convert.toString($state)];
}

# fullMatrix is a plain O(nm) Needleman-Wunsch over the whole matrix, written
# independently of alignBand: full score matrix, no band, no byte tricks. Same
# scores and tie order (diagonal, then up, then left). Returns [top, bottom,
# score].
func fullMatrix(ref as string, smp as string) {
    def r as list of string init strings.chars($ref);
    def s as list of string init strings.chars($smp);
    def m as int init len($r);
    def n as int init len($s);
    def score as list of list of int init [];
    def move as list of list of int init [];
    for (def i in 0..($m + 1)) {
        def srow as list of int init [$i * GAP];
        def mrow as list of int init [DIR_UP];
        for (def j in 1..($n + 1)) {
            if ($i == 0) {
                $srow[] = $j * GAP;
                $mrow[] = DIR_LEFT;
                continue;
            }
            def d as int init $score[$i - 1][$j - 1] + MISMATCH;
            if ($r[$i - 1] == $s[$j - 1]) {
                $d = $score[$i - 1][$j - 1] + MATCH;
            }
            def u as int init $score[$i - 1][$j] + GAP;
            def l as int init $srow[$j - 1] + GAP;
            def best as int init $d;
            def mv as int init DIR_DIAG;
            if ($u > $best) {
                $best = $u;
                $mv = DIR_UP;
            }
            if ($l > $best) {
                $best = $l;
                $mv = DIR_LEFT;
            }
            $srow[] = $best;
            $mrow[] = $mv;
        }
        $score[] = $srow;
        $move[] = $mrow;
    }
    def top as list of string init [];
    def bottom as list of string init [];
    def i as int init $m;
    def j as int init $n;
    while ($i > 0 or $j > 0) {
        def mv as int init $move[$i][$j];
        if ($mv == DIR_DIAG) {
            $top[] = $r[$i - 1];
            $bottom[] = $s[$j - 1];
            $i = $i - 1;
            $j = $j - 1;
        } elseif ($mv == DIR_UP) {
            $top[] = $r[$i - 1];
            $bottom[] = "-";
            $i = $i - 1;
        } else {
            $top[] = "-";
            $bottom[] = $s[$j - 1];
            $j = $j - 1;
        }
    }
    return [
        strings.join(lists.reverse($top), ""),
        strings.join(lists.reverse($bottom), ""),
        convert.toString($score[$m][$n])
    ];
}

func testBandedAlignmentIsExactOnAShiftedWindow() {
    # Reference = 40 bases + core, sample = core + 40 bases: the same 200-base
    # core read through a window shifted by 40. The path the band sees first
    # (all diagonal, about a quarter matches) never touches the band edge, but
    # the optimum leaves it: 40 deletions, 200 matches, 40 insertions score
    # 200 - 80 * 2 = 40, against about 60 - 180 = -120 on the diagonal.
    def state as int init 7;
    for (def trial in 0..3) {
        def core as list of string init seededSeq(200, $state);
        def head as list of string init seededSeq(40, convert.toInt($core[1]));
        def tail as list of string init seededSeq(40, convert.toInt($head[1]));
        $state = convert.toInt($tail[1]);
        def ref as string init $head[0] + $core[0];
        def smp as string init $core[0] + $tail[0];
        def full as list of string init fullMatrix($ref, $smp);
        def banded as list of string init alignPair($ref, $smp);
        testing.assertTrue(convert.toInt($full[2]) >= 40);
        testing.assertEqual($banded[0], $full[0]);
        testing.assertEqual($banded[1], $full[1]);
    }
}

func testBandedAlignmentMatchesTheFullMatrix() {
    # The same differential as testBandedAlignmentMatchesFullWidth, against
    # the independent full-matrix aligner instead of alignBand at full width.
    math.randSeed(4242);
    for (def trial in 0..20) {
        def ref as string init randomSeq(math.randInt(20, 80));
        def smp as string init mutate(
            $ref,
            math.randInt(0, 8),
            math.randInt(0, 6),
            math.randInt(0, 6));
        def full as list of string init fullMatrix($ref, $smp);
        def banded as list of string init alignPair($ref, $smp);
        testing.assertEqual($banded[0], $full[0]);
        testing.assertEqual($banded[1], $full[1]);
    }
}

func testBandBoundIsTheBestScoreOutsideTheBand() {
    # m = n = 10, w = 2: leaving the band visits diagonal +-3 and returns to
    # 0, so at least 6 gaps, leaving (20 - 6) / 2 = 7 diagonal moves:
    # 7 * 1 + 6 * -2 = -5.
    testing.assertEqual(bandBound(10, 10, 2), -5);
    # m = 12, n = 10, w = 3: g = 2 * 4 - 2 = 6, (22 - 6) / 2 = 8: 8 - 12 = -4.
    testing.assertEqual(bandBound(12, 10, 3), -4);
}

func alignTooDivergent() {
    # Unrelated 1400-base sequences need the full 1400 x 1401 matrix, over
    # MAX_ALIGN_CELLS (1250000); it is refused before that pass runs.
    def a as list of string init seededSeq(1400, 11);
    def b as list of string init seededSeq(1400, convert.toInt($a[1]));
    align($a[0], $b[0], 1);
}

func testAlignRefusesABandPastTheCellCap() {
    testing.assertThrows("alignTooDivergent", "lineage");
}

func testGapShiftIsLinear() {
    # 5000 A against one A: one A kept at position 1 and 4999 deletions slid
    # 3' to 2..5000. The swap-until-stable shift took 22 s here.
    def t0 as int init time.unixMillis(time.now());
    def ds as list of Difference init align(strings.repeat("A", 5000), "A", 1);
    def ms as int init time.unixMillis(time.now()) - $t0;
    testing.assertEqual(len($ds), 4999);
    testing.assertEqual(formatDifference($ds[0]), "2d");
    testing.assertEqual(formatDifference($ds[4998]), "5000d");
    testing.assertTrue($ms < 1000);
}

func testRepeatIndelsMoveThreePrimeAsABlock() {
    # GACACT -> GACACACT: one more AC unit. Every placement of the inserted
    # pair scores the same; the 3'-most is after the last C (position 5).
    testing.assertEqual(designations(align("GACACT", "GACACACT", 1)), ["5.1A", "5.2C"]);
    # A CAG unit gained: G T T C A G C A T, inserted after position 8.
    testing.assertEqual(
        designations(align("GTTCAGCAT", "GTTCAGCAGCAT", 1)),
        ["8.1G", "8.2C", "8.3A"]);
    # An AC unit lost from GACACACT: the 3'-most AC is positions 6-7.
    testing.assertEqual(designations(align("GACACACT", "GACACT", 1)), ["6d", "7d"]);
}

func testLeadingInsertionRoundTrips() {
    # A G before the first reference base is numbered after position 0.
    def ds as list of Difference init align("ACGTACGT", "GACGTACGT", 1);
    testing.assertEqual(designations($ds), ["0.1G"]);
    testing.assertEqual(parseDifference("0.1G"), $ds[0]);
    testing.assertEqual(
        parseDifference("0.2a"),
        Difference{position: 0, insertion: 2, base: "A", kind: "ins"});
    def h as MtHaplotype init mtHaplotype("H", [], "73G 0.1G");
    testing.assertEqual(mtDesignations($h), ["0.1G", "73G"]);
}

func deletionAtZero() {
    parseDifference("0d");
}

func testPositionZeroIsOnlyAnInsertion() {
    testing.assertThrows("positionZero", "lineage");
    testing.assertThrows("deletionAtZero", "lineage");
}

func testDIsADeletionNotAnIupacCode() {
    testing.assertEqual(parseDifference("16189D").kind, "del");
    testing.assertEqual(parseDifference("16189B").base, "B");
}

func testRangesAgreeIgnoresDuplicates() {
    # ["x", "x"] has the same length as ["x", "y"] and every element is in it.
    def a as MtHaplotype init mtHaplotype("a", ["73-340", "73-340"], "263G");
    def b as MtHaplotype init mtHaplotype("b", ["73-340", "16024-16365"], "263G");
    def c as MtHaplotype init mtHaplotype("c", ["73-340"], "263G");
    testing.assertFalse(mtRangesAgree($a, $b));
    testing.assertFalse(mtRangesAgree($b, $a));
    testing.assertTrue(mtRangesAgree($a, $c));
    testing.assertTrue(mtRangesAgree($b, mtHaplotype("d", ["16024-16365", "73-340"], "")));
}

func repeatedDesignation() {
    mtHaplotype("c", ["73-340"], "263G 263G");
}

func conflictingDesignation() {
    mtHaplotype("e", ["73-340"], "263G 263A");
}

func deletedAndSubstituted() {
    mtHaplotype("f", ["73-340"], "263G 263d");
}

func testRepeatedOrConflictingPositionsThrow() {
    testing.assertThrows("repeatedDesignation", "lineage");
    testing.assertThrows("conflictingDesignation", "lineage");
    testing.assertThrows("deletedAndSubstituted", "lineage");
    # Distinct insertion indices at one position are distinct designations.
    def h as MtHaplotype init mtHaplotype("g", [], "315.1C 315.2C 315T");
    testing.assertEqual(mtDesignations($h), ["315T", "315.1C", "315.2C"]);
}

func testMatchesWithinIsFalseOnACopyCountMismatch() {
    # yMatches is false here; yMatchesWithin(0.0) is documented as yMatches.
    def y1 as YHaplotype init yHaplotype("y1", {"DYS385": "11,14", "DYS19": "14"});
    def y2 as YHaplotype init yHaplotype("y2", {"DYS385": "11", "DYS19": "14"});
    testing.assertFalse(yMatches($y1, $y2));
    testing.assertFalse(yMatchesWithin($y1, $y2, 0.0));
    testing.assertFalse(yMatchesWithin($y1, $y2, 100.0));
}

func microvariantAgainstWhole() {
    yStepDistance(yHaplotype("a", {"DYS458": "17.2"}), yHaplotype("b", {"DYS458": "18"}));
}

func trailingDot() {
    yHaplotype("z", {"DYS19": "11."});
}

func leadingDot() {
    yHaplotype("z", {"DYS19": ".5"});
}

func testMicrovariantsStepOnlyWithinTheirClass() {
    # 17.2 -> 18 is not 0.8 of a step: no number of whole-repeat steps joins
    # them. 17.2 -> 18.2 is exactly one step.
    testing.assertThrows("microvariantAgainstWhole", "lineage");
    def a as YHaplotype init yHaplotype("a", {"DYS458": "17.2", "DYS19": "14"});
    testing.assertEqual(yStepDistance($a, yHaplotype("b", {"DYS458": "18.2", "DYS19": "15"})), 2.0);
    testing.assertFalse(yMatchesWithin($a, yHaplotype("c", {"DYS458": "18"}), 5.0));
    testing.assertTrue(yMatchesWithin($a, yHaplotype("d", {"DYS458": "17.20"}), 0.0));
    testing.assertThrows("trailingDot", "lineage");
    testing.assertThrows("leadingDot", "lineage");
}
