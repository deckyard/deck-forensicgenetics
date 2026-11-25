# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.28.0

# White-box test overlay for emp.j. Run with `jennifer test emp_test.j`.
#
# The files below are written by hand in the layout of EMPOP's template and
# published example files; no EMPOP data is copied.
use testing;
use strings;

# A small file: header comments, a control-region block, and a block typed
# over HVS-I and HVS-II.
func sample() {
    return strings.join(
        [
            "# Test population, 4 samples",
            "# contact: nobody@example.org",
            "#! 16024-576",
            "S1\tH1c\t1\t16519C\t263G\t315.1C\t477C",
            "S2\tC5b1\t1\t16093C\t16223T\t73G\t249DEL\t309.1C\t309.2C\t315.1C",
            "S3\tH2a2a1\t2\t",
            "#! 16024-16365 73-340",
            "S4\t?\t1\t16126C\t16294T\t73G\t263G\t315.1C"
        ],
        "\n") + "\n";
}

# --- reading ----------------------------------------------------------------

func testReadsCommentsRangesAndEntries() {
    def f as EmpFile init readEmp(sample());
    testing.assertEqual($f.comments, ["Test population, 4 samples", "contact: nobody@example.org"]);
    testing.assertEqual(len($f.entries), 4);
    def s1 as EmpEntry init $f.entries[0];
    testing.assertEqual($s1.haplotype.id, "S1");
    testing.assertEqual($s1.haplogroup, "H1c");
    testing.assertEqual($s1.count, 1);
    testing.assertEqual($s1.line, 4);
    testing.assertEqual($s1.haplotype.ranges, ["16024-576"]);
    testing.assertEqual(
        lineage.mtDesignations($s1.haplotype),
        ["263G", "315.1C", "477C", "16519C"]);
    testing.assertEqual($f.entries[3].haplotype.ranges, ["16024-16365", "73-340"]);
    testing.assertEqual($f.entries[3].haplogroup, "?");
}

func testReadsDeletionsAndRepeatedInsertions() {
    def f as EmpFile init readEmp(sample());
    testing.assertEqual(
        lineage.mtDesignations($f.entries[1].haplotype),
        ["73G", "249d", "309.1C", "309.2C", "315.1C", "16093C", "16223T"]);
}

func testAnEntryEqualToTheReferenceHasNoDifferences() {
    def f as EmpFile init readEmp(sample());
    testing.assertEqual(len($f.entries[2].haplotype.differences), 0);
    testing.assertEqual($f.entries[2].count, 2);
    # without the trailing empty cell as well
    $f = readEmp("#! 16024-576\nS1\tH2a2a1\t1\n");
    testing.assertEqual(len($f.entries[0].haplotype.differences), 0);
}

func testReadsCrlfABomAndBlankLines() {
    def text as string init "﻿# title\r\n\r\n#! 16024-576\r\nS1\tH\t1\t263G\t\r\n\r\n";
    def f as EmpFile init readEmp($text);
    testing.assertEqual($f.comments, ["title"]);
    testing.assertEqual(lineage.mtDesignations($f.entries[0].haplotype), ["263G"]);
}

func testCommentsBelowTheFirstRangeAreSkipped() {
    def text as string init "# head\n#! 16024-576\n# a note\n#S9\t?\t1\t263G\nS1\tH\t1\t263G\n";
    def f as EmpFile init readEmp($text);
    testing.assertEqual($f.comments, ["head"]);
    testing.assertEqual(len($f.entries), 1);
    testing.assertEqual($f.entries[0].haplotype.id, "S1");
}

func testCommentTextKeepsWhatFollowsOneSpace() {
    def f as EmpFile init readEmp("#\n#  indented\n#tight\n#! 1-10\n");
    testing.assertEqual($f.comments, ["", " indented", "tight"]);
    testing.assertEqual(len($f.entries), 0);
}

func testCountZeroIsKeptAndAddsNothing() {
    def f as EmpFile init readEmp("#! 16024-576\nS1\tH\t0\t263G\nS2\tH\t3\t263G\n");
    testing.assertEqual($f.entries[0].count, 0);
    testing.assertEqual(empSampleSize($f), 3);
}

func testDifferencesAcceptTheNotationsLineageAccepts() {
    def f as EmpFile init readEmp("#! 16024-576\nS1\tH\t1\ta263g\t249d\t16093y\n");
    testing.assertEqual(
        lineage.mtDesignations($f.entries[0].haplotype),
        ["249d", "263G", "16093Y"]);
}

func testAWrappingRangeCoversBothEnds() {
    def f as EmpFile init readEmp("#! 16024-576\nS1\tH\t1\t16569A\t1G\t576T\t16024C\n");
    testing.assertEqual(len($f.entries[0].haplotype.differences), 4);
}

# --- reading: refusals ------------------------------------------------------

func testRefusesAHaplotypeBeforeTheFirstRange() {
    testing.assertThrowsArgs("readEmp", ["# head\nS1\tH\t1\t263G\n"], "emp");
}

func testRefusesMalformedRanges() {
    testing.assertThrowsArgs("readEmp", ["#! SEQUENCE RANGE 1\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#!\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 16024-16570\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 0-576\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 16024\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 16024-576,73-340\n"], "emp");
}

func testRefusesSampleIdsEmpopRefuses() {
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS 1\tH\t1\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nSä1\tH\t1\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\n\tH\t1\t263G\n"], "emp");
    def f as EmpFile init readEmp("#! 1-576\nAUT-05/0038_b\tH\t1\t263G\n");
    testing.assertEqual($f.entries[0].haplotype.id, "AUT-05/0038_b");
}

func testRefusesARepeatedSampleId() {
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t263G\nS1\tH\t1\t73G\n"], "emp");
}

func testRefusesAMissingHaplogroupOrCount() {
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\t\t1\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1 H 1 263G\n"], "emp");
}

func testRefusesCountsOutsideTheRange() {
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t10000\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t-1\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\tone\t263G\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1.0\t263G\n"], "emp");
}

func testRefusesTwoVariantsInOneCell() {
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t263G 315.1C\n"], "emp");
}

func testRefusesWhatEmpcheckRefuses() {
    # not a designation
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t263Q\n"], "emp");
    # past the genome
    testing.assertThrowsArgs("readEmp", ["#! 16024-576\nS1\tH\t1\t16570A\n"], "emp");
    # outside the declared range
    testing.assertThrowsArgs("readEmp", ["#! 16024-576\nS1\tH\t1\t1000C\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 16024-16365 73-340\nS1\tH\t1\t16400C\n"], "emp");
    # a second insertion without a first
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t315.2C\n"], "emp");
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t309.1C\t309.3C\n"], "emp");
    # a position named twice
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t263G\t263A\n"], "emp");
    # an insertion before position 1, which lineage allows
    testing.assertThrowsArgs("readEmp", ["#! 1-576\nS1\tH\t1\t0.1A\n"], "emp");
}

func testARefusalNamesItsLine() {
    def msg as string init "";
    def kind as string init "";
    try {
        readEmp("# head\n#! 1-576\nS1\tH\t1\t263G\nS2\tH\t1\t263Q\n");
    } catch (e) {
        $kind = $e.kind;
        $msg = $e.message;
    }
    testing.assertEqual($kind, "emp");
    testing.assertTrue(strings.startsWith($msg, "line 4: "));
}

# --- building and writing ---------------------------------------------------

func testWritesTheSampleBack() {
    testing.assertEqual(writeEmp(readEmp(sample())), sample());
}

func testWritingAndReadingAgreeOnEveryEntry() {
    def f as EmpFile init readEmp(sample());
    def again as EmpFile init readEmp(writeEmp($f));
    testing.assertEqual($again, $f);
}

func testWritesDifferencesInRangeOrder() {
    def h as lineage.MtHaplotype init lineage.mtHaplotype(
        "S1",
        ["16024-576"],
        "73G 249d 16189C 315.1C");
    def f as EmpFile init EmpFile{comments: [], entries: [empEntry($h, "?", 1)]};
    testing.assertEqual(writeEmp($f), "#! 16024-576\nS1\t?\t1\t16189C\t73G\t249DEL\t315.1C\n");
}

func testWritesARangeLineOnlyWhereTheRangesChange() {
    def a as lineage.MtHaplotype init lineage.mtHaplotype("A", ["16024-576"], "263G");
    def b as lineage.MtHaplotype init lineage.mtHaplotype("B", ["16024-576"], "73G");
    def c as lineage.MtHaplotype init lineage.mtHaplotype("C", ["73-340"], "");
    def f as EmpFile init EmpFile{
        comments: ["t"],
        entries: [empEntry($a, "H", 1), empEntry($b, "H", 1), empEntry($c, "H", 1)]
    };
    testing.assertEqual(
        writeEmp($f),
        "# t\n#! 16024-576\nA\tH\t1\t263G\nB\tH\t1\t73G\n#! 73-340\nC\tH\t1\t\n");
}

func testEmpEntrySortsAndChecks() {
    def h as lineage.MtHaplotype init lineage.MtHaplotype{
        id: "S1",
        ranges: ["16024-576"],
        differences: [
            lineage.parseDifference("315.1C"),
            lineage.parseDifference("16189C"),
            lineage.parseDifference("73G")
        ]
    };
    def e as EmpEntry init empEntry($h, "H", 2);
    testing.assertEqual(lineage.mtDesignations($e.haplotype), ["73G", "315.1C", "16189C"]);
    testing.assertEqual($e.line, 0);
}

func badEntryRepeat() {
    def h as lineage.MtHaplotype init lineage.MtHaplotype{
        id: "S1",
        ranges: ["1-576"],
        differences: [lineage.parseDifference("263G"), lineage.parseDifference("263A")]
    };
    empEntry($h, "H", 1);
}

func badEntryId() {
    empEntry(lineage.mtHaplotype("S 1", ["1-576"], "263G"), "H", 1);
}

func badEntryCount() {
    empEntry(lineage.mtHaplotype("S1", ["1-576"], "263G"), "H", 10000);
}

func badEntryHaplogroup() {
    empEntry(lineage.mtHaplotype("S1", ["1-576"], "263G"), "H 1", 1);
}

func badEntryRange() {
    empEntry(lineage.mtHaplotype("S1", [], "263G"), "H", 1);
}

func testEmpEntryRefusesWhatReadEmpRefuses() {
    testing.assertThrows("badEntryRepeat", "emp");
    testing.assertThrows("badEntryId", "emp");
    testing.assertThrows("badEntryCount", "emp");
    testing.assertThrows("badEntryHaplogroup", "emp");
    testing.assertThrows("badEntryRange", "emp");
}

func writeTwice() {
    def e as EmpEntry init empEntry(lineage.mtHaplotype("S1", ["1-576"], "263G"), "H", 1);
    writeEmp(EmpFile{comments: [], entries: [$e, $e]});
}

func writeBrokenComment() {
    writeEmp(EmpFile{comments: ["a\nS1\tH\t1"], entries: []});
}

func writeUncheckedEntry() {
    def e as EmpEntry init EmpEntry{
        haplotype: lineage.mtHaplotype("S1", ["1-576"], "1000C"),
        haplogroup: "H",
        count: 1,
        line: 0
    };
    writeEmp(EmpFile{comments: [], entries: [$e]});
}

func testWriteEmpRefusesWhatReadEmpWouldRefuse() {
    testing.assertThrows("writeTwice", "emp");
    testing.assertThrows("writeBrokenComment", "emp");
    testing.assertThrows("writeUncheckedEntry", "emp");
}

func testEmpDesignation() {
    testing.assertEqual(empDesignation(lineage.parseDifference("16189c")), "16189C");
    testing.assertEqual(empDesignation(lineage.parseDifference("315.1C")), "315.1C");
    testing.assertEqual(empDesignation(lineage.parseDifference("249d")), "249DEL");
}

# --- counting ---------------------------------------------------------------

func testCountsMatchesOverTheSameRanges() {
    def text as string init "#! 16024-576\nA\tH\t2\t263G\t315.1C\nB\tH\t1\t263G\t315.1C\n" +
        "C\tH\t0\t263G\t315.1C\nD\tH\t1\t73G\n#! 16024-16365 73-340\nE\tH\t5\t263G\t315.1C\n";
    def f as EmpFile init readEmp($text);
    def h as lineage.MtHaplotype init lineage.mtHaplotype("q", ["16024-576"], "315.1C 263G");
    testing.assertEqual(empCount($f, $h), 3);
    testing.assertEqual(empSampleSize($f), 9);
    # over other ranges, only the entry typed over those ranges counts
    def hv as lineage.MtHaplotype init lineage.mtHaplotype(
        "q",
        ["73-340", "16024-16365"],
        "263G 315.1C");
    testing.assertEqual(empCount($f, $hv), 5);
    testing.assertEqual(empCount($f, lineage.mtHaplotype("q", ["16024-576"], "16189C")), 0);
}

func testCountsFeedTheLineageEstimates() {
    def f as EmpFile init readEmp(sample());
    def h as lineage.MtHaplotype init lineage.mtHaplotype("q", ["16024-576"], "");
    testing.assertEqual(empCount($f, $h), 2);
    testing.assertEqual(
        lineage.lineageMatchProbability(empCount($f, $h), empSampleSize($f)),
        lineage.haplotypeUpperBound(2, 5, lineage.DEFAULT_CONFIDENCE));
}
