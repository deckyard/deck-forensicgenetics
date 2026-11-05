# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for search.j. Run with `jennifer test search_test.j`.
#
# Direct-search expectations are counted by hand from the genotypes written
# beside them. Familial LRs use the kinship overlay's table - loci "L" and "M",
# alleles 10/11/12/13 at .1/.2/.3/.4, N = 10000 (floor below every allele) -
# and the textbook closed forms, and are cross-checked against
# kinship.kinshipLR candidate by candidate.
use testing;
use math;

def const TOL as float init 0.000000001;

func near(actual as float, expected as float) {
    def scale as float init math.max(1.0, math.abs($expected));
    testing.assertTrue(math.abs($actual - $expected) < TOL * $scale);
}

func gl(locus as string, a as string, b as string) {
    return forensics.genotype($locus, $a, $b);
}

func fixture() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("test"),
        "L",
        10000,
        {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
    return forensics.withLocus($out, "M", 10000, {"10": 0.1, "11": 0.2, "12": 0.3, "13": 0.4});
}

# The query: L 10/11, M 12/13, K 10/10.
func query() {
    return forensics.profile("Q", [gl("L", "10", "11"), gl("M", "12", "13"), gl("K", "10", "10")]);
}

# Six candidates; the comment on each gives what the query sees there.
func database() {
    return [
        # 0: no allele shared anywhere - 3 compared, 3 mismatches, IBS 0+0+0.
        forensics.profile("mis", [gl("L", "12", "13"), gl("M", "10", "11"), gl("K", "11", "12")]),
        # 1: identical at all 3 - IBS 2+2+2 = 6.
        forensics.profile("exact", [gl("L", "10", "11"), gl("M", "12", "13"), gl("K", "10", "10")]),
        # 2: K is query 10/10 vs 10/12 - drop-out-consistent, IBS 2+2+1 = 5.
        forensics.profile("drop", [gl("L", "10", "11"), gl("M", "12", "13"), gl("K", "10", "12")]),
        # 3: one shared locus only.
        forensics.profile("few", [gl("L", "10", "11")]),
        # 4: M is query 12/13 vs 12/12 - the query shows 13, which the
        # candidate lacks: a mismatch even with drop-out allowed. IBS 2+1+2.
        forensics.profile(
            "onemis",
            [gl("L", "10", "11"), gl("M", "12", "12"), gl("K", "10", "10")]),
        # 5: identical at the 2 loci it has - IBS 2+2 = 4.
        forensics.profile("exact2", [gl("L", "10", "11"), gl("M", "12", "13")])
    ];
}

func ids(rs as list of Comparison) {
    def out as list of string init [];
    for (def r in $rs) {
        $out[] = $r.id;
    }
    return $out;
}

# --- allele sharing ---------------------------------------------------------

func testAllelesSharedCountsTheMultisetIntersection() {
    testing.assertEqual(allelesShared(gl("L", "10", "11"), gl("L", "10", "11")), 2);
    testing.assertEqual(allelesShared(gl("L", "10", "10"), gl("L", "10", "10")), 2);
    testing.assertEqual(allelesShared(gl("L", "10", "10"), gl("L", "10", "11")), 1);
    testing.assertEqual(allelesShared(gl("L", "10", "11"), gl("L", "10", "10")), 1);
    testing.assertEqual(allelesShared(gl("L", "10", "11"), gl("L", "11", "12")), 1);
    testing.assertEqual(allelesShared(gl("L", "10", "11"), gl("L", "12", "13")), 0);
    # A literal out of canonical order still shares both.
    def rev as forensics.Genotype init forensics.Genotype{locus: "L", a: "11", b: "10"};
    testing.assertEqual(allelesShared(gl("L", "10", "11"), $rev), 2);
}

func sharedAcrossLoci() {
    allelesShared(gl("L", "10", "11"), gl("M", "10", "11"));
}

func testAllelesSharedRejectsDifferentLoci() {
    testing.assertThrows("sharedAcrossLoci", "search");
}

# --- options ----------------------------------------------------------------

func zeroMinLoci() {
    searchOptions(0, 0, false);
}

func negativeMismatches() {
    searchOptions(5, -1, false);
}

func badLiteral() {
    def bad as SearchOptions init SearchOptions{minLoci: 0, maxMismatches: 0, allowDropout: true};
    search(query(), database(), $bad);
}

func testOptionsAreValidated() {
    def o as SearchOptions init searchOptions(8, 2, true);
    testing.assertEqual($o.minLoci, 8);
    testing.assertEqual($o.maxMismatches, 2);
    testing.assertTrue($o.allowDropout);
    testing.assertThrows("zeroMinLoci", "search");
    testing.assertThrows("negativeMismatches", "search");
    testing.assertThrows("badLiteral", "search");
}

# --- one comparison ---------------------------------------------------------

func testCompareProfilesCounts() {
    def c as Comparison init compareProfiles(query(), database()[4], searchOptions(2, 1, true));
    testing.assertEqual($c.id, "onemis");
    testing.assertEqual($c.index, 0);
    testing.assertEqual($c.lociCompared, 3);
    testing.assertEqual($c.matchingLoci, 2);
    testing.assertEqual($c.dropoutLoci, 0);
    testing.assertEqual($c.mismatchedLoci, ["M"]);
    testing.assertEqual($c.ibs, {"L": 2, "M": 1, "K": 2});
    testing.assertEqual($c.ibs0, 0);
    testing.assertEqual($c.ibs1, 1);
    testing.assertEqual($c.ibs2, 2);
    testing.assertEqual($c.sharedAlleles, 5);
    testing.assertEqual($c.verdict, VERDICT_PARTIAL);
}

func testDropoutIsDirectionalAndOptional() {
    # Query K 10/10 vs candidate 10/12.
    def withDrop as Comparison init compareProfiles(
        query(),
        database()[2],
        searchOptions(2, 0, true));
    testing.assertEqual($withDrop.dropoutLoci, 1);
    testing.assertEqual(len($withDrop.mismatchedLoci), 0);
    testing.assertEqual($withDrop.verdict, VERDICT_PARTIAL);
    def without as Comparison init compareProfiles(
        query(),
        database()[2],
        searchOptions(2, 0, false));
    testing.assertEqual($without.dropoutLoci, 0);
    testing.assertEqual($without.mismatchedLoci, ["K"]);
    testing.assertEqual($without.verdict, VERDICT_MISMATCH);
    # Swapped roles: the query 10/12 shows 12, the candidate 10/10 lacks it.
    def swapped as Comparison init compareProfiles(
        database()[2],
        query(),
        searchOptions(2, 0, true));
    testing.assertEqual($swapped.mismatchedLoci, ["K"]);
}

func testVerdicts() {
    def o as SearchOptions init searchOptions(2, 1, true);
    testing.assertEqual(compareProfiles(query(), database()[0], $o).verdict, VERDICT_MISMATCH);
    testing.assertEqual(compareProfiles(query(), database()[1], $o).verdict, VERDICT_EXACT);
    testing.assertEqual(compareProfiles(query(), database()[3], $o).verdict, VERDICT_INSUFFICIENT);
    testing.assertEqual(compareProfiles(query(), database()[5], $o).verdict, VERDICT_EXACT);
    testing.assertEqual(VERDICT_INSUFFICIENT, "insufficient");
}

func emptyQuery() {
    compareProfiles(forensics.profile("E", []), query(), searchOptions(1, 0, false));
}

func emptyQuerySearch() {
    search(forensics.profile("E", []), database(), searchOptions(1, 0, false));
}

func testEmptyQueryThrows() {
    testing.assertThrows("emptyQuery", "search");
    testing.assertThrows("emptyQuerySearch", "search");
}

# --- ranking ----------------------------------------------------------------

func testSearchRanksByVerdictThenCloseness() {
    # exact (3 identical) > exact2 (2 identical); drop (0 mismatches) >
    # onemis (1); then mis; then few (insufficient).
    def rs as list of Comparison init search(query(), database(), searchOptions(2, 1, true));
    testing.assertEqual(ids($rs), ["exact", "exact2", "drop", "onemis", "mis", "few"]);
    testing.assertEqual($rs[0].index, 1);
    testing.assertEqual($rs[5].index, 3);
    testing.assertEqual(ids(hits($rs)), ["exact", "exact2", "drop", "onemis"]);
}

func testSearchTiesKeepDatabaseOrder() {
    # Without drop-out, "drop" and "onemis" both have 1 mismatch, 2 identical
    # loci and 5 shared alleles: database order (2 before 4) decides.
    def rs as list of Comparison init search(query(), database(), searchOptions(2, 1, false));
    testing.assertEqual(ids($rs), ["exact", "exact2", "drop", "onemis", "mis", "few"]);
    # A stricter mismatch allowance demotes both.
    def strict as list of Comparison init search(query(), database(), searchOptions(2, 0, false));
    testing.assertEqual(ids(hits($strict)), ["exact", "exact2"]);
}

func testSearchOfAnEmptyDatabase() {
    testing.assertEqual(len(search(query(), [], searchOptions(1, 0, false))), 0);
}

# --- familial search --------------------------------------------------------

func kinQuery() {
    return forensics.profile("Q", [gl("L", "10", "11"), gl("M", "10", "10")]);
}

func kinDatabase() {
    return [
        # 0: one allele shared at each locus.
        forensics.profile("a", [gl("L", "10", "12"), gl("M", "10", "11")]),
        # 1: identical heterozygote at L, identical homozygote at M.
        forensics.profile("b", [gl("L", "10", "11"), gl("M", "10", "10")]),
        # 2: nothing shared.
        forensics.profile("c", [gl("L", "12", "13"), gl("M", "12", "13")]),
        # 3: typed at L only.
        forensics.profile("d", [gl("L", "10", "11")]),
        # 4: a's genotype at M again, so the memoised value is reused.
        forensics.profile("e", [gl("L", "11", "13"), gl("M", "10", "11")])
    ];
}

func testFamilialSearchMatchesClosedForms() {
    def hs as list of FamilialHit init familialSearch(
        kinQuery(),
        kinDatabase(),
        kinship.FULL_SIBLINGS,
        fixture(),
        2);
    # Candidate b, full siblings, closed forms (docs/project/testing.md):
    # L identical heterozygotes 10/11: 1/4 + (1 + .1 + .2)/(8 * .1 * .2)
    #   = .25 + 1.3/.16 = 8.375;
    # M identical homozygotes 10/10: (1 + 1/.1)^2 / 4 = 121/4 = 30.25;
    # LR = 8.375 * 30.25 = 253.34375.
    testing.assertEqual($hs[0].id, "b");
    near($hs[0].lr, 253.34375);
    testing.assertEqual($hs[0].lociCompared, 2);
    testing.assertTrue($hs[0].scored);
    # d shares one locus, below minLoci 2: unscored and last.
    testing.assertEqual($hs[4].id, "d");
    testing.assertFalse($hs[4].scored);
    testing.assertEqual($hs[4].lociCompared, 1);
    near($hs[4].lr, 0.0);
    # Ranked by LR, descending, over the scored candidates.
    for (def i in 1..4) {
        testing.assertTrue($hs[$i - 1].lr >= $hs[$i].lr);
    }
}

func testFamilialSearchParentChild() {
    # Parent-child, one locus L: query 10/11 vs 10/12 shares only allele 10.
    # Joint under PC = P(10/11) * 1/2 * P(other allele is 12) = .04 * .5 * .3
    # = .006; unrelated = .04 * (2 * .1 * .3) = .0024; LR = 2.5 = 1/(4 p10).
    def q as forensics.Profile init forensics.profile("Q", [gl("L", "10", "11")]);
    def db as list of forensics.Profile init [
        forensics.profile("x", [gl("L", "12", "13")]),
        forensics.profile("y", [gl("L", "10", "12")])
    ];
    def hs as list of FamilialHit init familialSearch($q, $db, kinship.PARENT_CHILD, fixture(), 1);
    testing.assertEqual($hs[0].id, "y");
    testing.assertEqual($hs[0].index, 1);
    near($hs[0].lr, 2.5);
    # x shares no allele: excluded as a parent, LR 0, ranked after y.
    near($hs[1].lr, 0.0);
    testing.assertTrue($hs[1].scored);
}

func testFamilialSearchEqualsKinshipLR() {
    # The memoised per-locus product must reproduce kinship.kinshipLR exactly,
    # for every relationship and candidate.
    def rels as list of kinship.Relationship init [
        kinship.FULL_SIBLINGS,
        kinship.HALF_SIBLINGS,
        kinship.FIRST_COUSINS
    ];
    def db as list of forensics.Profile init kinDatabase();
    for (def rel in $rels) {
        def hs as list of FamilialHit init familialSearch(kinQuery(), $db, $rel, fixture(), 1);
        for (def h in $hs) {
            testing.assertEqual(
                $h.lr,
                kinship.kinshipLR(kinQuery(), $db[$h.index], $rel, fixture()));
        }
    }
}

func familialUncovered() {
    def q as forensics.Profile init forensics.profile(
        "Q",
        [gl("L", "10", "11"), gl("Z", "1", "2")]);
    def db as list of forensics.Profile init [
        forensics.profile("x", [gl("L", "10", "11"), gl("Z", "1", "2")])
    ];
    familialSearch($q, $db, kinship.FULL_SIBLINGS, fixture(), 1);
}

func familialZeroMin() {
    familialSearch(kinQuery(), kinDatabase(), kinship.FULL_SIBLINGS, fixture(), 0);
}

func familialEmptyQuery() {
    familialSearch(forensics.profile("E", []), kinDatabase(), kinship.FULL_SIBLINGS, fixture(), 1);
}

func testFamilialSearchErrors() {
    testing.assertThrows("familialUncovered", "kinship");
    testing.assertThrows("familialZeroMin", "search");
    testing.assertThrows("familialEmptyQuery", "search");
}
