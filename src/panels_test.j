# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# White-box test overlay for panels.j. Run with `jennifer test panels_test.j`.
#
# Panel expectations restate the published locus lists independently of the
# tables in panels.j: the set relations the sources state (CODIS 20 = CODIS 13
# + 7, GlobalFiler = CODIS 20 + SE33, Fusion = CODIS 20 + Penta D/E,
# Fusion 6C = Fusion + SE33, ESX 17 = ESS + D2S1338/D16S539/D19S433/SE33,
# Identifiler = CODIS 13 + D2S1338/D19S433) and the locus counts the manuals
# give.
use testing;
use lists;
use math;

# The 13 original CODIS loci, written out again from Budowle et al. (1998).
func codis13Again() {
    return [
        "CSF1PO",
        "FGA",
        "TH01",
        "TPOX",
        "vWA",
        "D3S1358",
        "D5S818",
        "D7S820",
        "D8S1179",
        "D13S317",
        "D16S539",
        "D18S51",
        "D21S11"
    ];
}

# The seven loci the 2009 ESS shares with the 1990s European set.
func oldEss() {
    return ["D3S1358", "vWA", "D8S1179", "D21S11", "D18S51", "TH01", "FGA"];
}

# sameSet compares two locus lists ignoring order, and rejects duplicates.
func sameSet(xs as list of string, ys as list of string) {
    testing.assertEqual(len($xs), len($ys));
    for (def x in $xs) {
        testing.assertTrue(lists.contains($ys, $x));
    }
}

# minus is xs without ys, order kept.
func minus(xs as list of string, ys as list of string) {
    def out as list of string init [];
    for (def x in $xs) {
        if (not lists.contains($ys, $x)) {
            $out[] = $x;
        }
    }
    return $out;
}

# --- catalogue --------------------------------------------------------------

func testPanelNamesListsEveryBuiltIn() {
    def names as list of string init panelNames();
    testing.assertEqual(len($names), 11);
    for (def n in $names) {
        testing.assertTrue(hasPanel($n));
        testing.assertEqual(panel($n).name, $n);
    }
}

func testLocusCountsMatchTheManuals() {
    testing.assertEqual(len(panel("CODIS13").loci), 13);
    testing.assertEqual(len(panel("CODIS20").loci), 20);
    testing.assertEqual(len(panel("ESS").loci), 12);
    testing.assertEqual(len(panel("GlobalFiler").loci), 21);
    testing.assertEqual(len(panel("PowerPlex Fusion").loci), 22);
    testing.assertEqual(len(panel("PowerPlex Fusion 6C").loci), 23);
    testing.assertEqual(len(panel("PowerPlex ESX 17").loci), 16);
    testing.assertEqual(len(panel("NGM SElect").loci), 16);
    testing.assertEqual(len(panel("Investigator ESSplex SE QS").loci), 16);
    testing.assertEqual(len(panel("Identifiler Plus").loci), 15);
    testing.assertEqual(len(panel("Investigator IDplex Plus").loci), 15);
}

func testCodisSetsAsPublished() {
    sameSet(panel("CODIS13").loci, codis13Again());
    def added as list of string init minus(panel("CODIS20").loci, codis13Again());
    sameSet($added, ["D1S1656", "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045"]);
}

func testEssAsPublished() {
    def added as list of string init minus(panel("ESS").loci, oldEss());
    sameSet($added, ["D1S1656", "D2S441", "D10S1248", "D12S391", "D22S1045"]);
}

func testKitsAsSupersetsOfTheCoreSets() {
    def codis20 as list of string init panel("CODIS20").loci;
    def ess as list of string init panel("ESS").loci;
    sameSet(minus(panel("GlobalFiler").loci, $codis20), ["SE33"]);
    sameSet(minus(panel("PowerPlex Fusion").loci, $codis20), ["Penta D", "Penta E"]);
    sameSet(minus(panel("PowerPlex Fusion 6C").loci, panel("PowerPlex Fusion").loci), ["SE33"]);
    def essPlus as list of string init ["D2S1338", "D16S539", "D19S433", "SE33"];
    sameSet(minus(panel("PowerPlex ESX 17").loci, $ess), $essPlus);
    sameSet(minus(panel("NGM SElect").loci, $ess), $essPlus);
    sameSet(minus(panel("Investigator ESSplex SE QS").loci, $ess), $essPlus);
    sameSet(minus(panel("Identifiler Plus").loci, codis13Again()), ["D2S1338", "D19S433"]);
    sameSet(panel("Investigator IDplex Plus").loci, panel("Identifiler Plus").loci);
}

func testPanelsHoldOnlyAutosomalCanonicalNames() {
    for (def n in panelNames()) {
        def p as Panel init panel($n);
        def seen as list of string init [];
        for (def locus in $p.loci) {
            testing.assertFalse(lists.contains($seen, $locus));
            testing.assertEqual(normalizeLocus($locus), $locus);
            testing.assertFalse($locus == "AMEL");
            testing.assertFalse(strings.startsWith($locus, "DYS"));
            $seen[] = $locus;
        }
        testing.assertTrue(len($p.source) > 0);
    }
    testing.assertTrue(lists.contains(panel("GlobalFiler").otherMarkers, "DYS391"));
    testing.assertTrue(lists.contains(panel("PowerPlex Fusion 6C").otherMarkers, "DYS576"));
}

func testPanelLookupFoldsCaseAndSeparators() {
    testing.assertEqual(panel("globalfiler").name, "GlobalFiler");
    testing.assertEqual(panel("PowerPlex-ESX-17").name, "PowerPlex ESX 17");
    testing.assertEqual(panel("powerplex_fusion_6c").name, "PowerPlex Fusion 6C");
    testing.assertFalse(hasPanel("PowerPlex 21"));
}

func unknownPanel() {
    panel("Profiler Plus");
}

func testUnknownPanelThrows() {
    testing.assertThrows("unknownPanel", "panels");
}

# --- intersections and custom panels ----------------------------------------

func testSharedLociOfIdentifilerAndEssIsTheOldEuropeanSet() {
    def got as list of string init sharedLoci(panel("Identifiler Plus"), panel("ESS"));
    sameSet($got, oldEss());
    # In the first panel's order: Identifiler lists D3S1358 before TH01.
    testing.assertTrue(lists.indexOf($got, "D3S1358") < lists.indexOf($got, "TH01"));
}

func testSharedLociOfFusionAndEsx17() {
    # ESX 17's sixteen less SE33, which Fusion does not type.
    def got as list of string init sharedLoci(panel("PowerPlex Fusion"), panel("PowerPlex ESX 17"));
    sameSet($got, minus(panel("PowerPlex ESX 17").loci, ["SE33"]));
}

func testIntersectIsAPanel() {
    def p as Panel init intersect(panel("GlobalFiler"), panel("ESS"));
    testing.assertEqual($p.name, "GlobalFiler & ESS");
    sameSet($p.loci, panel("ESS").loci);
}

func disjointIntersect() {
    intersect(customPanel("a", ["D3S1358"]), customPanel("b", ["TH01"]));
}

func testDisjointIntersectThrows() {
    testing.assertThrows("disjointIntersect", "panels");
}

func testCustomPanelNormalises() {
    def p as Panel init customPanel("mine", ["THO1", "vwa", "penta_e"]);
    testing.assertEqual($p.loci, ["TH01", "vWA", "Penta E"]);
    testing.assertEqual(len($p.otherMarkers), 0);
}

func customDuplicate() {
    customPanel("dup", ["TH01", "THO1"]);
}

func customEmpty() {
    customPanel("none", []);
}

func customUnnamed() {
    customPanel("", ["TH01"]);
}

func testCustomPanelValidation() {
    testing.assertThrows("customDuplicate", "panels");
    testing.assertThrows("customEmpty", "panels");
    testing.assertThrows("customUnnamed", "panels");
}

# --- locus names ------------------------------------------------------------

func testNormalizeLocusVariants() {
    testing.assertEqual(normalizeLocus("THO1"), "TH01");
    testing.assertEqual(normalizeLocus("th01"), "TH01");
    testing.assertEqual(normalizeLocus("HUMTH01"), "TH01");
    testing.assertEqual(normalizeLocus("VWA"), "vWA");
    testing.assertEqual(normalizeLocus("vWA"), "vWA");
    testing.assertEqual(normalizeLocus("HUMVWFA31/A"), "vWA");
    testing.assertEqual(normalizeLocus("VWFA31A"), "vWA");
    testing.assertEqual(normalizeLocus("Penta E"), "Penta E");
    testing.assertEqual(normalizeLocus("PentaE"), "Penta E");
    testing.assertEqual(normalizeLocus("PENTA_E"), "Penta E");
    testing.assertEqual(normalizeLocus("penta-d"), "Penta D");
    testing.assertEqual(normalizeLocus("d8s1179"), "D8S1179");
    testing.assertEqual(normalizeLocus(" D2S441 "), "D2S441");
    testing.assertEqual(normalizeLocus("dys391"), "DYS391");
    testing.assertEqual(normalizeLocus("CSF1P0"), "CSF1PO");
    testing.assertEqual(normalizeLocus("HUMFIBRA"), "FGA");
    testing.assertEqual(normalizeLocus("ACTBP2"), "SE33");
    testing.assertEqual(normalizeLocus("Amelogenin"), "AMEL");
    # An unfamiliar name passes through, trimmed, case kept.
    testing.assertEqual(normalizeLocus(" L "), "L");
    testing.assertEqual(normalizeLocus("Locus7"), "Locus7");
}

func testSegmentNameShape() {
    testing.assertTrue(isSegmentName("D21S11"));
    testing.assertTrue(isSegmentName("DXS10135"));
    testing.assertFalse(isSegmentName("DS11"));
    testing.assertFalse(isSegmentName("D21S"));
    testing.assertFalse(isSegmentName("D2AS11"));
    testing.assertTrue(isSegmentName("DYS391"));
    testing.assertFalse(isSegmentName("D1S2S3"));
    # not folded: the fast path in normalizeLocus must not take it
    testing.assertFalse(isSegmentName("d8s1179"));
    testing.assertFalse(isSegmentName("D8S1179 "));
}

func testCanonicalNamesAreFixedPointsOfTheFoldedPath() {
    # normalizeLocus returns a CANONICAL or upper-case D-segment name without
    # folding. Lower-casing forces the folded path, which must give the same
    # name back: for every CANONICAL entry and every built-in panel locus.
    for (def name in CANONICAL) {
        testing.assertEqual(normalizeLocus(strings.lower($name)), $name);
    }
    for (def pn in panelNames()) {
        for (def locus in panel($pn).loci) {
            testing.assertEqual(normalizeLocus(strings.lower($locus)), $locus);
            testing.assertEqual(normalizeLocus($locus), $locus);
        }
    }
}

func blankLocus() {
    normalizeLocus("  ");
}

func testBlankLocusThrows() {
    testing.assertThrows("blankLocus", "panels");
}

func testNormalizeProfileRenamesKeysAndGenotypes() {
    def p as forensics.Profile init forensics.profile(
        "P",
        [forensics.genotype("THO1", "6", "9.3"), forensics.genotype("vwa", "16", "17")]);
    def n as forensics.Profile init normalizeProfile($p);
    testing.assertEqual(forensics.profileLoci($n), ["TH01", "vWA"]);
    testing.assertEqual(forensics.genotypeAt($n, "TH01").locus, "TH01");
    testing.assertEqual(forensics.genotypeAt($n, "TH01").b, "9.3");
    testing.assertEqual($n.id, "P");
}

func profileCollision() {
    normalizeProfile(forensics.profile(
        "P",
        [forensics.genotype("TH01", "6", "6"), forensics.genotype("THO1", "7", "7")]));
}

func testNormalizeProfileCollisionThrows() {
    testing.assertThrows("profileCollision", "panels");
}

func messyDb() {
    def out as forensics.FrequencyDb init forensics.withLocus(
        forensics.db("messy"),
        "THO1",
        500,
        {"6": 0.25, "9.3": 0.75});
    return forensics.withLocus($out, "d8s1179", 0, {"12": 0.5, "13": 0.5});
}

func testNormalizeDbRenamesAndKeepsFrequencies() {
    def n as forensics.FrequencyDb init normalizeDb(messyDb());
    testing.assertEqual(forensics.loci($n), ["TH01", "D8S1179"]);
    testing.assertEqual($n.name, "messy");
    testing.assertEqual($n.loci["TH01"].locus, "TH01");
    testing.assertEqual(forensics.sampleSize($n, "TH01"), 500);
    testing.assertEqual(forensics.alleleFrequency($n, "TH01", "9.3"), 0.75);
}

func dbCollision() {
    normalizeDb(forensics.withLocus(messyDb(), "TH01", 10, {"6": 1.0}));
}

func testNormalizeDbCollisionThrows() {
    testing.assertThrows("dbCollision", "panels");
}

# --- narrowing --------------------------------------------------------------

func testRestrictProfileKeepsPanelLociInPanelOrder() {
    def p as forensics.Profile init forensics.profile(
        "P",
        [
            forensics.genotype("THO1", "6", "9.3"),
            forensics.genotype("D6S1043", "11", "12"),
            forensics.genotype("D8S1179", "12", "13")
        ]);
    def r as forensics.Profile init restrictProfile($p, panel("ESS"));
    # ESS order puts D8S1179 before TH01; D6S1043 is no ESS locus.
    testing.assertEqual(forensics.profileLoci($r), ["D8S1179", "TH01"]);
    testing.assertEqual(forensics.genotypeAt($r, "TH01").locus, "TH01");
    def gaps as list of string init missingInProfile($p, panel("ESS"));
    testing.assertEqual(len($gaps), 10);
    testing.assertFalse(lists.contains($gaps, "TH01"));
    testing.assertTrue(lists.contains($gaps, "vWA"));
}

func testRestrictedProfileScoresAsTheCoreRestrict() {
    # Cross-check: narrowing through a panel gives the same RMP as narrowing
    # by hand with forensicgenetics.restrict over the same loci.
    def fdb as forensics.FrequencyDb init normalizeDb(messyDb());
    def p as forensics.Profile init forensics.profile(
        "P",
        [
            forensics.genotype("TH01", "6", "9.3"),
            forensics.genotype("D8S1179", "12", "13"),
            forensics.genotype("D6S1043", "11", "12")
        ]);
    def mine as Panel init customPanel("two", ["TH01", "D8S1179"]);
    def viaPanel as float init forensics.randomMatchProbability(
        $fdb,
        restrictProfile($p, $mine),
        0.0);
    def viaCore as float init forensics.randomMatchProbability(
        $fdb,
        forensics.restrict($p, ["TH01", "D8S1179"]),
        0.0);
    # By hand, theta 0: TH01 6/9.3 is 2(.25)(.75) = .375; D8S1179 12/13 is
    # 2(.5)(.5) = .5; the product is .1875.
    testing.assertTrue(math.abs($viaPanel - 0.1875) < 0.000000001);
    testing.assertTrue(math.abs($viaPanel - $viaCore) < 0.000000001);
}

func testRestrictDbAndMissingInDb() {
    def r as forensics.FrequencyDb init restrictDb(messyDb(), panel("ESS"));
    testing.assertEqual(forensics.loci($r), ["D8S1179", "TH01"]);
    testing.assertEqual(len(missingInDb(messyDb(), panel("ESS"))), 10);
    testing.assertEqual(missingInDb(messyDb(), customPanel("x", ["TH01", "D8S1179"])), []);
}

# --- regressions from the casework review -----------------------------------

func missingWithProfileCollision() {
    # TH01 and THO1 are one locus; missingInProfile used to keep the last and
    # report the panel covered.
    def p as forensics.Profile init forensics.profile(
        "P",
        [forensics.genotype("TH01", "6", "7"), forensics.genotype("THO1", "8", "9")]);
    missingInProfile($p, customPanel("x", ["TH01"]));
}

func missingWithDbCollision() {
    def fdb as forensics.FrequencyDb init forensics.withLocus(
        forensics.withLocus(forensics.db("t"), "TH01", 100, {"6": 0.5, "7": 0.5}),
        "THO1",
        100,
        {"8": 0.5, "9": 0.5});
    missingInDb($fdb, customPanel("x", ["TH01"]));
}

func testMissingInThrowsOnACollisionLikeRestrict() {
    testing.assertThrows("missingWithProfileCollision", "panels");
    testing.assertThrows("missingWithDbCollision", "panels");
}
