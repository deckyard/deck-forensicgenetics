# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Weir & Cockerham F-statistics over three population samples, and pairwise
# genotypic linkage-disequilibrium tests within one of them.
#
# Run from the deck root:  jennifer run examples/popstructure.j
#
# The samples are synthetic: TH01 differs between the populations, D8S1179,
# FGA and VWA do not, and FGA and VWA genotypes were made to travel together.
use io;
use convert;

import "../src/forensicgenetics.j" as forensics;
import "../src/popstats.j" as popstats;

# The other three loci follow fixed cycles shared by every population, so they
# do not differ between populations. Person i gets the i mod 4 TH01 genotype,
# D8S1179 by (i div 4) mod 3, and FGA and VWA by i mod 5 - the same index for
# both, which ties them together.
def const D8 as list of list of string init [["13", "13"], ["13", "14"], ["14", "14"]];
def const FGA as list of list of string init [
    ["21", "21"],
    ["21", "22"],
    ["22", "22"],
    ["22", "24"],
    ["21", "24"]
];
def const VWA as list of list of string init [
    ["16", "16"],
    ["16", "17"],
    ["17", "17"],
    ["17", "18"],
    ["16", "18"]
];

# population builds `count` synthetic profiles with TH01 drawn in turn from
# `th01`.
func population(prefix as string, count as int, th01 as list of list of string) {
    def out as list of forensics.Profile init [];
    for (def i in 0..$count) {
        def t as list of string init $th01[$i % len($th01)];
        def d as list of string init D8[($i // 4) % 3];
        def f as list of string init FGA[$i % 5];
        def v as list of string init VWA[$i % 5];
        $out[] = forensics.profile(
            $prefix + "{$i}",
            [
                forensics.genotype("TH01", $t[0], $t[1]),
                forensics.genotype("D8S1179", $d[0], $d[1]),
                forensics.genotype("FGA", $f[0], $f[1]),
                forensics.genotype("VWA", $v[0], $v[1])
            ]);
    }
    return $out;
}

def popA as list of forensics.Profile init population(
    "A",
    120,
    [["6", "6"], ["6", "9.3"], ["9.3", "9.3"], ["7", "9.3"]]);
def popB as list of forensics.Profile init population(
    "B",
    100,
    [["9.3", "9.3"], ["8", "9.3"], ["9.3", "9.3"], ["6", "9.3"]]);
def popC as list of forensics.Profile init population(
    "C",
    80,
    [["6", "7"], ["7", "7"], ["6", "6"], ["6", "7"]]);

# --- F-statistics -----------------------------------------------------------

def loci as list of string init ["TH01", "D8S1179", "FGA", "VWA"];
def table as list of popstats.FStatistics init popstats.fStatisticsTable(
    [$popA, $popB, $popC],
    $loci);
io.printf("locus       theta     F         f\n");
for (def row in $table) {
    io.printf(
        "%s|pad=11 %f|prec=4 %f|prec=4 %f|prec=4\n",
        $row.locus,
        popstats.fst($row),
        popstats.fit($row),
        popstats.fis($row));
}
def all as popstats.FStatistics init popstats.combineFStatistics($table);
io.printf(
    "%s|pad=11 %f|prec=4 %f|prec=4 %f|prec=4\n",
    "all loci",
    popstats.fst($all),
    popstats.fit($all),
    popstats.fis($all));

def sampleTh01 as list of forensics.Genotype init [];
for (def i in 0..len($popA)) {
    $sampleTh01[] = $popA[$i].loci["TH01"];
}
io.printf(
    "\npopulation A, TH01 alone: f = %f|prec=4\n",
    popstats.inbreedingCoefficient($sampleTh01));

# --- linkage disequilibrium -------------------------------------------------

# The cycles make every pair but FGA-VWA exactly balanced, the most probable
# table, so those pairs get p = 1.
def pairs as list of popstats.LdTest init popstats.ldPermutationTests($popA, 1000, 7);
io.printf("\npopulation A, %d locus pairs (1000 permutations, seed 7)\n", len($pairs));
io.printf("Bonferroni threshold at 0.05: %f|prec=4\n", 0.05 / convert.toFloat(len($pairs)));
for (def t in $pairs) {
    io.printf("  %s|pad=8 %s|pad=8 n %d  p %f|prec=4\n", $t.first, $t.second, $t.n, $t.pValue);
}
