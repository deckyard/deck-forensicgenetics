# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# A stain searched against a small reference database: locus names
# normalised, narrowed to the loci the frequency table covers, a direct search
# allowing drop-out, and a familial search for a full sibling.
#
# Run from the deck root:  jennifer run examples/search.j
#
# The frequency table and every profile are synthetic; see
# examples/lib/exampledb.j.
use io;
use convert;

import "../src/forensicgenetics.j" as forensics;
import "./lib/exampledb.j" as exampledb;
import "../src/kinship.j" as kinship;
import "../src/panels.j" as panels;
import "../src/search.j" as search;

def db as forensics.FrequencyDb init exampledb.build();
def loci as list of string init forensics.loci($db);

# The stain, as another lab's software wrote it: "VWA", "THO1", and a ninth
# locus the table does not cover. TH01 shows 6/6 where the donor is 6/9.3.
def stain as forensics.Profile init forensics.profile(
    "stain",
    [
        forensics.genotype("D3S1358", "15", "16"),
        forensics.genotype("VWA", "17", "19"),
        forensics.genotype("FGA", "22", "24"),
        forensics.genotype("D8S1179", "13", "15"),
        forensics.genotype("D21S11", "30", "31.2"),
        forensics.genotype("D18S51", "15", "16"),
        forensics.genotype("THO1", "6", "6"),
        forensics.genotype("D16S539", "11", "13"),
        forensics.genotype("D2S1338", "17", "20")
    ]);

io.printf(
    "ESS loci the table lacks: %s\n",
    convert.toString(panels.missingInDb($db, panels.panel("ESS"))));
def table as panels.Panel init panels.customPanel("table", $loci);
def query as forensics.Profile init panels.restrictProfile($stain, $table);
io.printf("stain narrowed to %d loci\n\n", len($query.loci));

# Synthetic references: a linear congruential generator draws alleles from the
# table, so the run is reproducible.
def seed as int init 20261105;

func draw(locus as string) {
    def labels as list of string init [];
    for (def a in $db.loci[$locus].alleles) {
        $labels[] = $a;
    }
    $seed = ($seed * 1103515245 + 12345) % 2147483648;
    return $labels[($seed // 65536) % len($labels)];
}

func randomPerson(id as string) {
    def gs as list of forensics.Genotype init [];
    for (def locus in $loci) {
        $gs[] = forensics.genotype($locus, draw($locus), draw($locus));
    }
    return forensics.profile($id, $gs);
}

# A full sibling of the donor: at each locus one allele from the donor, the
# other drawn at random (IBD 1 throughout, a simplification).
func sibling(id as string, donor as forensics.Profile) {
    def gs as list of forensics.Genotype init [];
    for (def locus in $loci) {
        $gs[] = forensics.genotype($locus, $donor.loci[$locus].b, draw($locus));
    }
    return forensics.profile($id, $gs);
}

def donor as forensics.Profile init forensics.profile(
    "R-0042",
    [
        forensics.genotype("D3S1358", "15", "16"),
        forensics.genotype("vWA", "17", "19"),
        forensics.genotype("FGA", "22", "24"),
        forensics.genotype("D8S1179", "13", "15"),
        forensics.genotype("D21S11", "30", "31.2"),
        forensics.genotype("D18S51", "15", "16"),
        forensics.genotype("TH01", "6", "9.3"),
        forensics.genotype("D16S539", "11", "13")
    ]);

def database as list of forensics.Profile init [];
for (def i in 0..200) {
    if ($i == 42) {
        $database[] = $donor;
    } elseif ($i == 117) {
        $database[] = sibling("R-0117", $donor);
    } else {
        $database[] = randomPerson("R-{$i}");
    }
}

def opts as search.SearchOptions init search.searchOptions(8, 1, true);
def ranked as list of search.Comparison init search.search($query, $database, $opts);
io.printf("direct search, %d references\n", len($database));
io.printf(
    "%s|pad=8 %s|pad=13 %s|pad=9 %s|pad=10 %s|pad=9 %s\n",
    "id",
    "verdict",
    "compared",
    "identical",
    "drop-out",
    "shared");
for (def i in 0..3) {
    def c as search.Comparison init $ranked[$i];
    io.printf(
        "%s|pad=8 %s|pad=13 %d|pad=9 %d|pad=10 %d|pad=9 %d\n",
        $c.id,
        $c.verdict,
        $c.lociCompared,
        $c.matchingLoci,
        $c.dropoutLoci,
        $c.sharedAlleles);
}
io.printf("hits: %d\n\n", len(search.hits($ranked)));

# The donor tops the list (a person is maximally "sibling-like" to himself).
# With eight loci the planted sibling, R-0117, ranks below unrelated
# references that share alleles by chance: familial searching needs many loci.
def kin as list of search.FamilialHit init search.familialSearch(
    $query,
    $database,
    kinship.FULL_SIBLINGS,
    $db,
    8);
io.printf("familial search, full siblings vs unrelated\n");
for (def i in 0..len($kin)) {
    if ($i < 3 or $kin[$i].id == "R-0117") {
        io.printf("%d|pad=4 %s|pad=8 LR %f|prec=1\n", $i + 1, $kin[$i].id, $kin[$i].lr);
    }
}
