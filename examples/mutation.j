# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Familias mutation models and a silent allele on the trio of
# examples/paternity.j: the D8S1179 inconsistency under each model, the
# stationary (PM) variant, a maternal model, published per-locus rates, and an
# opposite-homozygote duo explained by a silent allele.
#
# Run from the deck root:  jennifer run examples/mutation.j
#
# The frequency table is synthetic; see examples/lib/exampledb.j.
use io;
use math;

import "../src/forensicgenetics.j" as forensics;
import "./lib/exampledb.j" as exampledb;
import "../src/kinship.j" as kinship;
import "../src/mutation.j" as mutation;

def db as forensics.FrequencyDb init exampledb.build();

def mother as forensics.Profile init forensics.profile(
    "mother",
    [
        forensics.genotype("D3S1358", "15", "16"),
        forensics.genotype("vWA", "16", "17"),
        forensics.genotype("FGA", "21", "22"),
        forensics.genotype("D8S1179", "12", "13"),
        forensics.genotype("D21S11", "29", "30"),
        forensics.genotype("D18S51", "14", "15"),
        forensics.genotype("TH01", "6", "7"),
        forensics.genotype("D16S539", "11", "12")
    ]);

def child as forensics.Profile init forensics.profile(
    "child",
    [
        forensics.genotype("D3S1358", "15", "17"),
        forensics.genotype("vWA", "17", "19"),
        forensics.genotype("FGA", "22", "24"),
        forensics.genotype("D8S1179", "13", "15"),
        forensics.genotype("D21S11", "30", "31.2"),
        forensics.genotype("D18S51", "15", "17"),
        forensics.genotype("TH01", "7", "9.3"),
        forensics.genotype("D16S539", "12", "13")
    ]);

# A 16/16 man cannot pass on the child's 15 at D8S1179 without a mutation.
def alleged as forensics.Profile init forensics.profile(
    "alleged father",
    [
        forensics.genotype("D3S1358", "17", "18"),
        forensics.genotype("vWA", "19", "19"),
        forensics.genotype("FGA", "24", "25"),
        forensics.genotype("D8S1179", "16", "16"),
        forensics.genotype("D21S11", "31.2", "32.2"),
        forensics.genotype("D18S51", "17", "18"),
        forensics.genotype("TH01", "9.3", "10"),
        forensics.genotype("D16S539", "13", "14")
    ]);

def m as forensics.Genotype init forensics.genotypeAt($mother, "D8S1179");
def c as forensics.Genotype init forensics.genotypeAt($child, "D8S1179");
def f as forensics.Genotype init forensics.genotypeAt($alleged, "D8S1179");

# One paternal rate, five models. The deck's own stepwise ladder for contrast.
def rate as float init 0.002;
def names as list of string init [
    "Equal",
    "Proportional",
    "Stepwise (range .1)",
    "Extended stepwise",
    "Extended stepwise, PM"
];
def specs as list of mutation.MutationSpec init [
    mutation.equalModel($rate),
    mutation.proportionalModel($rate),
    mutation.stepwiseModel($rate, 0.1),
    mutation.extendedStepwiseModel($rate, 0.1, 0.0001),
    mutation.stationaryModel(mutation.extendedStepwiseModel($rate, 0.1, 0.0001))
];
io.printf(
    "D8S1179  mother %s/%s  child %s/%s  alleged father %s/%s\n\n",
    $m.a,
    $m.b,
    $c.a,
    $c.b,
    $f.a,
    $f.b);
io.printf("paternal model                PI at D8S1179   CPI\n");
for (def i in 0..len($specs)) {
    def pm as mutation.ParentageModel init mutation.parentageModel(
        $specs[$i],
        mutation.NO_MUTATION);
    io.printf(
        "%s|pad=28  %f|prec=7       %f|prec=1\n",
        $names[$i],
        mutation.paternityIndexAt($db, $m, $c, $f, $pm, 0.0),
        mutation.paternityIndex($mother, $child, $alleged, $db, $pm, 0.0));
}
io.printf(
    "%s|pad=28  %f|prec=7\n",
    "kinship.j stepwise ladder",
    kinship.paternityIndexAt($db, $m, $c, $f, kinship.mutationModel($rate, 0.1)));

# The extended stepwise matrix is not stationary under the table's frequencies; PM fixes
# that and keeps the expected rate.
def la as mutation.LocusAlleles init mutation.locusAlleles($db, "D8S1179", [], 0.0);
for (def i in [3, 4]) {
    def mm as mutation.MutationMatrix init mutation.mutationMatrix($specs[$i], $la);
    io.printf(
        "\n%s|pad=28  stationary %t  expected rate %f|prec=6",
        $names[$i],
        mutation.isStationary($mm, $la),
        mutation.overallRate($mm, $la));
}
io.printf("\n");

# Published per-locus rates (AABB 2003 via STRBase): paternal rates for the
# father, a lower maternal model for the mother.
def paternal as mutation.MutationSpec init mutation.withLocusRates(
    mutation.stepwiseModel($rate, 0.1),
    mutation.aabb2003Rates("paternal"));
def maternal as mutation.MutationSpec init mutation.withLocusRates(
    mutation.stepwiseModel(0.0005, 0.1),
    mutation.aabb2003Rates("maternal"));
def sexed as mutation.ParentageModel init mutation.parentageModel($paternal, $maternal);
io.printf(
    "\nAABB 2003 paternal rate at D8S1179 = %f|prec=5;  CPI = %f|prec=1\n",
    mutation.aabb2003Rates("paternal")["D8S1179"],
    mutation.paternityIndex($mother, $child, $alleged, $db, $sexed, 0.0));

# Opposite homozygotes in a duo: excluded outright unless an allele can be
# silent.
def cg as forensics.Genotype init forensics.genotype("TH01", "7", "7");
def fg as forensics.Genotype init forensics.genotype("TH01", "9.3", "9.3");
io.printf("\nduo at TH01, child 7, alleged father 9.3 (single peaks)\n");
for (def q in [0.0, 0.001, 0.01]) {
    def pm as mutation.ParentageModel init mutation.withSilentFrequency(
        mutation.MENDELIAN,
        "TH01",
        $q);
    io.printf(
        "silent frequency %f|prec=3  PI %f|prec=5\n",
        $q,
        mutation.duoPaternityIndexAt($db, $cg, $fg, $pm, 0.0));
}
io.printf(
    "log10 of the parent/child LR with q = .01: %f|prec=3\n",
    math.log10(
        mutation.kinshipLRAt($db, $cg, $fg, kinship.PARENT_CHILD, kinship.UNRELATED, 0.01, 0.0)));
