# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# A mother / child / alleged-father trio, including the single locus where the
# child carries an allele the alleged father cannot have transmitted unchanged.
#
# Run from the deck root:  jennifer run examples/paternity.j
#
# The frequency table is synthetic; see examples/lib/exampledb.j.
use io;
use math;

import "../src/forensicgenetics.j" as forensics;
import "./lib/exampledb.j" as exampledb;
import "../src/kinship.j" as kinship;

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

# Consistent everywhere except D8S1179, where a 16/16 man cannot transmit the
# child's 15 without a one-repeat mutation.
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

def loci as list of string init kinship.trioLoci($db, $mother, $child, $alleged);
io.printf("trio over %d loci\n\n", len($loci));

io.printf("locus        mother  child    father     PI strict  PI mutation\n");
for (def locus in $loci) {
    def mg as forensics.Genotype init forensics.genotypeAt($mother, $locus);
    def cg as forensics.Genotype init forensics.genotypeAt($child, $locus);
    def fg as forensics.Genotype init forensics.genotypeAt($alleged, $locus);
    def strict as float init kinship.paternityIndexAt($db, $mg, $cg, $fg, kinship.NO_MUTATION);
    def relaxed as float init kinship.paternityIndexAt(
        $db,
        $mg,
        $cg,
        $fg,
        kinship.DEFAULT_MUTATION);
    io.printf(
        "%s|pad=12 %s|pad=8 %s|pad=8 %s|pad=10 %f|prec=4 %f|prec=4\n",
        $locus,
        $mg.a + "/" + $mg.b,
        $cg.a + "/" + $cg.b,
        $fg.a + "/" + $fg.b,
        $strict,
        $relaxed);
}

def bad as list of string init kinship.inconsistentLoci($mother, $child, $alleged, $db);
io.printf("\nMendelian-inconsistent loci: %v\n", $bad);

# Without a mutation model one inconsistency is an outright exclusion.
def strict as float init kinship.paternityIndexWith(
    $mother,
    $child,
    $alleged,
    $db,
    kinship.NO_MUTATION);
io.printf("\nstrict Mendelian CPI  = %f|prec=6\n", $strict);
if ($strict <= 0.0) {
    io.printf("  -> excluded: no allowance is made for mutation\n");
}

# With one, the trio survives - at a cost proportional to how unlikely the
# mutation is.
def cpi as float init kinship.paternityIndex($mother, $child, $alleged, $db);
def w as float init kinship.probabilityOfPaternity($cpi);
io.printf("\nCPI with mutation     = %f|prec=2\n", $cpi);
io.printf("log10 CPI             = %f|prec=3\n", math.log10($cpi));
io.printf("W (prior 0.5)         = %f|prec=8\n", $w);
io.printf("W (prior 0.1)         = %f|prec=8\n", kinship.probabilityOfPaternityPrior($cpi, 0.1));

# The same evidence read as a sibling question, for contrast.
io.printf(
    "\nfor contrast, child vs alleged father as full siblings: LR = %f|prec=4\n",
    kinship.kinshipLR($child, $alleged, kinship.FULL_SIBLINGS, $db));
