# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# Forensic-efficiency table for a frequency database, and the statistics of an
# observed genotype sample.
#
# Run from the deck root:  jennifer run examples/popstats.j
#
# The frequency table and the sample are synthetic; see examples/lib/exampledb.j.
use io;

import "../src/forensicgenetics.j" as forensics;
import "../src/popstats.j" as popstats;
import "./lib/exampledb.j" as exampledb;

def db as forensics.FrequencyDb init exampledb.build();

# --- the table --------------------------------------------------------------

def table as list of popstats.LocusStats init popstats.statsTable($db);
io.printf("locus       He      MP      PD      PIC     TPI     PE trio PE duo\n");
for (def row in $table) {
    io.printf(
        "%s|pad=11 %f|prec=4 %f|prec=4 %f|prec=4 %f|prec=4 %f|prec=3 %f|prec=4 %f|prec=4\n",
        $row.locus,
        $row.expectedHeterozygosity,
        $row.matchProbability,
        $row.powerOfDiscrimination,
        $row.pic,
        $row.typicalPaternityIndex,
        $row.powerOfExclusionTrio,
        $row.powerOfExclusionDuo);
}

def all as popstats.CombinedStats init popstats.combinedStats($table);
io.printf("\ncombined over %d loci\n", $all.loci);
io.printf("  MP        10^%f|prec=2\n", $all.log10MatchProbability);
io.printf("  PD        %f|prec=10\n", $all.powerOfDiscrimination);
io.printf("  PE trio   %f|prec=6\n", $all.powerOfExclusionTrio);
io.printf("  PE duo    %f|prec=6\n", $all.powerOfExclusionDuo);
io.printf("  TPI       %f|prec=1\n", $all.typicalPaternityIndex);

# --- an observed sample -----------------------------------------------------

# 200 people at TH01 with a heterozygote deficit built in: 60 homozygous 6/6,
# 60 homozygous 8/8, the rest spread over four heterozygotes.
def sample as list of forensics.Genotype init [];
for (def i in 0..60) {
    $sample[] = forensics.genotype("TH01", "6", "6");
    $sample[] = forensics.genotype("TH01", "8", "8");
}
for (def i in 0..20) {
    $sample[] = forensics.genotype("TH01", "6", "8");
    $sample[] = forensics.genotype("TH01", "6", "9.3");
    $sample[] = forensics.genotype("TH01", "7", "8");
    $sample[] = forensics.genotype("TH01", "8", "9.3");
}

io.printf("\nTH01 sample, n = %d\n", len($sample));
io.printf("  observed He   %f|prec=4\n", popstats.observedHeterozygosity($sample));
def own as forensics.FrequencyDb init popstats.withEstimatedLocus(forensics.db("sample"), $sample);
io.printf("  expected He   %f|prec=4 (Nei 1978)\n", popstats.unbiasedHeterozygosity($own, "TH01"));
io.printf("  observed PE   %f|prec=4\n", popstats.observedPowerOfExclusion($sample));
def p as float init popstats.hwePermutationTest($sample, 1000, 1);
io.printf("  HWE p-value   %f|prec=4 (1000 permutations, seed 1)\n", $p);
