# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Single-source match probability and a mixture's CPI / CPE.
#
# Run from the deck root:  jennifer run examples/matchprobability.j
#
# The frequency table is synthetic; see examples/lib/exampledb.j.
use io;
use math;

import "../src/forensicgenetics.j" as forensics;
import "./lib/exampledb.j" as exampledb;

def db as forensics.FrequencyDb init exampledb.build();

io.printf("reference table: %s, %d loci\n\n", $db.name, len(forensics.loci($db)));

# --- a single-source stain --------------------------------------------------

def suspect as forensics.Profile init forensics.profile(
    "suspect",
    [
        forensics.genotype("D3S1358", "15", "16"),
        forensics.genotype("vWA", "17", "17"),
        forensics.genotype("FGA", "21", "24"),
        forensics.genotype("D8S1179", "12", "14"),
        forensics.genotype("D21S11", "29", "31.2"),
        forensics.genotype("D18S51", "14", "17"),
        forensics.genotype("TH01", "6", "9.3"),
        forensics.genotype("D16S539", "11", "12")
    ]);

io.printf("locus       genotype    frequency\n");
for (def locus in forensics.profileLoci($suspect)) {
    def g as forensics.Genotype init forensics.genotypeAt($suspect, $locus);
    def f as float init forensics.genotypeFrequency($db, $g, forensics.THETA_GENERAL);
    io.printf("%s|pad=12 %s|pad=12 %f|prec=6\n", $locus, $g.a + "/" + $g.b, $f);
}

# Theta is the sub-population correction: a larger theta is more conservative.
io.printf("\ntheta   log10 RMP   the profile is about 1 in\n");
for (def theta in [forensics.THETA_NONE, forensics.THETA_GENERAL, forensics.THETA_ISOLATED]) {
    def logRmp as float init forensics.log10RandomMatchProbability($db, $suspect, $theta);
    def lr as float init forensics.matchLikelihoodRatio($db, $suspect, $theta);
    io.printf("%f|prec=2 %f|prec=3 %d\n", $theta, $logRmp, math.round($lr));
}

# --- a two-person mixture ---------------------------------------------------

def stain as forensics.Mixture init forensics.mixture(
    "stain",
    {
        "D3S1358": ["15", "16", "17"],
        "vWA": ["16", "17"],
        "FGA": ["21", "22", "24"],
        "TH01": ["6", "7", "9.3"]
    });

io.printf("\nmixture over %d loci\n", len(forensics.mixtureLoci($stain)));
for (def locus in forensics.mixtureLoci($stain)) {
    def pi as float init forensics.inclusionProbability($db, $locus, $stain.loci[$locus]);
    io.printf("  %s|pad=10 PI = %f|prec=4  PE = %f|prec=4\n", $locus, $pi, 1.0 - $pi);
}

def cpi as float init forensics.cpi($db, $stain);
io.printf("\nCPI = %f|prec=6  - 1 in %d could not be excluded\n", $cpi, math.round(1.0 / $cpi));
io.printf("CPE = %f|prec=6  - the share of the population excluded\n", forensics.cpe($db, $stain));

if (forensics.isIncluded($suspect, $stain)) {
    io.printf("\nthe suspect cannot be excluded as a contributor\n");
} else {
    io.printf("\nthe suspect is excluded as a contributor\n");
}
