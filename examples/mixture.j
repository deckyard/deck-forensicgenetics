# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# Semi-continuous mixture LR: victim + suspect against victim + unknown, with
# drop-out and drop-in. With d = 0 the suspect is excluded at D3S1358 (his 18
# is missing) and log10 LR is undefined.
#
# Run from the deck root:  jennifer run examples/mixture.j
#
# The frequency table is synthetic; see examples/lib/exampledb.j.
use io;

import "../src/forensicgenetics.j" as forensics;
import "../src/mixture.j" as mixture;
import "./lib/exampledb.j" as exampledb;

def db as forensics.FrequencyDb init exampledb.build();

def victim as forensics.Profile init forensics.profile(
    "victim",
    [
        forensics.genotype("D3S1358", "15", "17"),
        forensics.genotype("vWA", "16", "16"),
        forensics.genotype("TH01", "6", "7"),
        forensics.genotype("D16S539", "9", "11")
    ]);
def suspect as forensics.Profile init forensics.profile(
    "suspect",
    [
        forensics.genotype("D3S1358", "16", "18"),
        forensics.genotype("vWA", "17", "19"),
        forensics.genotype("TH01", "9.3", "9.3"),
        forensics.genotype("D16S539", "11", "13")
    ]);

# D3S1358 lost the suspect's 18 (drop-out); D16S539 shows a 14 nobody carries
# (drop-in, or a third person).
def stain as forensics.Mixture init forensics.mixture(
    "stain",
    {
        "D3S1358": ["15", "16", "17"],
        "vWA": ["16", "17", "19"],
        "TH01": ["6", "7", "9.3"],
        "D16S539": ["9", "11", "13", "14"]
    });

io.printf("stain needs at least %d contributors\n\n", mixture.minimumContributors($stain));

def theta as float init forensics.THETA_GENERAL;
def dropIn as float init 0.05;

io.printf("drop-out   log10 LR\n");
for (def d in [0.05, 0.1, 0.3]) {
    def hp as mixture.Hypothesis init mixture.hypothesis(
        "victim + suspect",
        [$victim, $suspect],
        [0.0, $d],
        0,
        $d);
    def hd as mixture.Hypothesis init mixture.hypothesis(
        "victim + unknown",
        [$victim],
        [0.0],
        1,
        $d);
    io.printf(
        "%f|prec=2     %f|prec=3\n",
        $d,
        mixture.log10MixtureLR($db, $stain, $hp, $hd, $dropIn, $theta));
}

def hp as mixture.Hypothesis init mixture.hypothesis(
    "victim + suspect",
    [$victim, $suspect],
    [0.0, 0.1],
    0,
    0.1);
def hd as mixture.Hypothesis init mixture.hypothesis("victim + unknown", [$victim], [0.0], 1, 0.1);
def per as map of string to float init mixture.locusLRs($db, $stain, $hp, $hd, $dropIn, $theta);
io.printf("\nper locus, d = 0.1\n");
for (def locus in $per) {
    io.printf("  %s|pad=10 LR = %f|prec=2\n", $locus, $per[$locus]);
}
