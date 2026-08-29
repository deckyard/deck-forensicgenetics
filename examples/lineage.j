# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# Lineage markers: rendering a sequence into rCRS-style difference
# nomenclature, the conservative counting estimate for a haplotype, and Y-STR
# haplotype comparison.
#
# Run from the deck root:  jennifer run examples/lineage.j
#
# The "reference" below is a short synthetic stand-in, NOT the rCRS. No
# reference sequence is bundled with this deck - pass the real region you
# sequenced (NC_012920.1) to `align` for real coordinates.
use io;

import "../src/lineage.j" as lineage;

# --- rendering a sequence into differences ----------------------------------

# A synthetic 40 bp "region", numbered as if it began at position 300, with a
# run of four Cs at 311-314 - the kind of homopolymer where indel placement
# actually matters.
def reference as string init "GTACGTTACGCCCCTAGGCATTACGGATCCAGTTACGATC";
def sample as string init "GTACGTTACGCCCCCTAGGCATTACGGATTCAGTTACGATC";

def diffs as list of lineage.Difference init lineage.align($reference, $sample, 300);

io.printf("reference : %s\n", $reference);
io.printf("sample    : %s\n\n", $sample);
io.printf("differences from the reference:\n");
for (def d in $diffs) {
    io.printf("  %s|pad=10 %s at %d\n", lineage.formatDifference($d), $d.kind, $d.position);
}

# The extra C lands at the 3' end of the homopolymer run, which is the
# convention that makes an extra C in the rCRS 311-315 stretch read as 315.1C.
def haplotype as lineage.MtHaplotype init lineage.mtHaplotypeFrom(
    "evidence",
    $reference,
    $sample,
    300);
io.printf(
    "\nhaplotype %s over %v: %v\n",
    $haplotype.id,
    $haplotype.ranges,
    lineage.mtDesignations($haplotype));

# --- what a database sighting is worth --------------------------------------

# Lineage markers do not multiply across loci: the whole haplotype is counted
# once against a reference database. The deck supplies the estimator; the
# count comes from whichever database was queried.
io.printf("\nk       N       k/N          95%% upper bound   LR\n");
for (def k in [0, 1, 2, 5]) {
    def n as int init 5000;
    def point as float init lineage.haplotypeFrequency($k, $n);
    def bound as float init lineage.lineageMatchProbability($k, $n);
    def lr as float init lineage.lineageLikelihoodRatio($k, $n);
    io.printf("%d|pad=8 %d|pad=8 %f|prec=8 %f|prec=8 %f|prec=1\n", $k, $n, $point, $bound, $lr);
}

io.printf("\nNote the k = 0 row: the point estimate is 0, which no database can\n");
io.printf("support. The exact-binomial upper bound is what a report quotes.\n");

# --- Y-STR haplotypes -------------------------------------------------------

def suspect as lineage.YHaplotype init lineage.yHaplotype(
    "suspect",
    {
        "DYS19": "14",
        "DYS385": "11,14",
        "DYS389I": "13",
        "DYS390": "24",
        "DYS391": "11",
        "DYS392": "13",
        "DYS393": "13"
    });

def relative as lineage.YHaplotype init lineage.yHaplotype(
    "paternal relative",
    {
        "DYS19": "14",
        "DYS385": "14,11",
        "DYS389I": "13",
        "DYS390": "25",
        "DYS391": "11",
        "DYS392": "13",
        "DYS393": "13"
    });

def unrelated as lineage.YHaplotype init lineage.yHaplotype(
    "unrelated man",
    {
        "DYS19": "15",
        "DYS385": "11,15",
        "DYS389I": "14",
        "DYS390": "23",
        "DYS391": "10",
        "DYS392": "11",
        "DYS393": "14"
    });

io.printf("\nY-STR comparison over %d loci\n", len(lineage.yLoci($suspect)));
for (def other in [$relative, $unrelated]) {
    def steps as float init lineage.yStepDistance($suspect, $other);
    io.printf(
        "  vs %s|pad=20 exact = %t  steps = %f|prec=1  within one step = %t\n",
        $other.id,
        lineage.yMatches($suspect, $other),
        $steps,
        lineage.yMatchesWithin($suspect, $other, 1.0));
}

io.printf("\nDYS385 is multi-copy: \"11,14\" and \"14,11\" are the same haplotype.\n");
io.printf("The relative differs by one repeat at DYS390 - a single-step match,\n");
io.printf("which is what keeps a paternal lineage together across one mutation.\n");
