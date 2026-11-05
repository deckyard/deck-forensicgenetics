# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# X-STRs in a deficiency case: are two women paternal half sisters? Then a
# daughter's paternity with and without her mother, and the same pair scored
# with the Argus X-12 linkage group LG1 as one haplotype.
#
# Run from the deck root:  jennifer run examples/xstr.j
#
# SYNTHETIC frequencies and haplotypes: made-up numbers for illustration, not
# population data. Real X-STR haplotype tables come from population studies.
use io;

import "../src/forensicgenetics.j" as forensics;
import "../src/kinship.j" as kinship;
import "../src/xstr.j" as xstr;

func table() {
    def db as forensics.FrequencyDb init forensics.db("Synthetic X-STR Population");
    $db = forensics.withLocus(
        $db,
        "DXS8378",
        500,
        {"10": 0.30, "11": 0.35, "12": 0.25, "13": 0.10});
    $db = forensics.withLocus(
        $db,
        "DXS10135",
        500,
        {"18": 0.10, "19": 0.15, "20": 0.20, "21": 0.25, "22": 0.30});
    $db = forensics.withLocus(
        $db,
        "DXS10148",
        500,
        {"22.1": 0.25, "23.1": 0.20, "24.1": 0.30, "25.1": 0.25});
    $db = forensics.withLocus(
        $db,
        "DXS7132",
        500,
        {"12": 0.10, "13": 0.35, "14": 0.35, "15": 0.20});
    return $db;
}

func woman(id as string, gs as list of list of string) {
    def loci as list of string init ["DXS8378", "DXS10135", "DXS10148", "DXS7132"];
    def out as list of forensics.Genotype init [];
    for (def i in 0..len($loci)) {
        $out[] = forensics.genotype($loci[$i], $gs[$i][0], $gs[$i][1]);
    }
    return xstr.femaleProfile($id, $out);
}

def independent as xstr.XFrequencies init xstr.independentModel(table());

def alice as xstr.XProfile init woman(
    "Alice",
    [["10", "13"], ["18", "21"], ["23.1", "24.1"], ["12", "14"]]);
def beth as xstr.XProfile init woman(
    "Beth",
    [["11", "13"], ["18", "22"], ["23.1", "25.1"], ["12", "15"]]);

io.printf("Paternal half sisters, X-STRs, independent loci\n");
def phs as float init xstr.kinshipLR($alice, $beth, xstr.PATERNAL_HALF_SISTERS, $independent);
def mhs as float init xstr.kinshipLR($alice, $beth, xstr.MATERNAL_HALF_SISTERS, $independent);
io.printf("  LR paternal half sisters vs unrelated: %f|prec=2\n", $phs);
io.printf("  LR maternal half sisters vs unrelated: %f|prec=2\n", $mhs);
io.printf("  autosomal STRs cannot tell these two apart; X-STRs can\n");

# Daughter paternity: the father passes his whole X, so the daughter's
# paternal allele must be his at every locus.
def mother as xstr.XProfile init woman(
    "mother",
    [["10", "12"], ["18", "21"], ["23.1", "24.1"], ["12", "14"]]);
def daughter as xstr.XProfile init woman(
    "daughter",
    [["11", "12"], ["21", "22"], ["24.1", "25.1"], ["13", "14"]]);
def af as xstr.XProfile init xstr.maleProfile(
    "alleged father",
    {"DXS8378": "11", "DXS10135": "22", "DXS10148": "25.1", "DXS7132": "13"});
def trio as float init xstr.daughterTrioIndex(
    $mother,
    $daughter,
    $af,
    $independent,
    kinship.NO_MUTATION);
def duo as float init xstr.daughterDuoIndex($daughter, $af, $independent, kinship.NO_MUTATION);
io.printf("Daughter paternity, independent loci\n");
io.printf("  trio PI (mother typed):   %f|prec=2\n", $trio);
io.printf("  duo PI (mother untyped):  %f|prec=2\n", $duo);

# LG1 as one marker, from a (synthetic) haplotype table; DXS7132 stays on its
# own. Unlisted haplotypes score at the 0.002 floor.
def lg1 as list of string init xstr.argusX12Groups()["LG1"];
def haps as map of string to float init {
    "10-19-22.1": 0.12,
    "11-21-24.1": 0.10,
    "11-22-25.1": 0.08,
    "10-21-24.1": 0.05,
    "13-18-23.1": 0.04
};
def linked as xstr.XFrequencies init xstr.withLinkageGroup(
    $independent,
    xstr.linkageGroup("LG1", $lg1, 0.002, $haps));
def phsLinked as float init xstr.kinshipLR($alice, $beth, xstr.PATERNAL_HALF_SISTERS, $linked);
io.printf("Paternal half sisters, LG1 as a haplotype\n");
io.printf("  scored over: %v\n", xstr.pairUnits($linked, $alice, $beth));
io.printf("  LR: %f|prec=2\n", $phsLinked);

