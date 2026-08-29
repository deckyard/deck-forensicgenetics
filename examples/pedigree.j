# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

# Deficient paternity: the alleged father is unavailable, so the case is made
# through his parents. This is the question `kinship.j` cannot answer - four
# people at once, joined through an untyped man - and what peeling is for.
#
#   H1                              H2
#     GF x GM                         GF   GM   (unrelated bystanders)
#        |                                 UN x MO
#       FA x MO                                |
#            |                                CH
#           CH
#
# Run from the deck root:  jennifer run examples/pedigree.j
#
# The frequency table is synthetic - see testdata/frequencies.xml.
use io;
use math;

import "../src/forensicgenetics.j" as forensics;
import "./lib/exampledb.j" as exampledb;
import "../src/pedigree.j" as pedigree;
import "../src/kinship.j" as kinship;

def db as forensics.FrequencyDb init exampledb.build();

# --- the four typed people --------------------------------------------------

def profiles as map of string to forensics.Profile init {
    "GF": forensics.profile(
        "GF",
        [
            forensics.genotype("D3S1358", "17", "18"),
            forensics.genotype("vWA", "19", "20"),
            forensics.genotype("FGA", "24", "25"),
            forensics.genotype("D8S1179", "15", "16"),
            forensics.genotype("D21S11", "31.2", "32.2"),
            forensics.genotype("D18S51", "17", "18"),
            forensics.genotype("TH01", "9.3", "10"),
            forensics.genotype("D16S539", "13", "14")
        ]),
    "GM": forensics.profile(
        "GM",
        [
            forensics.genotype("D3S1358", "12", "16"),
            forensics.genotype("vWA", "14", "15"),
            forensics.genotype("FGA", "18", "19"),
            forensics.genotype("D8S1179", "8", "9"),
            forensics.genotype("D21S11", "27", "28"),
            forensics.genotype("D18S51", "12", "13"),
            forensics.genotype("TH01", "5", "8"),
            forensics.genotype("D16S539", "8", "9")
        ]),
    "MO": forensics.profile(
        "MO",
        [
            forensics.genotype("D3S1358", "15", "16"),
            forensics.genotype("vWA", "16", "17"),
            forensics.genotype("FGA", "21", "22"),
            forensics.genotype("D8S1179", "12", "13"),
            forensics.genotype("D21S11", "29", "30"),
            forensics.genotype("D18S51", "14", "15"),
            forensics.genotype("TH01", "6", "7"),
            forensics.genotype("D16S539", "11", "12")
        ]),
    "CH": forensics.profile(
        "CH",
        [
            forensics.genotype("D3S1358", "15", "17"),
            forensics.genotype("vWA", "16", "19"),
            forensics.genotype("FGA", "21", "24"),
            forensics.genotype("D8S1179", "12", "15"),
            forensics.genotype("D21S11", "29", "31.2"),
            forensics.genotype("D18S51", "14", "17"),
            forensics.genotype("TH01", "6", "9.3"),
            forensics.genotype("D16S539", "11", "13")
        ])
};

# --- the two hypotheses -----------------------------------------------------

# H1: the child's father FA is the son of GF and GM. FA himself is untyped -
# the whole point of the exercise - so peeling sums over every genotype he
# could have had.
def h1 as pedigree.Pedigree init pedigree.pedigree("son of GF x GM");
$h1 = pedigree.withFounder($h1, "GF", "m");
$h1 = pedigree.withFounder($h1, "GM", "f");
$h1 = pedigree.withFounder($h1, "MO", "f");
$h1 = pedigree.withChild($h1, "FA", "m", "GF", "GM");
$h1 = pedigree.withChild($h1, "CH", "m", "FA", "MO");

# H2: the child's father is an unrelated man, and GF and GM are nobody's
# relatives. The same four people are typed, so the two likelihoods are
# comparable.
def h2 as pedigree.Pedigree init pedigree.pedigree("unrelated man");
$h2 = pedigree.withFounder($h2, "GF", "m");
$h2 = pedigree.withFounder($h2, "GM", "f");
$h2 = pedigree.withFounder($h2, "MO", "f");
$h2 = pedigree.withFounder($h2, "UN", "m");
$h2 = pedigree.withChild($h2, "CH", "m", "UN", "MO");

pedigree.validate($h1);
pedigree.validate($h2);
io.printf("typed in H1: %v\n", pedigree.typedMembers($h1, $profiles));
io.printf("typed in H2: %v\n\n", pedigree.typedMembers($h2, $profiles));

# --- the likelihood ratio, locus by locus -----------------------------------

io.printf("locus       L(H1)          L(H2)          LR\n");
for (def locus in pedigree.pedigreeLoci($h1, $profiles, $db)) {
    def g1 as map of string to forensics.Genotype init {};
    for (def id in $profiles) {
        if (pedigree.hasMember($h1, $id)) {
            $g1[$id] = forensics.genotypeAt($profiles[$id], $locus);
        }
    }
    def over as float init pedigree.locusLikelihood($h1, $g1, $db, $locus);
    def under as float init pedigree.locusLikelihood($h2, $g1, $db, $locus);
    io.printf("%s|pad=12 %f|prec=10 %f|prec=10 %f|prec=4\n", $locus, $over, $under, $over / $under);
}

def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
io.printf("\ncombined LR = %f|prec=2\n", $lr);
io.printf("log10 LR    = %f|prec=3\n", math.log10($lr));
io.printf("W (prior 0.5) = %f|prec=8\n", kinship.probabilityOfPaternity($lr));

# --- what peeling buys ------------------------------------------------------

# The pairwise machinery can only ask a much weaker question of the same data:
# is the child a grandchild of GF? That throws away GM entirely, and with her
# the constraint that the two grandparents must jointly explain the child's
# paternal allele.
io.printf(
    "\nfor contrast, the pairwise grandparent LR (GF vs CH alone) = %f|prec=2\n",
    kinship.kinshipLR($profiles["GF"], $profiles["CH"], kinship.GRANDPARENT, $db));
io.printf("Using both grandparents at once is worth far more than either alone.\n");
