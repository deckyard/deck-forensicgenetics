# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Pairwise kinship likelihood ratios and trio paternity.
 *
 * A non-inbred pair's relationship is its IBD coefficients `(k0, k1, k2)`: the
 * probabilities of sharing 0, 1 or 2 alleles identical by descent at a locus.
 * The joint genotype probability at one locus is
 *
 * ```
 *   k0 P(Ga) P(Gb) + k1 T1(Ga, Gb) + k2 [Ga == Gb] P(Ga)
 * ```
 *
 * where `T1` conditions on one of B's alleles being a copy of one of A's. A
 * kinship LR is that probability under one relationship over the same under
 * another, multiplied across loci.
 *
 * Theta is not applied: every probability here uses the plain product rule
 * (theta = 0), because the sub-population correction does not extend to the
 * IBD-1 term.
 *
 * `paternityIndex` conditions on the mother's genotype, which a pairwise
 * father/child LR cannot. Its stepwise mutation model keeps a single locus
 * where the child carries an allele the alleged father cannot have transmitted
 * from excluding him outright.
 * @module kinship
 * @example
 * import "@deckyard/forensicgenetics/kinship.j" as kinship;
 *
 * def pi as float init kinship.paternityIndex($mother, $child, $af, $db);
 * def w as float init kinship.probabilityOfPaternity($pi);
 *
 * def lr as float init kinship.kinshipLR($a, $b, kinship.FULL_SIBLINGS, $db);
 */
use math;
use maps;
use convert;

import "./forensicgenetics.j" as forensics;

/**
 * A relationship, as its IBD coefficients.
 *
 * `k0 + k1 + k2` is 1. The form covers non-inbred pairs only; an inbred
 * pedigree needs the nine condensed Jacquard coefficients.
 * @field name {string} a label for the relationship
 * @field k0 {float} the probability of sharing no allele identical by descent
 * @field k1 {float} the probability of sharing exactly one
 * @field k2 {float} the probability of sharing both
 */
export def struct Relationship {
    name as string,
    k0 as float,
    k1 as float,
    k2 as float
};

/**
 * Unrelated: the default alternative hypothesis, (1, 0, 0).
 */
export def const UNRELATED as Relationship init Relationship{
    name: "unrelated",
    k0: 1.0,
    k1: 0.0,
    k2: 0.0
};

/**
 * Parent and child, (0, 1, 0): exactly one allele shared at every locus.
 */
export def const PARENT_CHILD as Relationship init Relationship{
    name: "parent-child",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Full siblings, (1/4, 1/2, 1/4).
 */
export def const FULL_SIBLINGS as Relationship init Relationship{
    name: "full siblings",
    k0: 0.25,
    k1: 0.5,
    k2: 0.25
};

/**
 * Half siblings, (1/2, 1/2, 0).
 *
 * `GRANDPARENT` and `AVUNCULAR` have the same coefficients and so the same
 * likelihood ratio; autosomal STRs cannot separate the three.
 */
export def const HALF_SIBLINGS as Relationship init Relationship{
    name: "half siblings",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Grandparent and grandchild, (1/2, 1/2, 0); IBD-identical to half siblings.
 */
export def const GRANDPARENT as Relationship init Relationship{
    name: "grandparent",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Uncle/aunt and nephew/niece, (1/2, 1/2, 0); IBD-identical to half siblings.
 */
export def const AVUNCULAR as Relationship init Relationship{
    name: "avuncular",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * First cousins, (3/4, 1/4, 0).
 */
export def const FIRST_COUSINS as Relationship init Relationship{
    name: "first cousins",
    k0: 0.75,
    k1: 0.25,
    k2: 0.0
};

/**
 * The same person, or monozygotic twins, (0, 0, 1).
 */
export def const IDENTICAL as Relationship init Relationship{
    name: "identical",
    k0: 0.0,
    k1: 0.0,
    k2: 1.0
};

/**
 * A stepwise STR mutation model.
 *
 * The probability of an exactly-`k`-step change in one direction is
 * `rate * (1 - decay) * decay^(k-1) / 2`; over all non-zero changes in both
 * directions these sum to `rate`. `perLocus` overrides `rate` for named loci.
 * @field rate {float} the per-locus, per-meiosis mutation rate
 * @field decay {float} the geometric decay per additional repeat step, in [0, 1)
 * @field perLocus {map of string to float} locus -> rate, overriding `rate`
 */
export def struct MutationModel {
    rate as float,
    decay as float,
    perLocus as map of string to float
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "kinship", message: $msg, file: "", line: 0, col: 0};
}

# --- relationships ----------------------------------------------------------

/**
 * A relationship from its IBD coefficients.
 * @param name {string} a label for the relationship
 * @param k0 {float} the probability of sharing no allele IBD
 * @param k1 {float} the probability of sharing exactly one
 * @param k2 {float} the probability of sharing both
 * @return {Relationship} the relationship
 * @throws {Error} kind "kinship" if a coefficient is outside [0, 1] or the
 *   three do not sum to 1
 */
export func relationship(name as string, k0 as float, k1 as float, k2 as float) {
    for (def k in [$k0, $k1, $k2]) {
        if ($k < 0.0 or $k > 1.0) {
            fail("relationship '" + $name + "': IBD coefficients must be in [0, 1]");
        }
    }
    if (math.abs($k0 + $k1 + $k2 - 1.0) > 0.000000001) {
        fail("relationship '" + $name + "': IBD coefficients must sum to 1");
    }
    return Relationship{name: $name, k0: $k0, k1: $k1, k2: $k2};
}

# --- transmission and mutation ----------------------------------------------

/**
 * The Mendelian probability that a parent with genotype `g` transmits
 * `allele`: 1.0 from a homozygote for it, 0.5 from a heterozygote carrying it,
 * 0.0 otherwise. No mutation.
 * @param g {Genotype} the parent's genotype
 * @param allele {string} the allele label
 * @return {float} the transmission probability
 */
export func transmissionProbability(g as forensics.Genotype, allele as string) {
    def n as float init 0.0;
    if ($g.a == $allele) {
        $n = $n + 1.0;
    }
    if ($g.b == $allele) {
        $n = $n + 1.0;
    }
    return $n / 2.0;
}

/**
 * A stepwise mutation model with no per-locus overrides.
 * @param rate {float} the per-locus, per-meiosis mutation rate, in [0, 1)
 * @param decay {float} the geometric decay per extra repeat step, in [0, 1)
 * @return {MutationModel} the model
 * @throws {Error} kind "kinship" on an out-of-range rate or decay
 */
export func mutationModel(rate as float, decay as float) {
    if ($rate < 0.0 or $rate >= 1.0) {
        fail("mutation rate must be in [0, 1)");
    }
    if ($decay < 0.0 or $decay >= 1.0) {
        fail("mutation decay must be in [0, 1)");
    }
    return MutationModel{rate: $rate, decay: $decay, perLocus: {}};
}

/**
 * A copy of `model` with one locus's mutation rate overridden.
 * @param model {MutationModel} the model to copy
 * @param locus {string} the locus name
 * @param rate {float} that locus's per-meiosis mutation rate, in [0, 1)
 * @return {MutationModel} the updated copy
 * @throws {Error} kind "kinship" on an out-of-range rate
 */
export func withLocusRate(model as MutationModel, locus as string, rate as float) {
    if ($rate < 0.0 or $rate >= 1.0) {
        fail("mutation rate must be in [0, 1)");
    }
    def out as MutationModel init $model;
    $out.perLocus[$locus] = $rate;
    return $out;
}

/**
 * No mutation: Mendelian transmission only. A single inconsistent locus
 * excludes the alleged father.
 */
export def const NO_MUTATION as MutationModel init MutationModel{
    rate: 0.0,
    decay: 0.0,
    perLocus: {}
};

/**
 * A generic autosomal STR mutation model: 0.002 per meiosis, 0.1 decay.
 *
 * A placeholder of the right order of magnitude, not a published rate;
 * one-step changes are 90% of mutations. Observed rates vary about tenfold
 * between loci, so set published ones with `withLocusRate`.
 */
export def const DEFAULT_MUTATION as MutationModel init MutationModel{
    rate: 0.002,
    decay: 0.1,
    perLocus: {}
};

# rateFor returns the mutation rate in force at a locus.
func rateFor(model as MutationModel, locus as string) {
    if (maps.has($model.perLocus, $locus)) {
        return $model.perLocus[$locus];
    }
    return $model.rate;
}

/**
 * The probability that allele `source` is transmitted as `target`.
 *
 * `1 - rate` when they are equal, and `rate * (1 - decay) * decay^(k-1) / 2`
 * for a change of `k` repeat steps, `k` being the absolute difference in
 * repeat count rounded to the nearest integer and at least 1 (15 -> 15.2 is
 * one step). With a rate of 0 a change has probability 0.
 * @param model {MutationModel} the mutation model
 * @param locus {string} the locus, for a per-locus rate
 * @param source {string} the parent's allele
 * @param target {string} the transmitted allele
 * @return {float} the transmission probability under mutation
 * @throws {Error} kind "kinship" if the rate is above 0 and either allele of a
 *   change has no numeric repeat count
 */
export func mutationProbability(
    model as MutationModel,
    locus as string,
    source as string,
    target as string) {
    def rate as float init rateFor($model, $locus);
    if ($source == $target) {
        return 1.0 - $rate;
    }
    if ($rate <= 0.0) {
        return 0.0;
    }
    if (not forensics.isNumericAllele($source) or not forensics.isNumericAllele($target)) {
        fail("locus '" + $locus + "': cannot count mutation steps between '" + $source + "' and '" +
            $target + "' - neither is a repeat count");
    }
    def va as float init forensics.alleleValue($source);
    def vb as float init forensics.alleleValue($target);
    def delta as float init math.abs($va - $vb);
    def steps as int init math.round($delta);
    if ($steps < 1) {
        $steps = 1;
    }
    def geometric as float init (1.0 - $model.decay) *
        math.pow($model.decay, convert.toFloat($steps - 1));
    return $rate * $geometric / 2.0;
}

# paternalTransmission is `transmissionProbability` with mutation folded in:
# the chance that a parent with genotype `g` passes on `allele`, allowing for a
# stepwise change. With NO_MUTATION it reduces exactly to the Mendelian value.
func paternalTransmission(
    model as MutationModel,
    locus as string,
    g as forensics.Genotype,
    allele as string) {
    return 0.5 * mutationProbability($model, $locus, $g.a, $allele) +
        0.5 * mutationProbability($model, $locus, $g.b, $allele);
}

# --- joint genotype probabilities -------------------------------------------

# genotypeProbability is the plain Hardy-Weinberg genotype frequency (theta 0).
func genotypeProbability(db as forensics.FrequencyDb, g as forensics.Genotype) {
    return forensics.genotypeFrequency($db, $g, 0.0);
}

# conditionalOnAllele is P(genotype `g` | one of its alleles is fixed at
# `allele`, the other drawn from the population). Zero when `g` does not carry
# `allele` at all.
func conditionalOnAllele(db as forensics.FrequencyDb, g as forensics.Genotype, allele as string) {
    if ($g.a == $g.b) {
        if ($g.a != $allele) {
            return 0.0;
        }
        return forensics.alleleFrequency($db, $g.locus, $g.a);
    }
    if ($g.a == $allele) {
        return forensics.alleleFrequency($db, $g.locus, $g.b);
    }
    if ($g.b == $allele) {
        return forensics.alleleFrequency($db, $g.locus, $g.a);
    }
    return 0.0;
}

# ibdOneTerm is T1: the joint probability that A has genotype `ga` and B has
# `gb` given they share exactly one allele identical by descent. The IBD allele
# is one of A's two, each with probability 1/2; B's other allele comes from the
# population.
func ibdOneTerm(db as forensics.FrequencyDb, ga as forensics.Genotype, gb as forensics.Genotype) {
    def shared as float init 0.5 * conditionalOnAllele($db, $gb, $ga.a) +
        0.5 * conditionalOnAllele($db, $gb, $ga.b);
    return genotypeProbability($db, $ga) * $shared;
}

/**
 * The joint probability of two genotypes at one locus under a relationship:
 * `k0 P(Ga) P(Gb) + k1 T1(Ga, Gb) + k2 [Ga == Gb] P(Ga)`, theta 0.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the assumed relationship
 * @return {float} the joint probability
 * @throws {Error} kind "kinship" if the two genotypes are at different loci,
 *   or kind "forensicgenetics" if the database does not cover the locus
 */
export func jointGenotypeProbability(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship) {
    if ($ga.locus != $gb.locus) {
        fail("cannot pair a genotype at '" + $ga.locus + "' with one at '" + $gb.locus + "'");
    }
    def both as float init 0.0;
    if ($ga.a == $gb.a and $ga.b == $gb.b) {
        $both = genotypeProbability($db, $ga);
    }
    return $rel.k0 * genotypeProbability($db, $ga) * genotypeProbability($db, $gb) +
        $rel.k1 * ibdOneTerm($db, $ga, $gb) + $rel.k2 * $both;
}

/**
 * The kinship likelihood ratio at one locus, `rel` against `alt`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @return {float} the likelihood ratio at this locus
 * @throws {Error} kind "kinship" on mismatched loci, or if `alt` gives the
 *   pair probability zero
 */
export func kinshipLRAt(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship,
    alt as Relationship) {
    def under as float init jointGenotypeProbability($db, $ga, $gb, $alt);
    if ($under <= 0.0) {
        fail("locus '" + $ga.locus + "': the alternative hypothesis '" + $alt.name +
            "' gives these genotypes probability zero, so no ratio is defined");
    }
    return jointGenotypeProbability($db, $ga, $gb, $rel) / $under;
}

# sharedLoci returns the loci two profiles and the database all cover.
func sharedLoci(db as forensics.FrequencyDb, a as forensics.Profile, b as forensics.Profile) {
    def out as list of string init [];
    for (def locus in $a.loci) {
        if (maps.has($b.loci, $locus) and forensics.hasLocus($db, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * The loci a pair of profiles and the database all cover, in the first
 * profile's order.
 * @param db {FrequencyDb} the reference frequencies
 * @param a {Profile} one profile
 * @param b {Profile} the other
 * @return {list of string} the shared locus names
 */
export func pairLoci(db as forensics.FrequencyDb, a as forensics.Profile, b as forensics.Profile) {
    return sharedLoci($db, $a, $b);
}

/**
 * The kinship likelihood ratio for a relationship against `UNRELATED`.
 *
 * The product of `kinshipLRAt` over `pairLoci`; loci outside that set are
 * ignored. Above 1 favours the relationship.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship to test
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood ratio
 * @throws {Error} kind "kinship" if the profiles and database share no locus
 */
export func kinshipLR(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    db as forensics.FrequencyDb) {
    return kinshipLRAgainst($a, $b, $rel, UNRELATED, $db);
}

/**
 * The kinship likelihood ratio for one relationship against another, over
 * `pairLoci`. `kinshipLR` is this with `UNRELATED` as `alt`.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood ratio
 * @throws {Error} kind "kinship" if the profiles and database share no locus,
 *   or if the alternative gives an observed pair probability zero
 */
export func kinshipLRAgainst(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    alt as Relationship,
    db as forensics.FrequencyDb) {
    def loci as list of string init sharedLoci($db, $a, $b);
    if (len($loci) == 0) {
        fail("profiles '" + $a.id + "' and '" + $b.id +
            "' share no locus the frequency database covers");
    }
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = $out * kinshipLRAt($db, $a.loci[$locus], $b.loci[$locus], $rel, $alt);
    }
    return $out;
}

# --- paternity --------------------------------------------------------------

# childProbability is P(child genotype | mother's genotype, father's genotype),
# summed over which of the child's alleles came from which parent. `fatherOf`
# is the father-side transmission function, so the same code serves the alleged
# father (with mutation) and a random man (population frequencies).
func childProbability(
    mg as forensics.Genotype,
    cg as forensics.Genotype,
    fromFather as map of string to float) {
    if ($cg.a == $cg.b) {
        return transmissionProbability($mg, $cg.a) * $fromFather[$cg.a];
    }
    return transmissionProbability($mg, $cg.a) * $fromFather[$cg.b] +
        transmissionProbability($mg, $cg.b) * $fromFather[$cg.a];
}

/**
 * The paternity index at one locus.
 *
 * `X / Y`: the probability of the child's genotype given the mother and the
 * alleged father, over the same given the mother and a random man from the
 * reference population. The mutation model applies to the alleged father's
 * transmission only.
 * @param db {FrequencyDb} the reference frequencies
 * @param mother {Genotype} the mother's genotype
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {MutationModel} the alleged father's mutation model
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "kinship" if the three genotypes are not at one locus,
 *   or on a maternal exclusion (the mother can transmit neither of the child's
 *   alleles)
 */
export func paternityIndexAt(
    db as forensics.FrequencyDb,
    mother as forensics.Genotype,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as MutationModel) {
    if ($mother.locus != $child.locus or $child.locus != $father.locus) {
        fail("the mother, child and alleged father genotypes must be at one locus");
    }
    def locus as string init $child.locus;
    def alleged as map of string to float init {};
    def random as map of string to float init {};
    for (def a in [$child.a, $child.b]) {
        $alleged[$a] = paternalTransmission($model, $locus, $father, $a);
        $random[$a] = forensics.alleleFrequency($db, $locus, $a);
    }
    def y as float init childProbability($mother, $child, $random);
    if ($y <= 0.0) {
        fail("locus '" + $locus + "': the mother cannot have transmitted either of the child's" +
            " alleles - this is a maternal exclusion, not a paternity index");
    }
    return childProbability($mother, $child, $alleged) / $y;
}

/**
 * The loci a mother/child/alleged-father trio and the database all cover, in
 * the child's order.
 * @param db {FrequencyDb} the reference frequencies
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param father {Profile} the alleged father's profile
 * @return {list of string} the shared locus names
 */
export func trioLoci(
    db as forensics.FrequencyDb,
    mother as forensics.Profile,
    child as forensics.Profile,
    father as forensics.Profile) {
    def out as list of string init [];
    for (def locus in $child.loci) {
        if (maps.has($mother.loci, $locus) and maps.has($father.loci, $locus) and
            forensics.hasLocus($db, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * The combined paternity index (CPI) for a trio, under `DEFAULT_MUTATION`.
 *
 * The product of `paternityIndexAt` over `trioLoci`; loci outside that set
 * are ignored.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" if the trio and database share no locus, or on
 *   a maternal exclusion
 */
export func paternityIndex(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb) {
    return paternityIndexWith($mother, $child, $allegedFather, $db, DEFAULT_MUTATION);
}

/**
 * The combined paternity index for a trio under an explicit mutation model.
 * With `NO_MUTATION` one inconsistent locus makes the CPI zero.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {MutationModel} the alleged father's mutation model
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" if the trio and database share no locus, or on
 *   a maternal exclusion
 */
export func paternityIndexWith(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as MutationModel) {
    def loci as list of string init trioLoci($db, $mother, $child, $allegedFather);
    if (len($loci) == 0) {
        fail("the trio and the frequency database share no locus");
    }
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = $out *
            paternityIndexAt(
            $db,
            $mother.loci[$locus],
            $child.loci[$locus],
            $allegedFather.loci[$locus],
            $model);
    }
    return $out;
}

/**
 * The loci of `trioLoci` whose paternity index is zero under `NO_MUTATION`:
 * the inconsistencies a mutation model has to explain.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {list of string} the inconsistent locus names
 * @throws {Error} kind "kinship" on a maternal exclusion
 */
export func inconsistentLoci(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb) {
    def out as list of string init [];
    for (def locus in trioLoci($db, $mother, $child, $allegedFather)) {
        def pi as float init paternityIndexAt(
            $db,
            $mother.loci[$locus],
            $child.loci[$locus],
            $allegedFather.loci[$locus],
            NO_MUTATION);
        if ($pi <= 0.0) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * The probability of paternity under a prior of 0.5: `W = CPI / (CPI + 1)`.
 *
 * A restatement of the CPI, not further evidence; report the prior with it.
 * @param cpi {float} the combined paternity index
 * @return {float} the probability of paternity, in [0, 1)
 * @throws {Error} kind "kinship" on a negative CPI
 */
export func probabilityOfPaternity(cpi as float) {
    return probabilityOfPaternityPrior($cpi, 0.5);
}

/**
 * The probability of paternity under an explicit prior:
 * `W = CPI * prior / (CPI * prior + (1 - prior))`.
 * @param cpi {float} the combined paternity index
 * @param prior {float} the prior probability of paternity, in (0, 1)
 * @return {float} the posterior probability of paternity
 * @throws {Error} kind "kinship" on a negative CPI or a prior outside (0, 1)
 */
export func probabilityOfPaternityPrior(cpi as float, prior as float) {
    if ($cpi < 0.0) {
        fail("a combined paternity index cannot be negative");
    }
    if ($prior <= 0.0 or $prior >= 1.0) {
        fail("the prior probability of paternity must be in (0, 1)");
    }
    def weighted as float init $cpi * $prior;
    return $weighted / ($weighted + (1.0 - $prior));
}
