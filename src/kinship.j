# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Relationship testing: kinship likelihood ratios and paternity.
 *
 * Two people who are related share alleles **identical by descent** (IBD) more
 * often than chance allows, and for a pair of non-inbred individuals that
 * relationship is fully described by three numbers - the IBD coefficients
 * `(k0, k1, k2)`, the probability that the pair shares 0, 1 or 2 alleles by
 * descent at a locus. Parent and child always share exactly one (0, 1, 0);
 * full siblings share (1/4, 1/2, 1/4); half siblings share (1/2, 1/2, 0).
 *
 * From those coefficients the joint genotype probability at a locus is
 *
 * ```
 *   k0 pa pb          + k1 (pa δab + pb δaa)/2 + k2 δaa δab
 * ```
 *
 * where `T1` conditions on one of B's alleles being a copy of one of A's.
 * A kinship LR is that probability under one relationship over the same under
 * another, multiplied across loci.
 *
 * **Theta is not applied here.** The joint IBD probabilities use the plain
 * product rule (theta = 0), because the sub-population correction does not
 * extend to the IBD-1 term. Use `forensics.randomMatchProbability` with a
 * theta when a single-source match probability is what you need.
 *
 * Paternity gets its own trio calculation - `paternityIndex` over mother,
 * child and alleged father - because conditioning on the mother's genotype is
 * strictly more informative than treating father and child as an anonymous
 * pair. A stepwise mutation model lets a trio survive the single locus where the
 * child carries an allele the alleged father cannot have transmitted, which would
 * otherwise be a false exclusion.
 * @module kinship
 * @example
 * import "@jennifer/forensicgenetics/kinship.j" as kinship;
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
 * `k0 + k1 + k2` must be 1. The coefficients are for non-inbred, outbred
 * pairs; an inbred pedigree needs the nine condensed Jacquard coefficients,
 * which this pair-wise form does not carry.
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
 * Parent and child, (0, 1, 0) - exactly one allele shared at every locus.
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
 * Autosomal STRs cannot tell half siblings from an uncle/nephew or a
 * grandparent/grandchild pair: all three have the same IBD coefficients, and
 * so the same likelihood ratio. Distinguishing them needs lineage markers
 * (`lineage.j`) or a denser panel.
 */
export def const HALF_SIBLINGS as Relationship init Relationship{
    name: "half siblings",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Grandparent and grandchild, (1/2, 1/2, 0) - IBD-identical to half siblings.
 */
export def const GRANDPARENT as Relationship init Relationship{
    name: "grandparent",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Uncle/aunt and nephew/niece, (1/2, 1/2, 0) - IBD-identical to half siblings.
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
 * A mutation adds or removes whole repeat units, overwhelmingly one at a time.
 * The model gives a locus a per-meiosis mutation rate and a geometric decay
 * over the number of steps: the probability of an exactly-`k`-step change in a
 * given direction is `rate * (1 - decay) * decay^(k-1) / 2`, so the rates over
 * all non-zero changes sum to `rate` and one-step changes dominate.
 *
 * Mutation rates differ substantially between loci, so `perLocus` overrides the
 * blanket `rate` where a published locus-specific rate is available.
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
 * The probability that a parent with genotype `g` transmits `allele`.
 *
 * 1.0 from a homozygote for that allele, 0.5 from a heterozygote carrying it,
 * 0.0 otherwise - Mendel's first law, before any mutation.
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
 * A stepwise mutation model.
 * @param rate {float} the per-locus, per-meiosis mutation rate, in [0, 1)
 * @param decay {float} the geometric decay per extra repeat step, in [0, 1)
 * @return {MutationModel} the model, with no per-locus overrides
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
 * @param model {MutationModel} the model
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
 * No mutation - Mendelian transmission only.
 *
 * With this model a single inconsistent locus excludes the alleged father
 * outright, which is why it is not the default.
 */
export def const NO_MUTATION as MutationModel init MutationModel{
    rate: 0.0,
    decay: 0.0,
    perLocus: {}
};

/**
 * A generic autosomal STR mutation model: 0.002 per meiosis, 0.1 decay.
 *
 * A placeholder in the right order of magnitude for paternal autosomal STR
 * mutation, with one-step changes making up 90% of mutations. Casework should
 * substitute published per-locus rates through `withLocusRate` - observed rates
 * range over roughly an order of magnitude between loci.
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
 * The probability that an allele `source` is transmitted as `target` under a
 * stepwise mutation model.
 *
 * `1 - rate` when the two are the same allele, and
 * `rate * (1 - decay) * decay^(k-1) / 2` for a change of `k` repeat steps in
 * one direction. A microvariant difference below a full repeat (15 -> 15.2)
 * counts as one step.
 * @param model {MutationModel} the model
 * @param locus {string} the locus, for a per-locus rate
 * @param source {string} the parent's allele
 * @param target {string} the transmitted allele
 * @return {float} the transmission probability under mutation
 * @throws {Error} kind "kinship" if a changed allele has no numeric repeat
 *   count, so the number of steps is undefined
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
    def geometric as float init (1.0 - $model.decay) * math.pow(
        $model.decay,
        convert.toFloat($steps - 1));
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
 * The joint probability of two genotypes at one locus under a relationship.
 *
 * `k0 P(Ga) P(Gb) + k1 T1(Ga, Gb) + k2 [Ga == Gb] P(Ga)`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the relationship being assumed
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
 * The kinship likelihood ratio at one locus, for one relationship against
 * another.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the relationship under the numerator hypothesis
 * @param alt {Relationship} the relationship under the denominator hypothesis
 * @return {float} the likelihood ratio at this locus
 * @throws {Error} kind "kinship" on mismatched loci or a zero denominator
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
 * The kinship likelihood ratio for a relationship against being unrelated.
 *
 * The product of `kinshipLRAt` over every locus the two profiles and the
 * database share. An LR above 1 favours the relationship; below 1 favours
 * unrelatedness.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship being tested
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
 * The kinship likelihood ratio for one relationship against another.
 *
 * The general form: `kinshipLR` is this with `UNRELATED` as the alternative.
 * Testing full siblings against half siblings, say, is
 * `kinshipLRAgainst(a, b, FULL_SIBLINGS, HALF_SIBLINGS, db)`.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship under the numerator hypothesis
 * @param alt {Relationship} the relationship under the denominator hypothesis
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
 * `X / Y`, where `X` is the probability of the child's genotype given the
 * mother and the alleged father, and `Y` the same given the mother and a
 * random man from the reference population. The mutation model applies to the
 * alleged father's transmission only, which is what lets a trio survive a
 * single locus where the child carries an allele he cannot have transmitted
 * unchanged.
 * @param db {FrequencyDb} the reference frequencies
 * @param mother {Genotype} the mother's genotype
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {MutationModel} the mutation model
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "kinship" if the three genotypes are not at the same
 *   locus, or if the mother cannot have transmitted either of the child's
 *   alleles (a maternal exclusion, which no paternity index can express)
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
 * The product of `paternityIndexAt` over every locus the three profiles and
 * the database share - loci only some of them carry are simply not part of the
 * comparison. Use `trioLoci` to see which loci contributed.
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
 *
 * Pass `NO_MUTATION` for the strict Mendelian calculation, in which a single
 * inconsistent locus drives the CPI to zero and excludes the alleged father.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {MutationModel} the mutation model
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
        $out = $out * paternityIndexAt(
            $db,
            $mother.loci[$locus],
            $child.loci[$locus],
            $allegedFather.loci[$locus],
            $model);
    }
    return $out;
}

/**
 * The loci at which a trio is Mendelian-inconsistent with paternity.
 *
 * The loci whose paternity index is zero without a mutation model - the ones a
 * mutation model has to explain. One such locus is the classic single
 * inconsistency; several mean exclusion.
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
 * The probability of paternity, `W = CPI / (CPI + 1)`.
 *
 * The posterior probability under the conventional neutral prior of 0.5. It is
 * a restatement of the CPI, not extra evidence: quote both, and say which
 * prior it assumes.
 * @param cpi {float} the combined paternity index
 * @return {float} the probability of paternity, in [0, 1)
 * @throws {Error} kind "kinship" on a negative CPI
 */
export func probabilityOfPaternity(cpi as float) {
    return probabilityOfPaternityPrior($cpi, 0.5);
}

/**
 * The probability of paternity under an explicit prior.
 *
 * `W = CPI * prior / (CPI * prior + (1 - prior))`. At `prior = 0.5` this is
 * `probabilityOfPaternity`.
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
