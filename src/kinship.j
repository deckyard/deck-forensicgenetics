# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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
 * The original functions use the product rule (theta = 0). Their `...Theta`
 * variants draw every allele not identical by descent by the Balding-Nichols
 * sampling formula, conditioned on the alleles already drawn: four draws
 * under IBD 0, three under IBD 1, two under IBD 2. Putting the NRC II
 * genotype frequency into the theta-0 form instead would leave the IBD-1 term
 * on the product rule. At theta 0 every `...Theta` variant returns exactly what the
 * original returns. `duoPaternityIndex` covers an alleged father and child
 * without the mother.
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
 * def lr01 as float init kinship.kinshipLRTheta($a, $b, kinship.FULL_SIBLINGS, $db, 0.01);
 */
use math;
use maps;
use convert;
use decimal;

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

# FLOAT_MAX is the largest finite float. The interpreter raises a runtime error
# on a product past it, so the across-loci products check first.
def const FLOAT_MAX as float init 1.7976931348623157e308;

# timesLR multiplies one locus's factor into a running LR, throwing this
# module's error where the product would leave the float range. `acc <= 1`
# cannot overflow, and for `acc > 1` the quotient `FLOAT_MAX / acc` is finite.
# The message is assembled only on failure.
func timesLR(
    acc as float,
    factor as float,
    what as string,
    idA as string,
    idB as string,
    hint as string) {
    if ($acc > 1.0 and $factor > FLOAT_MAX / $acc) {
        fail("the " + $what + " of '" + $idA + "' and '" + $idB +
            "' exceeds the float range (above 1.8e308)" + $hint);
    }
    return $acc * $factor;
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
 * between loci, so set published ones with `withLocusRate`, or all of the
 * AABB 2003 rates bundled in `mutation.j` with `mutation.kinshipModelWithRates`.
 */
export def const DEFAULT_MUTATION as MutationModel init MutationModel{
    rate: 0.002,
    decay: 0.1,
    perLocus: {}
};

# rateFor returns the mutation rate in force at a locus, checking it and the
# decay: a `MutationModel` literal skips `mutationModel`'s checks, and a rate
# past 1 made a transmission probability negative.
func rateFor(model as MutationModel, locus as string) {
    def rate as float init $model.rate;
    if (maps.has($model.perLocus, $locus)) {
        $rate = $model.perLocus[$locus];
    }
    if ($rate < 0.0 or $rate >= 1.0) {
        fail("locus '" + $locus + "': mutation rate {$rate} must be in [0, 1)");
    }
    if ($model.decay < 0.0 or $model.decay >= 1.0) {
        fail("mutation decay {$model.decay} must be in [0, 1)");
    }
    return $rate;
}

/**
 * The probability that allele `source` is transmitted as `target`.
 *
 * `1 - rate` when they are equal, and `rate * (1 - decay) * decay^(k-1) / 2`
 * for a change of `k` repeat steps, `k` being the absolute difference in
 * repeat count rounded to the nearest integer (halves up) and at least 1
 * (15 -> 15.2 is one step). The difference is taken in exact decimal
 * arithmetic on the labels, so `10.15 -> 10.65` is exactly half a step and
 * rounds to 1; a difference too large for `decay^(k-1)` to be non-zero gives
 * 0. With a rate of 0 a change has probability 0.
 * @param model {MutationModel} the mutation model
 * @param locus {string} the locus, for a per-locus rate
 * @param source {string} the parent's allele
 * @param target {string} the transmitted allele
 * @return {float} the transmission probability under mutation
 * @throws {Error} kind "kinship" if the rate in force is outside [0, 1), the
 *   decay is outside [0, 1), or the rate is above 0 and either allele of a
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
            $target + "' - " + nonNumeric($source, $target) + " a repeat count");
    }
    # Exact decimal difference of the labels, rounded half up. A float
    # difference misplaces halves (10.65 - 10.15 is 0.5000000000000009), and
    # math.round throws once the count leaves the int range; as a float, a
    # count past 2^53 is still exact enough that decay^(k-1) is 0.
    def delta as decimal.Decimal init decimal.abs(decimal.sub($source, $target));
    def steps as float init decimal.toFloat(decimal.round($delta, 0, "half-up"));
    if ($steps < 1.0) {
        $steps = 1.0;
    }
    def geometric as float init (1.0 - $model.decay) * math.pow($model.decay, $steps - 1.0);
    return $rate * $geometric / 2.0;
}

# nonNumeric names the labels of a pair that are not repeat counts.
func nonNumeric(source as string, target as string) {
    if (forensics.isNumericAllele($source)) {
        return "'" + $target + "' is not";
    }
    if (forensics.isNumericAllele($target)) {
        return "'" + $source + "' is not";
    }
    return "neither is";
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

# Freqs is the floored frequencies of the alleles one locus's calculation
# reads, resolved once by `freqsOf`. A call into forensicgenetics.j rebuilds
# every argument holding one of its structs, so `forensics.alleleFrequency`
# costs the whole database per call (6.5 ms per 24-locus `kinshipLR`).
def struct Freqs {
    p as map of string to float
};

# freqsOf resolves the floored frequencies of `alleles` at `locus` without a
# call into the core: `forensics.alleleFrequency`'s floor and cap, term for
# term. An uncovered locus, or an unlisted allele with no floor, is passed to
# `forensics.alleleFrequency` to throw the core's own error.
func freqsOf(db as forensics.FrequencyDb, locus as string, alleles as list of string) {
    if (not maps.has($db.loci, $locus)) {
        forensics.alleleFrequency($db, $locus, $alleles[0]);
    }
    def n as int init $db.loci[$locus].n;
    def floor as float init 0.0;
    if ($n > 0) {
        $floor = forensics.MIN_ALLELE_COUNT / (2.0 * convert.toFloat($n));
    }
    def out as Freqs init Freqs{p: {}};
    for (def a in $alleles) {
        def p as float init 0.0;
        if (maps.has($db.loci[$locus].alleles, $a)) {
            $p = $db.loci[$locus].alleles[$a];
        }
        if ($p < $floor) {
            $p = $floor;
        }
        if ($p > 1.0) {
            $p = 1.0;
        }
        if ($p <= 0.0) {
            forensics.alleleFrequency($db, $locus, $a);
        }
        $out.p[$a] = $p;
    }
    return $out;
}

# pairFreqs is `freqsOf` over two genotypes' alleles.
func pairFreqs(db as forensics.FrequencyDb, ga as forensics.Genotype, gb as forensics.Genotype) {
    return freqsOf($db, $ga.locus, [$ga.a, $ga.b, $gb.a, $gb.b]);
}

# genotypeProbability is the plain Hardy-Weinberg genotype frequency (theta 0),
# capped at 1: `forensics.genotypeFrequency` at theta 0, term for term.
func genotypeProbability(f as Freqs, g as forensics.Genotype) {
    def pa as float init $f.p[$g.a];
    def out as float init 0.0;
    if ($g.a == $g.b) {
        $out = $pa * $pa;
    } else {
        $out = 2.0 * $pa * $f.p[$g.b];
    }
    if ($out > 1.0) {
        return 1.0;
    }
    return $out;
}

# conditionalOnAllele is P(genotype `g` | one of its alleles is fixed at
# `allele`, the other drawn from the population). Zero when `g` does not carry
# `allele` at all.
func conditionalOnAllele(f as Freqs, g as forensics.Genotype, allele as string) {
    if ($g.a == $g.b) {
        if ($g.a != $allele) {
            return 0.0;
        }
        return $f.p[$g.a];
    }
    if ($g.a == $allele) {
        return $f.p[$g.b];
    }
    if ($g.b == $allele) {
        return $f.p[$g.a];
    }
    return 0.0;
}

# Terms is one locus's joint-probability terms before the IBD weights, so a
# ratio computes them once for both hypotheses. Theta 0: `pa`, `pb` are
# P(Ga), P(Gb) (the IBD-0 term is their product); with theta, `none` is the
# IBD-0 term. `one` is T1, `both` the IBD-2 term.
def struct Terms {
    pa as float,
    pb as float,
    none as float,
    one as float,
    both as float
};

# productTerms is `Terms` at theta 0. T1: the IBD allele is one of A's two,
# each with probability 1/2; B's other allele comes from the population.
func productTerms(f as Freqs, ga as forensics.Genotype, gb as forensics.Genotype) {
    def pa as float init genotypeProbability($f, $ga);
    def shared as float init 0.5 * conditionalOnAllele($f, $gb, $ga.a) +
        0.5 * conditionalOnAllele($f, $gb, $ga.b);
    def both as float init 0.0;
    if (sameGenotype($ga, $gb)) {
        $both = $pa;
    }
    return Terms{
        pa: $pa,
        pb: genotypeProbability($f, $gb),
        none: 0.0,
        one: $pa * $shared,
        both: $both
    };
}

# productJoint weights theta-0 terms by a relationship.
func productJoint(t as Terms, rel as Relationship) {
    return $rel.k0 * $t.pa * $t.pb + $rel.k1 * $t.one + $rel.k2 * $t.both;
}

# checkRelationship repeats `relationship`'s checks: a `Relationship` literal
# skips them, and negative coefficients gave negative probabilities.
func checkRelationship(rel as Relationship) {
    if ($rel.k0 < 0.0 or $rel.k0 > 1.0 or $rel.k1 < 0.0 or $rel.k1 > 1.0 or $rel.k2 < 0.0 or
        $rel.k2 > 1.0) {
        fail("relationship '" + $rel.name + "': IBD coefficients must be in [0, 1]");
    }
    if (math.abs($rel.k0 + $rel.k1 + $rel.k2 - 1.0) > 0.000000001) {
        fail("relationship '" + $rel.name + "': IBD coefficients must sum to 1");
    }
}

# checkPair rejects genotypes at different loci.
func checkPair(ga as forensics.Genotype, gb as forensics.Genotype) {
    if ($ga.locus != $gb.locus) {
        fail("cannot pair a genotype at '" + $ga.locus + "' with one at '" + $gb.locus + "'");
    }
}

/**
 * The joint probability of two genotypes at one locus under a relationship:
 * `k0 P(Ga) P(Gb) + k1 T1(Ga, Gb) + k2 [Ga == Gb] P(Ga)`, theta 0.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the assumed relationship
 * @return {float} the joint probability
 * @throws {Error} kind "kinship" if the two genotypes are at different loci or
 *   `rel` fails `relationship`'s checks, or kind "forensicgenetics" if the
 *   database does not cover the locus
 */
export func jointGenotypeProbability(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship) {
    checkPair($ga, $gb);
    checkRelationship($rel);
    return productJoint(productTerms(pairFreqs($db, $ga, $gb), $ga, $gb), $rel);
}

/**
 * The kinship likelihood ratio at one locus, `rel` against `alt`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @return {float} the likelihood ratio at this locus
 * @throws {Error} kind "kinship" on mismatched loci, if `rel` or `alt` fails
 *   `relationship`'s checks, or if `alt` gives the pair probability zero
 */
export func kinshipLRAt(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship,
    alt as Relationship) {
    checkPair($ga, $gb);
    checkRelationship($alt);
    checkRelationship($rel);
    def t as Terms init productTerms(pairFreqs($db, $ga, $gb), $ga, $gb);
    def under as float init productJoint($t, $alt);
    if ($under <= 0.0) {
        failZero($ga.locus, $alt);
    }
    return productJoint($t, $rel) / $under;
}

func failZero(locus as string, alt as Relationship) {
    fail("locus '" + $locus + "': the alternative hypothesis '" + $alt.name +
        "' gives these genotypes probability zero, so no ratio is defined");
}

# sameGenotype compares two genotypes as unordered allele pairs. `genotype`
# canonicalises the order, but a Genotype literal need not, and the k2 term
# must not miss a match over the order its alleles were written in.
func sameGenotype(x as forensics.Genotype, y as forensics.Genotype) {
    return ($x.a == $y.a and $x.b == $y.b) or ($x.a == $y.b and $x.b == $y.a);
}

# commonLoci returns the loci both profiles carry, whether or not the database
# covers them.
func commonLoci(a as forensics.Profile, b as forensics.Profile) {
    def out as list of string init [];
    for (def locus in $a.loci) {
        if (maps.has($b.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

# sharedLoci returns the loci two profiles and the database all cover.
func sharedLoci(db as forensics.FrequencyDb, a as forensics.Profile, b as forensics.Profile) {
    def out as list of string init [];
    for (def locus in commonLoci($a, $b)) {
        if (maps.has($db.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

# requireCovered throws when the people are all typed at a locus the database
# does not cover. Dropping that locus would hide whatever it says - most
# dangerously an exclusion - behind a ratio built from the others, with no
# sign in the output. The caller narrows the profiles with forensics.restrict
# if the omission is deliberate, so it is a decision on record.
func requireCovered(db as forensics.FrequencyDb, loci as list of string, who as string) {
    for (def locus in $loci) {
        if (not maps.has($db.loci, $locus)) {
            fail($who + " are all typed at locus '" + $locus + "', which frequency database '" +
                $db.name + "' does not cover; it will not be dropped silently - narrow the" +
                " profiles with forensics.restrict if that is intended");
        }
    }
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
 * The product of `kinshipLRAt` over `pairLoci`. A locus only one profile
 * carries is skipped; a locus both carry but the database does not cover
 * throws, since dropping it could hide an exclusion. Above 1 favours the
 * relationship.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship to test
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood ratio
 * @throws {Error} kind "kinship" if the profiles share no locus, share one
 *   the database does not cover, or if the LR exceeds the float range (use
 *   `log10KinshipLR`)
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
 * @throws {Error} kind "kinship" if the profiles share no locus, share one the
 *   database does not cover, if the alternative gives an observed pair
 *   probability zero, or if the LR exceeds the float range (use
 *   `log10KinshipLRAgainst`)
 */
export func kinshipLRAgainst(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    alt as Relationship,
    db as forensics.FrequencyDb) {
    def loci as list of string init commonLoci($a, $b);
    if (len($loci) == 0) {
        fail("profiles '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    requireCovered($db, $loci, "profiles '" + $a.id + "' and '" + $b.id + "'");
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            kinshipLRAt($db, $a.loci[$locus], $b.loci[$locus], $rel, $alt),
            "kinship LR",
            $a.id,
            $b.id,
            "; use log10KinshipLR or log10KinshipLRAgainst");
    }
    return $out;
}

/**
 * The base-10 logarithm of `kinshipLR`, summed per locus so it neither
 * overflows nor underflows where the product would.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship to test
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "kinship" as for `kinshipLR`, or if a locus has ratio 0
 *   (an exclusion, which has no logarithm)
 */
export func log10KinshipLR(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    db as forensics.FrequencyDb) {
    return log10KinshipLRAgainstTheta($a, $b, $rel, UNRELATED, $db, 0.0);
}

/**
 * The base-10 logarithm of `kinshipLRAgainst`, summed per locus.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "kinship" as for `kinshipLRAgainst`, or if a locus has
 *   ratio 0 (an exclusion, which has no logarithm)
 */
export func log10KinshipLRAgainst(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    alt as Relationship,
    db as forensics.FrequencyDb) {
    return log10KinshipLRAgainstTheta($a, $b, $rel, $alt, $db, 0.0);
}

# --- theta: the Balding-Nichols sampling formula ----------------------------

# checkTheta rejects a theta outside [0, 1), with this module's error kind. At
# theta = 1 the first draw's denominator (1 - theta) is zero.
func checkTheta(theta as float) {
    if ($theta < 0.0 or $theta >= 1.0) {
        fail("theta must be in [0, 1)");
    }
}

# drawProbability is the Balding-Nichols (1994) sampling formula: the
# probability that the next allele drawn from the sub-population is `allele`,
# given the alleles already drawn (`sampled`), `m` of which are `allele`:
# (m theta + (1 - theta) p) / (1 + (n - 1) theta), n = len(sampled).
func drawProbability(f as Freqs, sampled as list of string, allele as string, theta as float) {
    def m as int init 0;
    for (def s in $sampled) {
        if ($s == $allele) {
            $m = $m + 1;
        }
    }
    def n as float init convert.toFloat(len($sampled));
    def p as float init $f.p[$allele];
    return (convert.toFloat($m) * $theta + (1.0 - $theta) * $p) / (1.0 + ($n - 1.0) * $theta);
}

# sequenceProbability is the probability of drawing `alleles` from the
# sub-population in exactly that order: the product of drawProbability over
# each prefix.
func sequenceProbability(f as Freqs, alleles as list of string, theta as float) {
    def out as float init 1.0;
    def seen as list of string init [];
    for (def a in $alleles) {
        $out = $out * drawProbability($f, $seen, $a, $theta);
        $seen[] = $a;
    }
    return $out;
}

# orderings is the number of ordered allele pairs an unordered genotype stands
# for: 1 for a homozygote, 2 for a heterozygote.
func orderings(g as forensics.Genotype) {
    if ($g.a == $g.b) {
        return 1.0;
    }
    return 2.0;
}

# conditionalOnAlleleTheta is conditionalOnAllele with the free allele drawn by
# the sampling formula given `sampled` (the other person's two alleles).
func conditionalOnAlleleTheta(
    f as Freqs,
    g as forensics.Genotype,
    allele as string,
    sampled as list of string,
    theta as float) {
    if ($g.a == $g.b) {
        if ($g.a != $allele) {
            return 0.0;
        }
        return drawProbability($f, $sampled, $g.a, $theta);
    }
    if ($g.a == $allele) {
        return drawProbability($f, $sampled, $g.b, $theta);
    }
    if ($g.b == $allele) {
        return drawProbability($f, $sampled, $g.a, $theta);
    }
    return 0.0;
}

# thetaTerms is `Terms` with every allele that is not identical by descent
# drawn from one sub-population by the sampling formula: IBD 0 draws four
# alleles, IBD 1 three (B's IBD allele is a copy, not a draw), IBD 2 two. At
# theta 0 it is the product-rule form term by term; the exported functions
# still route theta 0 through that form so it stays byte-identical.
func thetaTerms(f as Freqs, ga as forensics.Genotype, gb as forensics.Genotype, theta as float) {
    def sa as list of string init [$ga.a, $ga.b];
    def pa as float init orderings($ga) * sequenceProbability($f, $sa, $theta);
    def none as float init orderings($ga) *
        orderings($gb) *
        sequenceProbability($f, [$ga.a, $ga.b, $gb.a, $gb.b], $theta);
    def one as float init $pa *
        (0.5 * conditionalOnAlleleTheta($f, $gb, $ga.a, $sa, $theta) +
        0.5 * conditionalOnAlleleTheta($f, $gb, $ga.b, $sa, $theta));
    def both as float init 0.0;
    if (sameGenotype($ga, $gb)) {
        $both = $pa;
    }
    return Terms{pa: $pa, pb: 0.0, none: $none, one: $one, both: $both};
}

# thetaJoint weights theta terms by a relationship.
func thetaJoint(t as Terms, rel as Relationship) {
    return $rel.k0 * $t.none + $rel.k1 * $t.one + $rel.k2 * $t.both;
}

/**
 * `jointGenotypeProbability` with the sub-population correction.
 *
 * Every allele not identical by descent is drawn from one sub-population by
 * the Balding-Nichols (1994) sampling formula,
 * `(m θ + (1-θ) p) / (1 + (n-1) θ)` for the next allele after `n` drawn of
 * which `m` match: four draws under IBD 0, three under IBD 1, two under IBD 2
 * (Fung & Hu 2008, ch. 5; Buckleton et al. 2016, ch. 4). At `theta == 0` it
 * returns exactly `jointGenotypeProbability`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the assumed relationship
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the joint probability
 * @throws {Error} kind "kinship" on an out-of-range theta, mismatched loci or
 *   a `rel` failing `relationship`'s checks, or kind "forensicgenetics" if
 *   the database does not cover the locus
 */
export func jointGenotypeProbabilityTheta(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship,
    theta as float) {
    checkTheta($theta);
    if ($theta == 0.0) {
        return jointGenotypeProbability($db, $ga, $gb, $rel);
    }
    checkPair($ga, $gb);
    checkRelationship($rel);
    return thetaJoint(thetaTerms(pairFreqs($db, $ga, $gb), $ga, $gb, $theta), $rel);
}

/**
 * `kinshipLRAt` with the sub-population correction: the ratio of
 * `jointGenotypeProbabilityTheta` under `rel` and `alt`. At `theta == 0` it
 * returns exactly `kinshipLRAt`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio at this locus
 * @throws {Error} kind "kinship" on an out-of-range theta, mismatched loci, a
 *   `rel` or `alt` failing `relationship`'s checks, or if `alt` gives the pair
 *   probability zero
 */
export func kinshipLRAtTheta(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as Relationship,
    alt as Relationship,
    theta as float) {
    checkTheta($theta);
    if ($theta == 0.0) {
        return kinshipLRAt($db, $ga, $gb, $rel, $alt);
    }
    checkPair($ga, $gb);
    checkRelationship($alt);
    checkRelationship($rel);
    def t as Terms init thetaTerms(pairFreqs($db, $ga, $gb), $ga, $gb, $theta);
    def under as float init thetaJoint($t, $alt);
    if ($under <= 0.0) {
        failZero($ga.locus, $alt);
    }
    return thetaJoint($t, $rel) / $under;
}

/**
 * `kinshipLR` with the sub-population correction: a relationship against
 * `UNRELATED`, the product of `kinshipLRAtTheta` over the shared loci. At
 * `theta == 0` it returns exactly `kinshipLR`.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship to test
 * @param db {FrequencyDb} the reference frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "kinship" on an out-of-range theta, if the profiles
 *   share no locus, share one the database does not cover, or if the LR
 *   exceeds the float range (use `log10KinshipLRTheta`)
 */
export func kinshipLRTheta(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    db as forensics.FrequencyDb,
    theta as float) {
    return kinshipLRAgainstTheta($a, $b, $rel, UNRELATED, $db, $theta);
}

/**
 * `kinshipLRAgainst` with the sub-population correction: one relationship
 * against another, the product of `kinshipLRAtTheta` over the shared loci. At
 * `theta == 0` it returns exactly `kinshipLRAgainst`.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "kinship" on an out-of-range theta, if the profiles
 *   share no locus, share one the database does not cover, if the
 *   alternative gives an observed pair probability zero, or if the LR exceeds
 *   the float range (use `log10KinshipLRAgainstTheta`)
 */
export func kinshipLRAgainstTheta(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    alt as Relationship,
    db as forensics.FrequencyDb,
    theta as float) {
    checkTheta($theta);
    if ($theta == 0.0) {
        return kinshipLRAgainst($a, $b, $rel, $alt, $db);
    }
    def loci as list of string init commonLoci($a, $b);
    if (len($loci) == 0) {
        fail("profiles '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    requireCovered($db, $loci, "profiles '" + $a.id + "' and '" + $b.id + "'");
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            kinshipLRAtTheta($db, $a.loci[$locus], $b.loci[$locus], $rel, $alt, $theta),
            "kinship LR",
            $a.id,
            $b.id,
            "; use log10KinshipLRTheta or log10KinshipLRAgainstTheta");
    }
    return $out;
}

/**
 * The base-10 logarithm of `kinshipLRTheta`, summed per locus.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the relationship to test
 * @param db {FrequencyDb} the reference frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "kinship" as for `kinshipLRTheta`, or if a locus has
 *   ratio 0 (an exclusion, which has no logarithm)
 */
export func log10KinshipLRTheta(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    db as forensics.FrequencyDb,
    theta as float) {
    return log10KinshipLRAgainstTheta($a, $b, $rel, UNRELATED, $db, $theta);
}

/**
 * The base-10 logarithm of `kinshipLRAgainstTheta`: the sum over the shared
 * loci of `log10(kinshipLRAtTheta)`. The other three `log10KinshipLR...`
 * functions are this with `UNRELATED` and / or theta 0.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "kinship" as for `kinshipLRAgainstTheta`, or if a locus
 *   has ratio 0 (an exclusion, which has no logarithm)
 */
export func log10KinshipLRAgainstTheta(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as Relationship,
    alt as Relationship,
    db as forensics.FrequencyDb,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init commonLoci($a, $b);
    if (len($loci) == 0) {
        fail("profiles '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    requireCovered($db, $loci, "profiles '" + $a.id + "' and '" + $b.id + "'");
    def out as float init 0.0;
    for (def locus in $loci) {
        def lr as float init kinshipLRAtTheta(
            $db,
            $a.loci[$locus],
            $b.loci[$locus],
            $rel,
            $alt,
            $theta);
        if ($lr <= 0.0) {
            fail("profiles '" + $a.id + "' and '" + $b.id + "': the ratio at locus '" + $locus +
                "' is 0 (an exclusion), which has no logarithm");
        }
        $out = $out + math.log10($lr);
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
 *   as `mutationProbability` on the model, or on a maternal exclusion (the
 *   mother can transmit neither of the child's alleles)
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
    for (def a in [$child.a, $child.b]) {
        $alleged[$a] = paternalTransmission($model, $locus, $father, $a);
    }
    def random as map of string to float init freqsOf($db, $locus, [$child.a, $child.b]).p;
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
    for (def locus in commonTrioLoci($mother, $child, $father)) {
        if (maps.has($db.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

# commonTrioLoci returns the loci all three profiles carry, in the child's
# order, whether or not the database covers them.
func commonTrioLoci(
    mother as forensics.Profile,
    child as forensics.Profile,
    father as forensics.Profile) {
    def out as list of string init [];
    for (def locus in $child.loci) {
        if (maps.has($mother.loci, $locus) and maps.has($father.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

# coveredTrioLoci is trioLoci with the guard that every locus the three share
# is in the database - a shared locus the table lacks is never dropped.
func coveredTrioLoci(
    db as forensics.FrequencyDb,
    mother as forensics.Profile,
    child as forensics.Profile,
    father as forensics.Profile) {
    def loci as list of string init commonTrioLoci($mother, $child, $father);
    if (len($loci) == 0) {
        fail("the mother, child and alleged father share no locus");
    }
    requireCovered(
        $db,
        $loci,
        "profiles '" + $mother.id + "', '" + $child.id + "' and '" + $father.id + "'");
    return $loci;
}

/**
 * The combined paternity index (CPI) for a trio, under `DEFAULT_MUTATION`.
 *
 * The product of `paternityIndexAt` over `trioLoci`. A locus not all three
 * profiles carry is skipped; a locus all three carry but the database does
 * not cover throws, since dropping it could hide an exclusion.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" if the trio shares no locus, shares one the
 *   database does not cover, on a maternal exclusion, or if the CPI exceeds
 *   the float range
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
 * @throws {Error} kind "kinship" if the trio shares no locus, shares one the
 *   database does not cover, on a maternal exclusion, or if the CPI exceeds
 *   the float range
 */
export func paternityIndexWith(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as MutationModel) {
    def loci as list of string init coveredTrioLoci($db, $mother, $child, $allegedFather);
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            paternityIndexAt(
                $db,
                $mother.loci[$locus],
                $child.loci[$locus],
                $allegedFather.loci[$locus],
                $model),
            "combined paternity index",
            $child.id,
            $allegedFather.id,
            "");
    }
    return $out;
}

# --- paternity with theta ---------------------------------------------------

# trioIndexTheta is the trio index with the random man's allele drawn by the
# sampling formula given the four typed alleles of the mother and the alleged
# father. X does not involve the population, so theta leaves it as it is.
func trioIndexTheta(
    db as forensics.FrequencyDb,
    mother as forensics.Genotype,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as MutationModel,
    theta as float) {
    def locus as string init $child.locus;
    def sampled as list of string init [$mother.a, $mother.b, $father.a, $father.b];
    def f as Freqs init freqsOf($db, $locus, [$child.a, $child.b]);
    def alleged as map of string to float init {};
    def random as map of string to float init {};
    for (def a in [$child.a, $child.b]) {
        $alleged[$a] = paternalTransmission($model, $locus, $father, $a);
        $random[$a] = drawProbability($f, $sampled, $a, $theta);
    }
    def y as float init childProbability($mother, $child, $random);
    if ($y <= 0.0) {
        fail("locus '" + $locus + "': the mother cannot have transmitted either of the child's" +
            " alleles - this is a maternal exclusion, not a paternity index");
    }
    return childProbability($mother, $child, $alleged) / $y;
}

/**
 * `paternityIndexAt` with the sub-population correction.
 *
 * `X` is unchanged: under paternity the child's paternal allele comes from
 * the alleged father, with mutation. In `Y` the random man's allele `c` is
 * drawn from the sub-population the typed mother and alleged father belong
 * to, `(m θ + (1-θ) p_c) / (1 + 3θ)` with `m` the copies of `c` among their
 * four alleles (Balding & Nichols 1995; Fung & Hu 2008, ch. 4; Buckleton et
 * al. 2016, ch. 11). At `theta == 0` it returns exactly `paternityIndexAt`.
 * @param db {FrequencyDb} the reference frequencies
 * @param mother {Genotype} the mother's genotype
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {MutationModel} the alleged father's mutation model
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "kinship" on an out-of-range theta, if the three
 *   genotypes are not at one locus, as `mutationProbability` on the model, or
 *   on a maternal exclusion
 */
export func paternityIndexAtTheta(
    db as forensics.FrequencyDb,
    mother as forensics.Genotype,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as MutationModel,
    theta as float) {
    checkTheta($theta);
    if ($theta == 0.0) {
        return paternityIndexAt($db, $mother, $child, $father, $model);
    }
    if ($mother.locus != $child.locus or $child.locus != $father.locus) {
        fail("the mother, child and alleged father genotypes must be at one locus");
    }
    return trioIndexTheta($db, $mother, $child, $father, $model, $theta);
}

/**
 * `paternityIndex` with the sub-population correction: the trio CPI under
 * `DEFAULT_MUTATION`. At `theta == 0` it returns exactly `paternityIndex`.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" on an out-of-range theta, if the trio shares
 *   no locus, shares one the database does not cover, on a maternal
 *   exclusion, or if the CPI exceeds the float range
 */
export func paternityIndexTheta(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    theta as float) {
    return paternityIndexWithTheta($mother, $child, $allegedFather, $db, DEFAULT_MUTATION, $theta);
}

/**
 * `paternityIndexWith` with the sub-population correction: the trio CPI, the
 * product of `paternityIndexAtTheta` over the trio's loci. At `theta == 0` it
 * returns exactly `paternityIndexWith`.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {MutationModel} the alleged father's mutation model
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" on an out-of-range theta, if the trio shares
 *   no locus, shares one the database does not cover, on a maternal
 *   exclusion, or if the CPI exceeds the float range
 */
export func paternityIndexWithTheta(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as MutationModel,
    theta as float) {
    checkTheta($theta);
    if ($theta == 0.0) {
        return paternityIndexWith($mother, $child, $allegedFather, $db, $model);
    }
    def loci as list of string init coveredTrioLoci($db, $mother, $child, $allegedFather);
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            paternityIndexAtTheta(
                $db,
                $mother.loci[$locus],
                $child.loci[$locus],
                $allegedFather.loci[$locus],
                $model,
                $theta),
            "combined paternity index",
            $child.id,
            $allegedFather.id,
            "");
    }
    return $out;
}

/**
 * The paternity index at one locus for a duo: alleged father and child, the
 * mother untyped.
 *
 * `X = P(child | alleged father)`: one allele from him (with mutation), the
 * other drawn from the sub-population given his two alleles.
 * `Y = P(child | alleged father, unrelated)`: both drawn, given his two and
 * then the first. With the sampling formula `D(c | S)`:
 *
 * ```
 *   X = F(c1) D(c2 | f) + F(c2) D(c1 | f)     (one term if c1 == c2)
 *   Y = k D(c1 | f) D(c2 | f, c1)             k = 2 if c1 != c2, else 1
 * ```
 *
 * With `NO_MUTATION` this is `kinshipLRAtTheta` for `PARENT_CHILD` against
 * `UNRELATED`. Theta 0 is allowed and gives the product-rule duo index.
 * @param db {FrequencyDb} the reference frequencies
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {MutationModel} the alleged father's mutation model
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "kinship" on an out-of-range theta, if the two
 *   genotypes are at different loci, as `mutationProbability` on the model,
 *   or if the child's genotype has probability zero under the alternative
 */
export func duoPaternityIndexAt(
    db as forensics.FrequencyDb,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as MutationModel,
    theta as float) {
    checkTheta($theta);
    if ($child.locus != $father.locus) {
        fail("the child and alleged father genotypes must be at one locus");
    }
    def locus as string init $child.locus;
    def sf as list of string init [$father.a, $father.b];
    def f as Freqs init freqsOf($db, $locus, [$child.a, $child.b]);
    def d1 as float init drawProbability($f, $sf, $child.a, $theta);
    def d2 as float init drawProbability($f, $sf, $child.b, $theta);
    def x as float init paternalTransmission($model, $locus, $father, $child.a) * $d2;
    def y as float init $d1 *
        drawProbability($f, [$father.a, $father.b, $child.a], $child.b, $theta);
    if ($child.a != $child.b) {
        $x = $x + paternalTransmission($model, $locus, $father, $child.b) * $d1;
        $y = 2.0 * $y;
    }
    if ($y <= 0.0) {
        fail("locus '" + $locus + "': the child's genotype has probability zero in the" +
            " reference population, so no ratio is defined");
    }
    return $x / $y;
}

/**
 * The combined paternity index for a duo (alleged father and child, mother
 * untyped): the product of `duoPaternityIndexAt` over the loci both carry. A
 * locus only one carries is skipped; one both carry but the database does not
 * cover throws.
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {MutationModel} the alleged father's mutation model
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the combined paternity index
 * @throws {Error} kind "kinship" on an out-of-range theta, if the two share no
 *   locus, share one the database does not cover, or if the index exceeds the
 *   float range
 */
export func duoPaternityIndex(
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as MutationModel,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init commonLoci($child, $allegedFather);
    if (len($loci) == 0) {
        fail("profiles '" + $child.id + "' and '" + $allegedFather.id + "' share no locus");
    }
    requireCovered($db, $loci, "profiles '" + $child.id + "' and '" + $allegedFather.id + "'");
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            duoPaternityIndexAt(
                $db,
                $child.loci[$locus],
                $allegedFather.loci[$locus],
                $model,
                $theta),
            "duo paternity index",
            $child.id,
            $allegedFather.id,
            "");
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
 * @throws {Error} kind "kinship" if the trio shares no locus, shares one the
 *   database does not cover, or on a maternal exclusion
 */
export func inconsistentLoci(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb) {
    def out as list of string init [];
    for (def locus in coveredTrioLoci($db, $mother, $child, $allegedFather)) {
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
 * @return {float} the probability of paternity, in [0, 1]; a CPI above about
 *   1e16 rounds to 1.0
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
