# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * Familias mutation models, silent alleles, and the trio, duo and pairwise
 * kinship indices that use them.
 *
 * A `MutationSpec` names a model: Equal, Proportional, Stepwise or Extended
 * stepwise, as the Familias R package `FamiliasLocus` builds them, optionally
 * made stationary by the PM transformation (Simonsson & Mostad 2016). At a
 * locus it becomes a `MutationMatrix` over the table's alleles plus the typed
 * ones: row `i`, column `j` is the probability that allele `i` is transmitted
 * as `j`, every row summing to 1. Unlike `kinship.MutationModel` (an
 * unbounded ladder, never normalised over a table), the matrix depends on
 * which alleles the locus has.
 *
 * A silent (null) allele is an allele with frequency `q` that no typing
 * shows: a single-peak phenotype `a` is genotype `aa` or `a/silent`, as in
 * Familias (`pers::execute` in its C++ source). Listed frequencies are read as
 * frequencies among detectable alleles and scaled by `1 - q`. No model
 * mutates to or from the silent allele.
 *
 * A `ParentageModel` holds a paternal and a maternal `MutationSpec` and the
 * per-locus silent frequencies. `paternityIndex` (trio) applies the maternal
 * model to the mother, the paternal one to the alleged father; the random man
 * of the alternative does not mutate. `duoPaternityIndex` uses the paternal
 * model only. `kinshipLR` adds silent alleles to the pairwise IBD form.
 * With no mutation and no silent allele every function here returns exactly
 * what the `kinship.j` function it extends returns.
 *
 * Sources: the Familias R package (Egeland, Kling & Mostad, "Relationship
 * Inference with Familias and R", Academic Press 2016) for the models and
 * silent alleles; Simonsson I, Mostad P (2016), "Stationary mutation models",
 * Forensic Sci. Int. Genet. 23:217-225, for the PM transformation, as
 * implemented in Familias' `stabilize` and pedmut's `makeStationary`.
 * `aabb2003Counts` is the AABB 2003 annual report summary as published by
 * NIST STRBase.
 * @module mutation
 * @example
 * import "@deckyard/forensicgenetics/mutation.j" as mutation;
 *
 * def male as mutation.MutationSpec init mutation.withLocusRates(
 *     mutation.stepwiseModel(0.002, 0.1),
 *     mutation.aabb2003Rates("total"));
 * def pm as mutation.ParentageModel init mutation.parentageModel($male, mutation.NO_MUTATION);
 * $pm = mutation.withSilentFrequency($pm, "D8S1179", 0.002);
 * def pi as float init mutation.paternityIndex($mother, $child, $af, $db, $pm, 0.0);
 */
use math;
use maps;
use linalg;
use convert;
use strings;

import "./forensicgenetics.j" as forensics;
import "./kinship.j" as kinship;

/**
 * The label of the silent allele in a `LocusAlleles` and `MutationMatrix`. A
 * typed genotype never carries it.
 */
export def const SILENT_ALLELE as string init "silent";

/**
 * Equal: every other detectable allele equally likely after a mutation.
 */
export def const MODEL_EQUAL as string init "equal";

/**
 * Proportional: a mutation's target drawn in proportion to its frequency.
 */
export def const MODEL_PROPORTIONAL as string init "proportional";

/**
 * Stepwise and extended stepwise: probability falling geometrically with the
 * repeat difference, microvariant changes at a separate rate.
 */
export def const MODEL_STEPWISE as string init "stepwise";

/**
 * A mutation model before it meets a locus's alleles.
 *
 * `rate` is the mutation rate (for Stepwise, of integer-step changes),
 * overridden per locus by `perLocus`. `range` and `rate2` are used by
 * Stepwise only. `stationary` applies the PM transformation.
 * @field model {string} `MODEL_EQUAL`, `MODEL_PROPORTIONAL` or `MODEL_STEPWISE`
 * @field rate {float} the per-meiosis mutation rate, in [0, 1)
 * @field range {float} Stepwise: the probability of `n+1` steps relative to `n`, above 0
 * @field rate2 {float} Stepwise: the rate of non-integer-step (microvariant) changes
 * @field stationary {bool} whether to make the matrix stationary (PM)
 * @field perLocus {map of string to float} locus -> rate, overriding `rate`
 */
export def struct MutationSpec {
    model as string,
    rate as float,
    range as float,
    rate2 as float,
    stationary as bool,
    perLocus as map of string to float
};

/**
 * No mutation: the identity matrix at every locus.
 */
export def const NO_MUTATION as MutationSpec init MutationSpec{
    model: "equal",
    rate: 0.0,
    range: 1.0,
    rate2: 0.0,
    stationary: false,
    perLocus: {}
};

/**
 * The alleles of one locus and their frequencies, in matrix order: the
 * table's alleles in insertion order, then typed alleles the table lacks,
 * then `SILENT_ALLELE` when the silent frequency is above 0.
 * @field locus {string} the locus name
 * @field alleles {list of string} the allele labels
 * @field frequencies {list of float} floored frequencies times `1 - q`; `q` for the silent allele
 */
export def struct LocusAlleles {
    locus as string,
    alleles as list of string,
    frequencies as list of float
};

/**
 * A mutation matrix at one locus: `rows[i][j]` is the probability that
 * `alleles[i]` is transmitted as `alleles[j]`. Every row sums to 1.
 * @field locus {string} the locus name
 * @field alleles {list of string} the row and column labels
 * @field rows {list of list of float} the transition probabilities
 */
export def struct MutationMatrix {
    locus as string,
    alleles as list of string,
    rows as list of list of float
};

/**
 * Sex-specific mutation models and per-locus silent allele frequencies for
 * parentage and kinship.
 * @field paternal {MutationSpec} the model for a father's transmissions
 * @field maternal {MutationSpec} the model for a mother's transmissions
 * @field silent {map of string to float} locus -> silent allele frequency; absent means 0
 */
export def struct ParentageModel {
    paternal as MutationSpec,
    maternal as MutationSpec,
    silent as map of string to float
};

/**
 * Mendelian transmission and no silent allele: the assumptions of
 * `kinship.paternityIndexAt` under `kinship.NO_MUTATION`.
 */
export def const MENDELIAN as ParentageModel init ParentageModel{
    paternal: NO_MUTATION,
    maternal: NO_MUTATION,
    silent: {}
};

/**
 * Mutations counted at one locus in paternity testing, as published.
 *
 * Maternal and paternal counts are mutations assigned to one parent;
 * `unassignedMutations` could be from either. Totals include all three.
 * @field locus {string} the locus, named as `panels.normalizeLocus` returns it
 * @field maternalMutations {int} mutations assigned to the mother
 * @field maternalMeioses {int} maternal meioses observed
 * @field paternalMutations {int} mutations assigned to the father
 * @field paternalMeioses {int} paternal meioses observed
 * @field unassignedMutations {int} mutations from either parent
 * @field totalMutations {int} all mutations
 * @field totalMeioses {int} all meioses
 */
export def struct MutationCount {
    locus as string,
    maternalMutations as int,
    maternalMeioses as int,
    paternalMutations as int,
    paternalMeioses as int,
    unassignedMutations as int,
    totalMutations as int,
    totalMeioses as int
};

/**
 * The most detectable alleles a locus may have under a stationary (PM) model;
 * `mutationMatrix` throws above it. One PM matrix takes about 0.19 s at 120
 * alleles and 0.47 s at 200. Real tables (SE33 included) list under 100.
 */
export def const MAX_STATIONARY_ALLELES as int init 200;

# ROW_TOLERANCE bounds how far a built row may sum from 1, as rounding.
def const ROW_TOLERANCE as float init 0.000000001;

# MAX_WEIGHT_LOG10 bounds |log10| of a Stepwise weight range^k (and of the
# row's weight sum), inside the normal float range.
def const MAX_WEIGHT_LOG10 as float init 300.0;

# MAX_REPEAT_VALUE bounds a Stepwise allele's repeat value: one decimal stays
# exact and step counts fit an int far below it.
def const MAX_REPEAT_VALUE as float init 1000000.0;

# FLOAT_MAX is the largest finite float; the across-loci products check
# against it before multiplying (the interpreter errors on overflow).
def const FLOAT_MAX as float init 1.7976931348623157e308;

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "mutation", message: $msg, file: "", line: 0, col: 0};
}

func checkRate(rate as float, what as string) {
    if ($rate < 0.0 or $rate >= 1.0) {
        fail($what + " must be in [0, 1)");
    }
}

func checkTheta(theta as float) {
    if ($theta < 0.0 or $theta >= 1.0) {
        fail("theta must be in [0, 1)");
    }
}

func checkSilent(q as float) {
    if ($q < 0.0 or $q >= 1.0) {
        fail("silent allele frequency must be in [0, 1)");
    }
}

# timesLR multiplies one locus's factor into a running product, throwing
# where it would leave the float range (kinship.j's guard). `acc <= 1` cannot
# overflow; for `acc > 1`, `FLOAT_MAX / acc` is finite.
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

# --- model specifications ---------------------------------------------------

func spec(model as string, rate as float, range as float, rate2 as float) {
    return MutationSpec{
        model: $model,
        rate: $rate,
        range: $range,
        rate2: $rate2,
        stationary: false,
        perLocus: {}
    };
}

/**
 * The Equal model: `1 - rate` on the diagonal, `rate / (n - 1)` to each of the
 * other `n - 1` detectable alleles.
 * @param rate {float} the mutation rate, in [0, 1)
 * @return {MutationSpec} the model
 * @throws {Error} kind "mutation" on an out-of-range rate
 */
export func equalModel(rate as float) {
    checkRate($rate, "mutation rate");
    return spec(MODEL_EQUAL, $rate, 1.0, 0.0);
}

/**
 * The Proportional model: `alpha * p_j` to allele `j`, diagonal
 * `1 - alpha + alpha * p_i`, with `p` normalised over the detectable alleles
 * and `alpha = rate / S / Σ p_j (1 - p_j)`, `S` their summed frequency, so the
 * expected rate is `rate`. Stationary by construction.
 * @param rate {float} the expected mutation rate, in [0, 1)
 * @return {MutationSpec} the model
 * @throws {Error} kind "mutation" on an out-of-range rate
 */
export func proportionalModel(rate as float) {
    checkRate($rate, "mutation rate");
    return spec(MODEL_PROPORTIONAL, $rate, 1.0, 0.0);
}

/**
 * The Stepwise model: from allele `i`, an integer-step change to `j` has
 * probability `rate * range^|i - j|`, normalised over the other alleles of
 * `i`'s microvariant class. With only integer alleles this is Familias 2's
 * "probability decreasing with range (equal)". Microvariants never mutate
 * across classes (`rate2` 0).
 * @param rate {float} the integer-step mutation rate, in [0, 1)
 * @param range {float} the probability of `n+1` steps relative to `n`, above 0
 * @return {MutationSpec} the model
 * @throws {Error} kind "mutation" on an out-of-range rate or range
 */
export func stepwiseModel(rate as float, range as float) {
    return extendedStepwiseModel($rate, $range, 0.0);
}

/**
 * The Extended stepwise model (Familias 3): Stepwise, plus changes between
 * microvariant classes (`15` -> `15.2`) at total rate `rate2`, split equally
 * over the alleles of the other classes. An allele alone in its class keeps
 * `1 - rate2`; one in a class with every allele keeps `1 - rate`; otherwise
 * `1 - rate - rate2`.
 * @param rate {float} the integer-step mutation rate, in [0, 1)
 * @param range {float} the probability of `n+1` steps relative to `n`, above 0
 * @param rate2 {float} the non-integer-step mutation rate, in [0, 1)
 * @return {MutationSpec} the model
 * @throws {Error} kind "mutation" on an out-of-range rate, range or rate2, or
 *   if `rate + rate2` exceeds 1
 */
export func extendedStepwiseModel(rate as float, range as float, rate2 as float) {
    checkRate($rate, "mutation rate");
    checkRate($rate2, "microvariant mutation rate");
    if ($range <= 0.0) {
        fail("mutation range must be above 0");
    }
    if ($rate + $rate2 > 1.0) {
        fail("mutation rate plus microvariant rate must be at most 1");
    }
    return spec(MODEL_STEPWISE, $rate, $range, $rate2);
}

/**
 * A copy of `model` made stationary by the PM transformation of Simonsson &
 * Mostad (2016): each row `i` off the diagonal is scaled by
 * `d_i = R v_i / (c p_i)`, `v` the matrix's stationary distribution, `R` and
 * `c` the expected mutation rates under `p` and `v`. The result keeps `p`
 * (the locus frequencies) stationary and the expected rate `R`. Building the
 * matrix throws when a scaled diagonal would go negative or the locus has
 * more than `MAX_STATIONARY_ALLELES` detectable alleles.
 * @param model {MutationSpec} the model to transform
 * @return {MutationSpec} the stationary copy
 */
export func stationaryModel(model as MutationSpec) {
    def out as MutationSpec init $model;
    $out.stationary = true;
    return $out;
}

/**
 * A copy of `model` with one locus's rate overridden (for Stepwise, the
 * integer-step rate).
 * @param model {MutationSpec} the model to copy
 * @param locus {string} the locus name
 * @param rate {float} that locus's mutation rate, in [0, 1)
 * @return {MutationSpec} the updated copy
 * @throws {Error} kind "mutation" on an out-of-range rate, or a Stepwise rate
 *   that with `rate2` exceeds 1
 */
export func withLocusRate(model as MutationSpec, locus as string, rate as float) {
    checkRate($rate, "mutation rate");
    if ($rate + $model.rate2 > 1.0) {
        fail("locus '" + $locus + "': mutation rate plus microvariant rate must be at most 1");
    }
    def out as MutationSpec init $model;
    $out.perLocus[$locus] = $rate;
    return $out;
}

/**
 * A copy of `model` with the rates of several loci overridden, for example
 * from `aabb2003Rates`.
 * @param model {MutationSpec} the model to copy
 * @param rates {map of string to float} locus -> mutation rate
 * @return {MutationSpec} the updated copy
 * @throws {Error} kind "mutation" on an out-of-range rate, or a Stepwise rate
 *   that with `rate2` exceeds 1
 */
export func withLocusRates(model as MutationSpec, rates as map of string to float) {
    def out as MutationSpec init $model;
    for (def locus in $rates) {
        def rate as float init $rates[$locus];
        checkRate($rate, "locus '" + $locus + "': mutation rate");
        if ($rate + $model.rate2 > 1.0) {
            fail("locus '" + $locus + "': mutation rate plus microvariant rate must be at most 1");
        }
        $out.perLocus[$locus] = $rate;
    }
    return $out;
}

# rateAt returns the rate in force at a locus.
func rateAt(model as MutationSpec, locus as string) {
    if (maps.has($model.perLocus, $locus)) {
        return $model.perLocus[$locus];
    }
    return $model.rate;
}

# --- matrices ---------------------------------------------------------------

/**
 * The allele set and frequencies of one locus: the table's alleles, typed
 * alleles it lacks (at the `5/(2N)` floor), and the silent allele when
 * `silent` is above 0. Detectable frequencies are floored table values times
 * `1 - silent`.
 * @param db {FrequencyDb} the reference frequencies
 * @param locus {string} the locus name
 * @param typed {list of string} alleles seen in the case
 * @param silent {float} the silent allele frequency, in [0, 1)
 * @return {LocusAlleles} the alleles and frequencies
 * @throws {Error} kind "mutation" on an out-of-range silent frequency, a
 *   locus the table lacks, or a typed `SILENT_ALLELE`; kind
 *   "forensicgenetics" for an allele with no frequency and no floor
 */
export func locusAlleles(
    db as forensics.FrequencyDb,
    locus as string,
    typed as list of string,
    silent as float) {
    checkSilent($silent);
    if (not forensics.hasLocus($db, $locus)) {
        fail("frequency database '" + $db.name + "' does not cover locus '" + $locus + "'");
    }
    def labels as list of string init [];
    def seen as map of string to bool init {};
    for (def a in $db.loci[$locus].alleles) {
        $labels[] = $a;
        $seen[$a] = true;
    }
    for (def a in $typed) {
        if ($a == SILENT_ALLELE) {
            fail("locus '" + $locus + "': a typed genotype cannot carry the silent allele");
        }
        if (not maps.has($seen, $a)) {
            $labels[] = $a;
            $seen[$a] = true;
        }
    }
    def freqs as list of float init [];
    for (def a in $labels) {
        $freqs[] = (1.0 - $silent) * forensics.alleleFrequency($db, $locus, $a);
    }
    if ($silent > 0.0) {
        $labels[] = SILENT_ALLELE;
        $freqs[] = $silent;
    }
    return LocusAlleles{locus: $locus, alleles: $labels, frequencies: $freqs};
}

# detectableCount is the number of alleles before a trailing silent allele.
func detectableCount(alleles as list of string) {
    def n as int init len($alleles);
    if ($n > 0 and $alleles[$n - 1] == SILENT_ALLELE) {
        return $n - 1;
    }
    return $n;
}

func identityRows(n as int) {
    def rows as list of list of float init [];
    for (def i in 0..$n) {
        def row as list of float init [];
        for (def j in 0..$n) {
            if ($i == $j) {
                $row[] = 1.0;
            } else {
                $row[] = 0.0;
            }
        }
        $rows[] = $row;
    }
    return $rows;
}

/**
 * The mutation matrix of `model` over a locus's alleles. The silent allele's
 * row and column are those of the identity. A stationary `model` is
 * transformed by PM over the detectable alleles.
 * @param model {MutationSpec} the mutation model
 * @param alleles {LocusAlleles} the locus's alleles and frequencies
 * @return {MutationMatrix} the matrix
 * @throws {Error} kind "mutation" on an unknown model or a rate in force at
 *   the locus outside [0, 1) (a `MutationSpec` literal); for Stepwise, a
 *   `rate2` outside [0, 1), `rate + rate2` above 1, a `range` not above 0, an
 *   allele that is not a repeat count with at most one decimal, a `.0` label,
 *   a repeat value above 1e6, two labels with one repeat value, or a
 *   `range^k` outside the float range; an entry outside [0, 1] (a
 *   Proportional rate too high for the frequencies), a row that does not sum
 *   to 1; for a stationary model, more than `MAX_STATIONARY_ALLELES`
 *   detectable alleles, or a PM transformation that does not exist or has no
 *   unique stationary distribution
 */
export func mutationMatrix(model as MutationSpec, alleles as LocusAlleles) {
    def n as int init len($alleles.alleles);
    def nd as int init detectableCount($alleles.alleles);
    def rate as float init rateAt($model, $alleles.locus);
    checkSpec($model, $rate, $alleles.locus);
    def m as MutationMatrix init MutationMatrix{
        locus: $alleles.locus,
        alleles: $alleles.alleles,
        rows: []
    };
    if ($nd <= 1 or ($rate <= 0.0 and $model.rate2 <= 0.0)) {
        $m.rows = identityRows($n);
        return $m;
    }
    if ($model.stationary and $nd > MAX_STATIONARY_ALLELES) {
        fail("locus '" + $alleles.locus + "': {$nd} detectable alleles under a stationary" +
            " model; MAX_STATIONARY_ALLELES is {MAX_STATIONARY_ALLELES}");
    }
    if ($model.model == MODEL_EQUAL) {
        $m.rows = equalRows($n, $nd, $rate);
    } elseif ($model.model == MODEL_PROPORTIONAL) {
        $m.rows = proportionalRows($n, $nd, $rate, $alleles.frequencies);
    } else {
        $m.rows = stepwiseRows($n, $nd, $rate, $model, $alleles);
    }
    checkRows($m);
    if ($model.stationary) {
        $m.rows = stabilisePM($m, $nd, $alleles);
        checkRows($m);
    }
    return $m;
}

# checkSpec validates what a `MutationSpec` literal can carry past the
# constructors: the model name and the parameters in force at the locus.
func checkSpec(model as MutationSpec, rate as float, locus as string) {
    if ($model.model != MODEL_EQUAL and $model.model != MODEL_PROPORTIONAL and
        $model.model != MODEL_STEPWISE) {
        fail("unknown mutation model '" + $model.model + "'");
    }
    checkRate($rate, "locus '" + $locus + "': mutation rate");
    if ($model.model != MODEL_STEPWISE) {
        return;
    }
    checkRate($model.rate2, "locus '" + $locus + "': microvariant mutation rate");
    if ($model.range <= 0.0) {
        fail("locus '" + $locus + "': mutation range must be above 0");
    }
    if ($rate + $model.rate2 > 1.0) {
        fail("locus '" + $locus + "': mutation rate plus microvariant rate must be at most 1");
    }
}

func equalRows(n as int, nd as int, rate as float) {
    def rows as list of list of float init identityRows($n);
    def off as float init $rate / convert.toFloat($nd - 1);
    for (def i in 0..$nd) {
        for (def j in 0..$nd) {
            if ($i == $j) {
                $rows[$i][$j] = 1.0 - $rate;
            } else {
                $rows[$i][$j] = $off;
            }
        }
    }
    return $rows;
}

# proportionalRows follows FamiliasLocus: frequencies normalised over the
# detectable alleles, alpha divided by their sum so the expected rate over all
# alleles (silent included, which never mutates) is `rate`.
func proportionalRows(n as int, nd as int, rate as float, freqs as list of float) {
    def rows as list of list of float init identityRows($n);
    def sum as float init 0.0;
    for (def i in 0..$nd) {
        $sum = $sum + $freqs[$i];
    }
    def frq as list of float init [];
    def spread as float init 0.0;
    for (def i in 0..$nd) {
        def f as float init $freqs[$i] / $sum;
        $frq[] = $f;
        $spread = $spread + $f * (1.0 - $f);
    }
    def alpha as float init $rate / $sum / $spread;
    for (def i in 0..$nd) {
        for (def j in 0..$nd) {
            if ($i == $j) {
                $rows[$i][$j] = 1.0 - $alpha + $alpha * $frq[$j];
            } else {
                $rows[$i][$j] = $alpha * $frq[$j];
            }
        }
    }
    return $rows;
}

# microClass is FamiliasLocus's `microgroup`, `round((v - round(v)) * 10)`
# with R's round (half to even): 15 -> 0, 15.2 -> 2, 15.7 -> -3, and
# 14.5 / 16.5 -> 5 but 15.5 -> -5.
func microClass(value as float) {
    return math.round(($value - convert.toFloat(roundHalfEven($value))) * 10.0);
}

# roundHalfEven is R's round(x): halves go to the even neighbour.
# math.round takes them away from zero.
func roundHalfEven(value as float) {
    def lo as int init math.floor($value);
    if ($value - convert.toFloat($lo) == 0.5) {
        if ($lo % 2 == 0) {
            return $lo;
        }
        return $lo + 1;
    }
    return math.round($value);
}

func stepwiseValue(locus as string, allele as string) {
    if (not forensics.isNumericAllele($allele)) {
        fail("locus '" + $locus + "': the stepwise model needs repeat counts, not '" +
            $allele + "'");
    }
    def parts as list of string init strings.split($allele, ".");
    if (len($parts) == 2 and len($parts[1]) != 1) {
        fail("locus '" + $locus + "': allele '" + $allele +
            "' - a microvariant takes exactly one decimal");
    }
    if (len($parts) == 2 and $parts[1] == "0") {
        fail("locus '" + $locus + "': allele '" + $allele + "' - write the integer label '" +
            $parts[0] + "', not a '.0' decimal");
    }
    def v as float init forensics.alleleValue($allele);
    if (math.abs($v) > MAX_REPEAT_VALUE) {
        fail("locus '" + $locus + "': allele '" + $allele +
            "' - the stepwise model takes repeat values up to 1e6");
    }
    return $v;
}

# checkStepwiseAlleles rejects two labels with one repeat value ("10" and
# "010"): the matrix would treat them as distinct alleles zero steps apart.
# It also rejects a `range` whose largest within-class weight `range^k` (or
# the row's weight sum) leaves the finite, normal float range.
func checkStepwiseAlleles(
    alleles as LocusAlleles,
    nd as int,
    range as float,
    values as list of float,
    classes as list of int) {
    if ($range <= 0.0) {
        fail("locus '" + $alleles.locus + "': mutation range must be above 0");
    }
    def far as int init 0;
    def farI as int init 0;
    def farJ as int init 0;
    for (def i in 0..$nd) {
        for (def j in $i + 1..$nd) {
            if ($values[$i] == $values[$j]) {
                fail("locus '" + $alleles.locus + "': alleles '" + $alleles.alleles[$i] +
                    "' and '" + $alleles.alleles[$j] + "' have the same repeat value");
            }
            if ($classes[$i] != $classes[$j]) {
                continue;
            }
            def steps as int init math.round(math.abs($values[$i] - $values[$j]));
            if ($steps > $far) {
                $far = $steps;
                $farI = $i;
                $farJ = $j;
            }
        }
    }
    if ($far == 0) {
        return;
    }
    def lw as float init convert.toFloat($far) * math.log10($range);
    if ($lw < -MAX_WEIGHT_LOG10 or $lw + math.log10(convert.toFloat($nd)) > MAX_WEIGHT_LOG10) {
        fail("locus '" + $alleles.locus + "': range^{$far} between alleles '" +
            $alleles.alleles[$farI] + "' and '" + $alleles.alleles[$farJ] + "' is 10^{$lw}," +
            " outside the float range - the step weight would be 0 or infinite");
    }
}

# stepwiseRows follows FamiliasLocus's 'Stepwise' model: weights
# range^steps within a microvariant class, scaled to sum to `rate`; rate2
# spread equally over the alleles of other classes.
func stepwiseRows(
    n as int,
    nd as int,
    rate as float,
    model as MutationSpec,
    alleles as LocusAlleles) {
    def rows as list of list of float init [];
    def values as list of float init [];
    def classes as list of int init [];
    for (def i in 0..$nd) {
        def v as float init stepwiseValue($alleles.locus, $alleles.alleles[$i]);
        $values[] = $v;
        $classes[] = microClass($v);
    }
    checkStepwiseAlleles($alleles, $nd, $model.range, $values, $classes);
    for (def i in 0..$nd) {
        $rows[] = stepwiseRow($n, $nd, $i, $rate, $model, $values, $classes);
    }
    for (def i in $nd..$n) {
        def row as list of float init [];
        for (def j in 0..$n) {
            $row[] = 0.0;
        }
        $row[$i] = 1.0;
        $rows[] = $row;
    }
    return $rows;
}

# stepwiseRow is row i over all n alleles: the in-class weights computed once
# into the row, then scaled; 0 in the silent allele's column.
func stepwiseRow(
    n as int,
    nd as int,
    i as int,
    rate as float,
    model as MutationSpec,
    values as list of float,
    classes as list of int) {
    def same as int init 0;
    def weight as float init 0.0;
    def row as list of float init [];
    for (def j in 0..$nd) {
        def w as float init 0.0;
        if ($classes[$j] == $classes[$i]) {
            $same = $same + 1;
            if ($j != $i) {
                $w = stepWeight($model.range, $values[$i], $values[$j]);
                $weight = $weight + $w;
            }
        }
        $row[] = $w;
    }
    def across as float init 0.0;
    if ($same < $nd) {
        $across = $model.rate2 / convert.toFloat($nd - $same);
    }
    for (def j in 0..$nd) {
        if ($j == $i) {
            $row[$j] = stepwiseDiagonal($same, $nd, $rate, $model.rate2);
        } elseif ($classes[$j] == $classes[$i]) {
            $row[$j] = $row[$j] / $weight * $rate;
        } else {
            $row[$j] = $across;
        }
    }
    for (def j in $nd..$n) {
        $row[] = 0.0;
    }
    return $row;
}

# stepWeight is range^k for a change of k repeats within one class (k an
# integer up to rounding of the decimal labels).
func stepWeight(range as float, a as float, b as float) {
    def steps as int init math.round(math.abs($a - $b));
    return math.pow($range, convert.toFloat($steps));
}

func stepwiseDiagonal(same as int, nd as int, rate as float, rate2 as float) {
    if ($same == $nd) {
        return 1.0 - $rate;
    }
    if ($same == 1) {
        return 1.0 - $rate2;
    }
    return 1.0 - $rate - $rate2;
}

# checkRows enforces what FamiliasLocus's validateMutation does: entries in
# [0, 1], rows summing to 1.
func checkRows(m as MutationMatrix) {
    for (def i in 0..len($m.rows)) {
        def sum as float init 0.0;
        for (def x in $m.rows[$i]) {
            if ($x < 0.0 or $x > 1.0) {
                fail("locus '" + $m.locus +
                    "': the mutation matrix has an entry outside [0, 1] - the rate is too high" +
                    " for this model and these frequencies");
            }
            $sum = $sum + $x;
        }
        if (math.abs($sum - 1.0) > ROW_TOLERANCE) {
            fail("locus '" + $m.locus + "': a mutation matrix row does not sum to 1");
        }
    }
}

# stabilisePM is pedmut's .stabilizePM (Familias' stabilize, method "PM") on
# the detectable block, with p the detectable frequencies normalised to 1. It
# returns fresh rows; the silent allele's row and column are copied.
func stabilisePM(m as MutationMatrix, nd as int, alleles as LocusAlleles) {
    def n as int init len($m.rows);
    def sum as float init 0.0;
    for (def i in 0..$nd) {
        $sum = $sum + $alleles.frequencies[$i];
    }
    def r as float init 0.0;
    for (def i in 0..$nd) {
        $r = $r + $m.rows[$i][$i] * $alleles.frequencies[$i] / $sum;
    }
    $r = 1.0 - $r;
    def out as list of list of float init [];
    if ($r <= 0.0) {
        for (def i in 0..$n) {
            $out[] = $m.rows[$i];
        }
        return $out;
    }
    def v as list of float init stationaryVector($m, $nd, $alleles);
    def c as float init 0.0;
    for (def i in 0..$nd) {
        $c = $c + $m.rows[$i][$i] * $v[$i];
    }
    $c = 1.0 - $c;
    for (def i in 0..$n) {
        if ($i >= $nd) {
            $out[] = $m.rows[$i];
            continue;
        }
        def p as float init $alleles.frequencies[$i] / $sum;
        def d as float init $r * $v[$i] / ($c * $p);
        if ($d * (1.0 - $m.rows[$i][$i]) > 1.0) {
            fail("locus '" + $alleles.locus +
                "': the PM stationary transformation does not exist for this model");
        }
        def row as list of float init [];
        for (def j in 0..$n) {
            if ($j == $i) {
                $row[] = 1.0 - $d * (1.0 - $m.rows[$i][$j]);
            } elseif ($j < $nd) {
                $row[] = $d * $m.rows[$i][$j];
            } else {
                $row[] = $m.rows[$i][$j];
            }
        }
        $out[] = $row;
    }
    return $out;
}

# stationaryVector solves v M = v, Σ v = 1 over the first nd alleles: the
# system (M^T - I) v = 0 with its last equation replaced by Σ v = 1, by
# linalg.solve (Gauss-Jordan with partial pivoting, in Go: the O(n³) loop
# interpreted took 2.5 s at 120 alleles). Uniqueness is checked on the
# matrix's zero pattern first (`reachesFrom`), since the solve only rejects
# an exactly zero pivot.
func stationaryVector(m as MutationMatrix, nd as int, alleles as LocusAlleles) {
    def hub as int init 0;
    for (def i in 1..$nd) {
        if ($alleles.frequencies[$i] > $alleles.frequencies[$hub]) {
            $hub = $i;
        }
    }
    if (not reachesFrom($m, $nd, $hub)) {
        noStationary($m.locus);
    }
    def a as list of list of float init [];
    def b as list of float init [];
    for (def r in 0..$nd) {
        def eq as list of float init [];
        for (def c in 0..$nd) {
            if ($r == $nd - 1) {
                $eq[] = 1.0;
            } elseif ($r == $c) {
                $eq[] = $m.rows[$c][$r] - 1.0;
            } else {
                $eq[] = $m.rows[$c][$r];
            }
        }
        $a[] = $eq;
        if ($r == $nd - 1) {
            $b[] = 1.0;
        } else {
            $b[] = 0.0;
        }
    }
    def v as list of float init [];
    try {
        $v = linalg.solve($a, $b);
    } catch (e) {
        noStationary($m.locus);
    }
    for (def x in $v) {
        if ($x <= 0.0) {
            noStationary($m.locus);
        }
    }
    return $v;
}

func noStationary(locus as string) {
    fail("locus '" + $locus + "': the mutation matrix has no unique stationary" +
        " distribution, so PM cannot be applied");
}

# reachesFrom reports whether every detectable allele can mutate, in some
# number of steps, to allele `hub`. Then every closed class contains `hub`,
# so there is one and the stationary distribution is unique. The builders'
# zero patterns are symmetric for positive frequencies (no transient allele),
# so with `hub` the most frequent allele the converse holds too: false means
# two closed classes (Stepwise with microvariants and `rate2` 0).
func reachesFrom(m as MutationMatrix, nd as int, hub as int) {
    def reached as list of bool init [];
    for (def i in 0..$nd) {
        $reached[] = $i == $hub;
    }
    def queue as list of int init [$hub];
    def head as int init 0;
    while ($head < len($queue)) {
        def j as int init $queue[$head];
        $head = $head + 1;
        for (def i in 0..$nd) {
            if (not $reached[$i] and $m.rows[$i][$j] > 0.0) {
                $reached[$i] = true;
                $queue[] = $i;
            }
        }
    }
    return len($queue) == $nd;
}

# indexOfAllele maps labels to matrix positions.
func indexOfAllele(alleles as list of string) {
    def out as map of string to int init {};
    for (def i in 0..len($alleles)) {
        $out[$alleles[$i]] = $i;
    }
    return $out;
}

/**
 * One entry of a mutation matrix: the probability that `source` is
 * transmitted as `target`.
 * @param m {MutationMatrix} the matrix
 * @param source {string} the parent's allele
 * @param target {string} the transmitted allele
 * @return {float} the transition probability
 * @throws {Error} kind "mutation" if either allele is not in the matrix
 */
export func transitionProbability(m as MutationMatrix, source as string, target as string) {
    def idx as map of string to int init indexOfAllele($m.alleles);
    for (def a in [$source, $target]) {
        if (not maps.has($idx, $a)) {
            fail("locus '" + $m.locus + "': allele '" + $a + "' is not in the mutation matrix");
        }
    }
    return $m.rows[$idx[$source]][$idx[$target]];
}

func requireSameAlleles(m as MutationMatrix, alleles as LocusAlleles) {
    if ($m.alleles != $alleles.alleles) {
        fail("locus '" + $m.locus +
            "': the matrix and the allele frequencies list different alleles");
    }
}

/**
 * The expected mutation rate `Σ π_i (1 - M_ii)` under the locus frequencies.
 * @param m {MutationMatrix} the matrix
 * @param alleles {LocusAlleles} the alleles it was built over
 * @return {float} the expected rate
 * @throws {Error} kind "mutation" if the two list different alleles
 */
export func overallRate(m as MutationMatrix, alleles as LocusAlleles) {
    requireSameAlleles($m, $alleles);
    def out as float init 0.0;
    for (def i in 0..len($m.rows)) {
        $out = $out + $alleles.frequencies[$i] * (1.0 - $m.rows[$i][$i]);
    }
    return $out;
}

/**
 * Whether the locus frequencies are stationary under the matrix,
 * `Σ_i π_i M_ij = π_j` for every `j` to within 1e-9.
 * @param m {MutationMatrix} the matrix
 * @param alleles {LocusAlleles} the alleles it was built over
 * @return {bool} true when one generation of mutation keeps the frequencies
 * @throws {Error} kind "mutation" if the two list different alleles
 */
export func isStationary(m as MutationMatrix, alleles as LocusAlleles) {
    requireSameAlleles($m, $alleles);
    def n as int init len($m.rows);
    for (def j in 0..$n) {
        def after as float init 0.0;
        for (def i in 0..$n) {
            $after = $after + $alleles.frequencies[$i] * $m.rows[$i][$j];
        }
        if (math.abs($after - $alleles.frequencies[$j]) > ROW_TOLERANCE) {
            return false;
        }
    }
    return true;
}

# --- parentage models -------------------------------------------------------

/**
 * A parentage model with separate paternal and maternal mutation models and
 * no silent allele.
 * @param paternal {MutationSpec} the model for a father's transmissions
 * @param maternal {MutationSpec} the model for a mother's transmissions
 * @return {ParentageModel} the model
 */
export func parentageModel(paternal as MutationSpec, maternal as MutationSpec) {
    return ParentageModel{paternal: $paternal, maternal: $maternal, silent: {}};
}

/**
 * A copy of `model` with a silent allele of frequency `q` at one locus.
 * @param model {ParentageModel} the model to copy
 * @param locus {string} the locus name
 * @param q {float} the silent allele frequency, in [0, 1)
 * @return {ParentageModel} the updated copy
 * @throws {Error} kind "mutation" on an out-of-range frequency
 */
export func withSilentFrequency(model as ParentageModel, locus as string, q as float) {
    checkSilent($q);
    def out as ParentageModel init $model;
    $out.silent[$locus] = $q;
    return $out;
}

func silentAt(silent as map of string to float, locus as string) {
    if (maps.has($silent, $locus)) {
        return $silent[$locus];
    }
    return 0.0;
}

# expand lists the true genotypes a typed genotype can be: itself, and for a
# single peak `a` also `a/silent` when there is a silent allele.
func expand(g as forensics.Genotype, q as float) {
    def out as list of forensics.Genotype init [$g];
    if ($q > 0.0 and $g.a == $g.b) {
        $out[] = forensics.Genotype{locus: $g.locus, a: $g.a, b: SILENT_ALLELE};
    }
    return $out;
}

func orderings(g as forensics.Genotype) {
    if ($g.a == $g.b) {
        return 1.0;
    }
    return 2.0;
}

# draw is the Balding-Nichols sampling formula over the locus frequencies, in
# the same arithmetic as kinship.j's drawProbability.
func draw(
    alleles as LocusAlleles,
    idx as map of string to int,
    sampled as list of string,
    allele as string,
    theta as float) {
    def m as int init 0;
    for (def s in $sampled) {
        if ($s == $allele) {
            $m = $m + 1;
        }
    }
    def n as float init convert.toFloat(len($sampled));
    def p as float init $alleles.frequencies[$idx[$allele]];
    return (convert.toFloat($m) * $theta + (1.0 - $theta) * $p) / (1.0 + ($n - 1.0) * $theta);
}

func sequence(
    alleles as LocusAlleles,
    idx as map of string to int,
    xs as list of string,
    theta as float) {
    def out as float init 1.0;
    def seen as list of string init [];
    for (def a in $xs) {
        $out = $out * draw($alleles, $idx, $seen, $a, $theta);
        $seen[] = $a;
    }
    return $out;
}

# transmit is a parent's transmission probability with mutation:
# ½ M[g.a][allele] + ½ M[g.b][allele].
func transmit(
    m as MutationMatrix,
    idx as map of string to int,
    g as forensics.Genotype,
    allele as string) {
    def j as int init $idx[$allele];
    return 0.5 * $m.rows[$idx[$g.a]][$j] + 0.5 * $m.rows[$idx[$g.b]][$j];
}

# childProbability is kinship.j's: P(child | mother's and father's transmission
# maps), the same arithmetic.
func childProbability(
    fromMother as map of string to float,
    cg as forensics.Genotype,
    fromFather as map of string to float) {
    if ($cg.a == $cg.b) {
        return $fromMother[$cg.a] * $fromFather[$cg.a];
    }
    return $fromMother[$cg.a] * $fromFather[$cg.b] + $fromMother[$cg.b] * $fromFather[$cg.a];
}

func requireNotSilent(gs as list of forensics.Genotype) {
    for (def g in $gs) {
        if ($g.a == SILENT_ALLELE or $g.b == SILENT_ALLELE) {
            fail("locus '" + $g.locus + "': a typed genotype cannot carry the silent allele");
        }
    }
}

# --- trio -------------------------------------------------------------------

# trioTerms adds one (mother, father) true-genotype pair's contribution to X
# and Y, returned as [x, y].
func trioTerms(
    alleles as LocusAlleles,
    idx as map of string to int,
    mm as MutationMatrix,
    fm as MutationMatrix,
    gm as forensics.Genotype,
    gf as forensics.Genotype,
    children as list of forensics.Genotype,
    theta as float) {
    def sampled as list of string init [$gm.a, $gm.b, $gf.a, $gf.b];
    def x as float init 0.0;
    def y as float init 0.0;
    for (def gc in $children) {
        def maternal as map of string to float init {};
        def alleged as map of string to float init {};
        def random as map of string to float init {};
        for (def a in [$gc.a, $gc.b]) {
            $maternal[$a] = transmit($mm, $idx, $gm, $a);
            $alleged[$a] = transmit($fm, $idx, $gf, $a);
            $random[$a] = draw($alleles, $idx, $sampled, $a, $theta);
        }
        $x = $x + childProbability($maternal, $gc, $alleged);
        $y = $y + childProbability($maternal, $gc, $random);
    }
    return [$x, $y];
}

/**
 * The trio paternity index at one locus under a parentage model.
 *
 * `X = P(child | mother, alleged father)`, the mother transmitting through
 * the maternal matrix and the alleged father through the paternal one;
 * `Y = P(child | mother, random man)`, the random man's allele drawn from the
 * sub-population given the four parental alleles (Balding-Nichols), without
 * mutation. With a silent allele, single-peak phenotypes are summed over
 * `aa` and `a/silent`, each parental pair weighted by its joint probability.
 * With `MENDELIAN` (or both rates 0 and no silent allele) it returns exactly
 * `kinship.paternityIndexAtTheta` under `kinship.NO_MUTATION`.
 * @param db {FrequencyDb} the reference frequencies
 * @param mother {Genotype} the mother's genotype
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {ParentageModel} the mutation models and silent frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "mutation" on an out-of-range theta, genotypes at
 *   different loci, a typed `SILENT_ALLELE`, or a child genotype impossible
 *   given the mother (a maternal exclusion when the maternal model has no
 *   mutation); see `mutationMatrix` for model errors
 */
export func paternityIndexAt(
    db as forensics.FrequencyDb,
    mother as forensics.Genotype,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as ParentageModel,
    theta as float) {
    checkTheta($theta);
    if ($mother.locus != $child.locus or $child.locus != $father.locus) {
        fail("the mother, child and alleged father genotypes must be at one locus");
    }
    requireNotSilent([$mother, $child, $father]);
    def locus as string init $child.locus;
    def q as float init silentAt($model.silent, $locus);
    def alleles as LocusAlleles init locusAlleles(
        $db,
        $locus,
        [$mother.a, $mother.b, $child.a, $child.b, $father.a, $father.b],
        $q);
    def idx as map of string to int init indexOfAllele($alleles.alleles);
    def mm as MutationMatrix init mutationMatrix($model.maternal, $alleles);
    def fm as MutationMatrix init mutationMatrix($model.paternal, $alleles);
    def children as list of forensics.Genotype init expand($child, $q);
    def x as float init 0.0;
    def y as float init 0.0;
    for (def gm in expand($mother, $q)) {
        for (def gf in expand($father, $q)) {
            def w as float init 1.0;
            if ($q > 0.0) {
                $w = orderings($gm) *
                    orderings($gf) *
                    sequence($alleles, $idx, [$gm.a, $gm.b, $gf.a, $gf.b], $theta);
            }
            def t as list of float init trioTerms(
                $alleles,
                $idx,
                $mm,
                $fm,
                $gm,
                $gf,
                $children,
                $theta);
            $x = $x + $w * $t[0];
            $y = $y + $w * $t[1];
        }
    }
    if ($y <= 0.0) {
        fail("locus '" + $locus + "': the mother cannot have transmitted either of the child's" +
            " alleles under the maternal model - a maternal exclusion, not a paternity index");
    }
    return $x / $y;
}

# commonLoci and coveredLoci mirror kinship.j: a locus every person carries
# but the table lacks is never dropped silently.
func commonLoci(profiles as list of forensics.Profile) {
    def out as list of string init [];
    def first as int init len($profiles) - 1;
    for (def locus in $profiles[$first].loci) {
        def all as bool init true;
        for (def i in 0..$first) {
            if (not maps.has($profiles[$i].loci, $locus)) {
                $all = false;
            }
        }
        if ($all) {
            $out[] = $locus;
        }
    }
    return $out;
}

func coveredLoci(db as forensics.FrequencyDb, profiles as list of forensics.Profile) {
    def loci as list of string init commonLoci($profiles);
    def ids as list of string init [];
    for (def p in $profiles) {
        $ids[] = "'" + $p.id + "'";
    }
    def who as string init "profiles " + strings.join($ids, ", ");
    if (len($loci) == 0) {
        fail($who + " share no locus");
    }
    for (def locus in $loci) {
        if (not forensics.hasLocus($db, $locus)) {
            fail($who + " are all typed at locus '" + $locus + "', which frequency database '" +
                $db.name + "' does not cover; it will not be dropped silently - narrow the" +
                " profiles with forensics.restrict if that is intended");
        }
    }
    return $loci;
}

/**
 * The combined trio paternity index under a parentage model: the product of
 * `paternityIndexAt` over the loci all three profiles carry, in the child's
 * order. A locus all three carry but the table lacks throws.
 * @param mother {Profile} the mother's profile
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {ParentageModel} the mutation models and silent frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the combined paternity index
 * @throws {Error} kind "mutation" on an out-of-range theta, if the trio shares
 *   no locus, shares one the table does not cover, on a maternal exclusion, or
 *   if the CPI exceeds the float range
 */
export func paternityIndex(
    mother as forensics.Profile,
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as ParentageModel,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init coveredLoci($db, [$mother, $allegedFather, $child]);
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            paternityIndexAt(
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

# --- duo --------------------------------------------------------------------

/**
 * The duo paternity index at one locus (mother untyped) under a parentage
 * model: `kinship.duoPaternityIndexAt` with the paternal matrix for the
 * alleged father's transmission and, with a silent allele, single-peak
 * phenotypes summed over `aa` and `a/silent`, each of the father's true
 * genotypes weighted by its probability. The maternal model is not used: the
 * untyped mother's allele is drawn from the population. With no mutation and
 * no silent allele it returns exactly `kinship.duoPaternityIndexAt` under
 * `kinship.NO_MUTATION`.
 * @param db {FrequencyDb} the reference frequencies
 * @param child {Genotype} the child's genotype
 * @param father {Genotype} the alleged father's genotype
 * @param model {ParentageModel} the mutation models and silent frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the paternity index at this locus
 * @throws {Error} kind "mutation" on an out-of-range theta, genotypes at
 *   different loci, a typed `SILENT_ALLELE`, or a child genotype of
 *   probability zero in the population; see `mutationMatrix` for model errors
 */
export func duoPaternityIndexAt(
    db as forensics.FrequencyDb,
    child as forensics.Genotype,
    father as forensics.Genotype,
    model as ParentageModel,
    theta as float) {
    checkTheta($theta);
    if ($child.locus != $father.locus) {
        fail("the child and alleged father genotypes must be at one locus");
    }
    requireNotSilent([$child, $father]);
    def locus as string init $child.locus;
    def q as float init silentAt($model.silent, $locus);
    def alleles as LocusAlleles init locusAlleles(
        $db,
        $locus,
        [$child.a, $child.b, $father.a, $father.b],
        $q);
    def idx as map of string to int init indexOfAllele($alleles.alleles);
    def fm as MutationMatrix init mutationMatrix($model.paternal, $alleles);
    def children as list of forensics.Genotype init expand($child, $q);
    def x as float init 0.0;
    def y as float init 0.0;
    for (def gf in expand($father, $q)) {
        def sf as list of string init [$gf.a, $gf.b];
        def w as float init 1.0;
        if ($q > 0.0) {
            $w = orderings($gf) * sequence($alleles, $idx, $sf, $theta);
        }
        for (def gc in $children) {
            def d1 as float init draw($alleles, $idx, $sf, $gc.a, $theta);
            def d2 as float init draw($alleles, $idx, $sf, $gc.b, $theta);
            def xg as float init transmit($fm, $idx, $gf, $gc.a) * $d2;
            def yg as float init $d1 *
                draw($alleles, $idx, [$gf.a, $gf.b, $gc.a], $gc.b, $theta);
            if ($gc.a != $gc.b) {
                $xg = $xg + transmit($fm, $idx, $gf, $gc.b) * $d1;
                $yg = 2.0 * $yg;
            }
            $x = $x + $w * $xg;
            $y = $y + $w * $yg;
        }
    }
    if ($y <= 0.0) {
        fail("locus '" + $locus + "': the child's genotype has probability zero in the" +
            " reference population, so no ratio is defined");
    }
    return $x / $y;
}

/**
 * The combined duo paternity index under a parentage model: the product of
 * `duoPaternityIndexAt` over the loci both profiles carry.
 * @param child {Profile} the child's profile
 * @param allegedFather {Profile} the alleged father's profile
 * @param db {FrequencyDb} the reference frequencies
 * @param model {ParentageModel} the mutation models and silent frequencies
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the combined paternity index
 * @throws {Error} kind "mutation" on an out-of-range theta, if the two share
 *   no locus, share one the table does not cover, or if the CPI exceeds the
 *   float range
 */
export func duoPaternityIndex(
    child as forensics.Profile,
    allegedFather as forensics.Profile,
    db as forensics.FrequencyDb,
    model as ParentageModel,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init coveredLoci($db, [$allegedFather, $child]);
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

# --- pairwise kinship with a silent allele ----------------------------------

# conditional is kinship.j's conditionalOnAlleleTheta over the locus
# frequencies: P(g | one allele fixed at `allele`, the other drawn).
func conditional(
    alleles as LocusAlleles,
    idx as map of string to int,
    g as forensics.Genotype,
    allele as string,
    sampled as list of string,
    theta as float) {
    if ($g.a == $g.b) {
        if ($g.a != $allele) {
            return 0.0;
        }
        return draw($alleles, $idx, $sampled, $g.a, $theta);
    }
    if ($g.a == $allele) {
        return draw($alleles, $idx, $sampled, $g.b, $theta);
    }
    if ($g.b == $allele) {
        return draw($alleles, $idx, $sampled, $g.a, $theta);
    }
    return 0.0;
}

# trueJoint is kinship.j's jointTheta for two true genotypes over the locus
# frequencies (silent allele included).
func trueJoint(
    alleles as LocusAlleles,
    idx as map of string to int,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as kinship.Relationship,
    theta as float) {
    def sa as list of string init [$ga.a, $ga.b];
    def pa as float init orderings($ga) * sequence($alleles, $idx, $sa, $theta);
    def none as float init orderings($ga) *
        orderings($gb) *
        sequence($alleles, $idx, [$ga.a, $ga.b, $gb.a, $gb.b], $theta);
    def one as float init $pa *
        (0.5 * conditional($alleles, $idx, $gb, $ga.a, $sa, $theta) +
        0.5 * conditional($alleles, $idx, $gb, $ga.b, $sa, $theta));
    def both as float init 0.0;
    if (($ga.a == $gb.a and $ga.b == $gb.b) or ($ga.a == $gb.b and $ga.b == $gb.a)) {
        $both = $pa;
    }
    return $rel.k0 * $none + $rel.k1 * $one + $rel.k2 * $both;
}

# phenotypeJointTheta sums trueJoint over the true genotypes behind two typed
# ones.
func phenotypeJointTheta(
    alleles as LocusAlleles,
    idx as map of string to int,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as kinship.Relationship,
    q as float,
    theta as float) {
    def out as float init 0.0;
    for (def ta in expand($ga, $q)) {
        for (def tb in expand($gb, $q)) {
            $out = $out + trueJoint($alleles, $idx, $ta, $tb, $rel, $theta);
        }
    }
    return $out;
}

/**
 * The kinship likelihood ratio at one locus, `rel` against `alt`, with a
 * silent allele of frequency `silent`: each single-peak phenotype is summed
 * over `aa` and `a/silent` in `kinship.jointGenotypeProbabilityTheta`'s form.
 * At `silent == 0` it returns exactly `kinship.kinshipLRAtTheta`.
 * @param db {FrequencyDb} the reference frequencies
 * @param ga {Genotype} one person's typed genotype
 * @param gb {Genotype} the other's, at the same locus
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param silent {float} the silent allele frequency, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio at this locus
 * @throws {Error} kind "mutation" on an out-of-range silent frequency or
 *   theta, mismatched loci, or if `alt` gives the pair probability zero;
 *   kind "kinship" from `kinship.kinshipLRAtTheta` at `silent == 0`
 */
export func kinshipLRAt(
    db as forensics.FrequencyDb,
    ga as forensics.Genotype,
    gb as forensics.Genotype,
    rel as kinship.Relationship,
    alt as kinship.Relationship,
    silent as float,
    theta as float) {
    checkSilent($silent);
    checkTheta($theta);
    requireNotSilent([$ga, $gb]);
    if ($silent == 0.0) {
        return kinship.kinshipLRAtTheta($db, $ga, $gb, $rel, $alt, $theta);
    }
    if ($ga.locus != $gb.locus) {
        fail("cannot pair a genotype at '" + $ga.locus + "' with one at '" + $gb.locus + "'");
    }
    def alleles as LocusAlleles init locusAlleles(
        $db,
        $ga.locus,
        [$ga.a, $ga.b, $gb.a, $gb.b],
        $silent);
    def idx as map of string to int init indexOfAllele($alleles.alleles);
    def under as float init phenotypeJointTheta($alleles, $idx, $ga, $gb, $alt, $silent, $theta);
    if ($under <= 0.0) {
        fail("locus '" + $ga.locus + "': the alternative hypothesis '" + $alt.name +
            "' gives these genotypes probability zero, so no ratio is defined");
    }
    return phenotypeJointTheta($alleles, $idx, $ga, $gb, $rel, $silent, $theta) / $under;
}

/**
 * The kinship likelihood ratio, `rel` against `alt`, with per-locus silent
 * allele frequencies: the product of `kinshipLRAt` over the loci both
 * profiles carry. With an empty `silent` map it returns exactly
 * `kinship.kinshipLRAgainstTheta`.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @param silent {map of string to float} locus -> silent allele frequency; absent means 0
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "mutation" on an out-of-range theta or silent
 *   frequency, if the profiles share no locus, share one the table does not
 *   cover, if `alt` gives an observed pair probability zero, or if the LR
 *   exceeds the float range (use `log10KinshipLR`); at a locus with silent
 *   frequency 0, kind "kinship" from `kinship.kinshipLRAtTheta`
 */
export func kinshipLR(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as kinship.Relationship,
    alt as kinship.Relationship,
    db as forensics.FrequencyDb,
    silent as map of string to float,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init coveredLoci($db, [$b, $a]);
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = timesLR(
            $out,
            kinshipLRAt(
                $db,
                $a.loci[$locus],
                $b.loci[$locus],
                $rel,
                $alt,
                silentAt($silent, $locus),
                $theta),
            "kinship LR",
            $a.id,
            $b.id,
            "; use log10KinshipLR");
    }
    return $out;
}

/**
 * The base-10 logarithm of `kinshipLR`: the sum over the shared loci of
 * `log10(kinshipLRAt)`. For panels whose LR leaves the float range.
 * @param a {Profile} one person's profile
 * @param b {Profile} the other's
 * @param rel {Relationship} the numerator hypothesis
 * @param alt {Relationship} the denominator hypothesis
 * @param db {FrequencyDb} the reference frequencies
 * @param silent {map of string to float} locus -> silent allele frequency; absent means 0
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "mutation" as `kinshipLR`, or if a locus has ratio 0
 *   (an exclusion, which has no logarithm); kind "kinship" as `kinshipLR`
 */
export func log10KinshipLR(
    a as forensics.Profile,
    b as forensics.Profile,
    rel as kinship.Relationship,
    alt as kinship.Relationship,
    db as forensics.FrequencyDb,
    silent as map of string to float,
    theta as float) {
    checkTheta($theta);
    def loci as list of string init coveredLoci($db, [$b, $a]);
    def out as float init 0.0;
    for (def locus in $loci) {
        def lr as float init kinshipLRAt(
            $db,
            $a.loci[$locus],
            $b.loci[$locus],
            $rel,
            $alt,
            silentAt($silent, $locus),
            $theta);
        if ($lr <= 0.0) {
            fail("profiles '" + $a.id + "' and '" + $b.id + "': the ratio at locus '" + $locus +
                "' is 0 (an exclusion), which has no logarithm");
        }
        $out = $out + math.log10($lr);
    }
    return $out;
}

# --- published mutation counts ----------------------------------------------

func count(
    locus as string,
    maternal as list of int,
    paternal as list of int,
    unassigned as int,
    total as list of int) {
    return MutationCount{
        locus: $locus,
        maternalMutations: $maternal[0],
        maternalMeioses: $maternal[1],
        paternalMutations: $paternal[0],
        paternalMeioses: $paternal[1],
        unassignedMutations: $unassigned,
        totalMutations: $total[0],
        totalMeioses: $total[1]
    };
}

/**
 * Apparent mutations observed at 18 STR loci in paternity testing: the AABB
 * 2003 annual report summary (a compilation of several years, 44
 * laboratories) as published by NIST STRBase, "Mutations Observed at STR
 * Loci", <https://strbase.nist.gov/mutation.htm>, read from the Internet
 * Archive copy of 18 December 2020 on 5 November 2026 (the live page has
 * moved). Counts as printed; locus names as `panels.normalizeLocus` returns
 * them (`VWA` -> `vWA`, `SE33 (ACTBP2)` -> `SE33`); "none reported" is 0.
 * @return {list of MutationCount} one entry per locus, in the published order
 */
export func aabb2003Counts() {
    return [
        count("CSF1PO", [95, 304307], [982, 643118], 410, [1487, 947425]),
        count("FGA", [205, 408230], [2210, 692776], 710, [3125, 1101006]),
        count("TH01", [31, 327172], [41, 452382], 28, [100, 779554]),
        count("TPOX", [18, 400061], [54, 457420], 28, [100, 857481]),
        count("vWA", [184, 564398], [1482, 873547], 814, [2480, 1437945]),
        count("D3S1358", [60, 405452], [713, 558836], 379, [1152, 964288]),
        count("D5S818", [111, 451736], [763, 655603], 385, [1259, 1107339]),
        count("D7S820", [59, 440562], [745, 644743], 285, [1089, 1085305]),
        count("D8S1179", [96, 409869], [779, 489968], 364, [1239, 899837]),
        count("D13S317", [192, 482136], [881, 621146], 485, [1558, 1103282]),
        count("D16S539", [129, 467774], [540, 494465], 372, [1041, 962239]),
        count("D18S51", [186, 296244], [1094, 494098], 466, [1746, 790342]),
        count("D21S11", [464, 435388], [772, 526708], 580, [1816, 962096]),
        count("Penta D", [12, 18701], [21, 22501], 24, [57, 41202]),
        count("Penta E", [29, 44311], [75, 55719], 59, [163, 100030]),
        count("D2S1338", [15, 72830], [157, 152310], 90, [262, 225140]),
        count("D19S433", [38, 70001], [78, 103489], 71, [187, 173490]),
        count("SE33", [0, 330], [330, 51610], 0, [330, 51940])
    ];
}

/**
 * Per-locus mutation rates from `aabb2003Counts`: `"total"` is all mutations
 * over all meioses (STRBase's "mutation rate" column); `"paternal"` and
 * `"maternal"` are the mutations assigned to that parent over that parent's
 * meioses, which leave out the unassigned ones and so understate the rate.
 * @param which {string} `"total"`, `"paternal"` or `"maternal"`
 * @return {map of string to float} locus -> rate per meiosis
 * @throws {Error} kind "mutation" on any other `which`
 */
export func aabb2003Rates(which as string) {
    if ($which != "total" and $which != "paternal" and $which != "maternal") {
        fail("rates are 'total', 'paternal' or 'maternal', not '" + $which + "'");
    }
    def out as map of string to float init {};
    for (def c in aabb2003Counts()) {
        if ($which == "total") {
            $out[$c.locus] = convert.toFloat($c.totalMutations) / convert.toFloat($c.totalMeioses);
        } elseif ($which == "paternal") {
            $out[$c.locus] = convert.toFloat($c.paternalMutations) /
                convert.toFloat($c.paternalMeioses);
        } else {
            $out[$c.locus] = convert.toFloat($c.maternalMutations) /
                convert.toFloat($c.maternalMeioses);
        }
    }
    return $out;
}

/**
 * A copy of a `kinship.MutationModel` with several loci's rates overridden,
 * for `kinship.paternityIndexWith` and the other `kinship.j` functions.
 * @param model {MutationModel} the kinship.j model to copy
 * @param rates {map of string to float} locus -> mutation rate
 * @return {MutationModel} the updated copy
 * @throws {Error} kind "kinship" on an out-of-range rate
 */
export func kinshipModelWithRates(model as kinship.MutationModel, rates as map of string to float) {
    def out as kinship.MutationModel init $model;
    for (def locus in $rates) {
        $out = kinship.withLocusRate($out, $locus, $rates[$locus]);
    }
    return $out;
}
