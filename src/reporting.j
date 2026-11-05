# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Reporting a single-source match: allowances for the finite size of the
 * frequency database, and the verbal equivalent of a likelihood ratio.
 *
 * Every allele frequency in a `FrequencyDb` is an estimate from `2N` sampled
 * alleles; this module reads the locus's counts back as `x = p · 2N` (not
 * rounded: a published frequency rarely gives a whole count). Two published
 * allowances for that sampling error:
 *
 * - the size-bias correction (Balding & Nichols 1994, as stated in the UK
 *   Forensic Science Regulator's FSR-GUI-0012, 2026, section 11.3): the suspect's and
 *   the offender's alleles are added to the database before the frequency is
 *   taken, `(x + 2) / (2N + 4)` for each allele of a heterozygote and
 *   `(x + 4) / (2N + 4)` for a homozygote;
 * - a posterior quantile of the match probability (the highest-posterior-
 *   density method of Curran, Buckleton, Triggs & Weir 2002): each locus's
 *   allele frequencies get a Dirichlet posterior, `Dir(x_i + α)` with a
 *   uniform (`α = 1`) or `1/k` prior (Curran & Buckleton 2011), and the
 *   random match probability's posterior is simulated. The upper `level`
 *   quantile is a one-sided bound on the RMP, its reciprocal a lower bound on
 *   the likelihood ratio.
 *
 * Both evaluate the NRC II 4.2 genotype frequency (the formula of
 * `forensics.genotypeFrequency`) at the adjusted frequencies, so theta applies
 * as elsewhere in the deck. Neither
 * applies the 5/(2N) floor: the added alleles and the prior are the allowance
 * for rare alleles, and an unlisted allele counts as `x = 0`. Every locus needs
 * its sample size (`n > 0`).
 *
 * The simulation uses its own Park-Miller generator (the one `popstats.j`
 * uses), normal deviates by the polar method (Marsaglia & Bray 1964) and gamma
 * deviates by Marsaglia & Tsang (2000), so a seed fixes the answer and
 * `math.rand` is not touched.
 *
 * `verbalEquivalent` places a likelihood ratio on the verbal scale of the ENFSI
 * Guideline for Evaluative Reporting in Forensic Science (approved version
 * 3.0, 2015, guidance note 4, p. 17). The wording is the guideline's; it
 * gives the scale as an illustration and calls its ranges indicative.
 * @module reporting
 * @example
 * import "@deckyard/forensicgenetics/" as forensics;
 * import "@deckyard/forensicgenetics/reporting.j" as reporting;
 *
 * def lr as float init reporting.sizeBiasLikelihoodRatio($db, $suspect, 0.01);
 * def bound as float init reporting.likelihoodRatioLowerBound(
 *     $db, $suspect, 0.01, reporting.PRIOR_ONE_OVER_K, 0.99, reporting.DEFAULT_DRAWS, 1);
 * def words as reporting.VerbalEquivalent init reporting.verbalEquivalent($bound);
 * def say as string init $words.statement;
 */
use math;
use maps;
use lists;
use convert;

import "./forensicgenetics.j" as forensics;

/**
 * Dirichlet prior `α = 1` on every allele (the uniform prior).
 */
export def const PRIOR_UNIFORM as int init 1;

/**
 * Dirichlet prior `α = 1/k` on every allele, `k` the alleles at the locus
 * (listed, plus any of the profile's the table lacks).
 */
export def const PRIOR_ONE_OVER_K as int init 2;

/**
 * A simulation size for the posterior quantiles. The Monte Carlo error of a
 * quantile is about `sqrt(q (1 - q) / draws)` on the probability scale.
 */
export def const DEFAULT_DRAWS as int init 10000;

/**
 * The largest simulation size the posterior functions accept.
 */
export def const MAX_DRAWS as int init 1000000;

/**
 * `VerbalEquivalent.supports` when the LR exceeds 1.
 */
export def const SUPPORTS_FIRST as string init "first";

/**
 * `VerbalEquivalent.supports` when the LR is below 1.
 */
export def const SUPPORTS_ALTERNATIVE as string init "alternative";

/**
 * `VerbalEquivalent.supports` when the LR is exactly 1.
 */
export def const SUPPORTS_NEITHER as string init "neither";

/**
 * A likelihood ratio's place on the ENFSI (2015) verbal scale.
 *
 * Bounds are on the LR scale. Bands supporting the first proposition are
 * `(lower, upper]`; their mirror images, supporting the alternative, are
 * `[lower, upper)`; LR = 1 is `[1, 1]`. The top band has no upper bound
 * (`hasUpper` false, `upper` 0.0) and its mirror image runs down to 0
 * (exclusive).
 * @field label {string} the guideline's term: "weak support" to "extremely
 *   strong support", or "no support" at LR = 1
 * @field supports {string} `SUPPORTS_FIRST`, `SUPPORTS_ALTERNATIVE` or
 *   `SUPPORTS_NEITHER`
 * @field statement {string} the guideline's first phrasing for the band,
 *   with the propositions exchanged when the alternative is supported
 * @field lower {float} the band's lower LR bound
 * @field lowerInclusive {bool} whether `lower` belongs to the band
 * @field upper {float} the band's upper LR bound; 0.0 when `hasUpper` is false
 * @field upperInclusive {bool} whether `upper` belongs to the band
 * @field hasUpper {bool} false only for the open-ended top band
 */
export def struct VerbalEquivalent {
    label as string,
    supports as string,
    statement as string,
    lower as float,
    lowerInclusive as bool,
    upper as float,
    upperInclusive as bool,
    hasUpper as bool
};

# Park-Miller "minimal standard" generator with the Park, Miller & Stockmeyer
# (1993) multiplier, as in popstats.j (private there): x <- 48271 x mod
# (2^31 - 1). Products stay below 2^47, so int64 never overflows.
def const RAND_MODULUS as int init 2147483647;
def const RAND_MULTIPLIER as int init 48271;

# A remainder count (2N minus the listed counts) below this share of 2N is
# rounding in the table, not unlisted alleles.
def const REST_TOLERANCE as float init 0.000000001;

# The generator's state, the deviate just drawn, and the polar method's
# unused second normal. Every sampler takes one and returns the next, since
# arguments are copies.
def struct Stream {
    state as int,
    value as float,
    spare as float,
    hasSpare as bool
};

# One locus of the simulation: its Dirichlet parameters, aggregated to the
# profile's alleles and one class for every other allele (the Dirichlet
# aggregation property makes that exact).
def struct LocusPrior {
    homozygous as bool,
    alphaA as float,
    alphaB as float,
    alphaRest as float
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "reporting", message: $msg, file: "", line: 0, col: 0};
}

func checkTheta(theta as float) {
    if ($theta < 0.0 or $theta >= 1.0) {
        fail("theta must be in [0, 1)");
    }
}

func checkProfile(p as forensics.Profile) {
    if (len($p.loci) == 0) {
        fail("profile '" + $p.id + "' has no genotypes");
    }
}

# locusCounts fetches a locus for a count-based method: it must exist and
# carry its sample size, since counts are p * 2N.
func locusCounts(fdb as forensics.FrequencyDb, locus as string) {
    if (not maps.has($fdb.loci, $locus)) {
        fail("frequency database '" + $fdb.name + "' has no locus '" + $locus + "'");
    }
    def lf as forensics.LocusFrequencies init $fdb.loci[$locus];
    if ($lf.n <= 0) {
        fail("locus '" + $locus + "' has no sample size; allele counts need N");
    }
    return $lf;
}

# count is an allele's database count, p * 2N, or 0 if unlisted.
func count(lf as forensics.LocusFrequencies, allele as string) {
    if (maps.has($lf.alleles, $allele)) {
        return $lf.alleles[$allele] * 2.0 * convert.toFloat($lf.n);
    }
    return 0.0;
}

# nrcGenotype is NRC II 4.2 (the Balding-Nichols conditional genotype
# probability) at given frequencies. No cap is needed: a homozygote's
# frequency is at most 1 (x <= 2N), and a heterozygote's pair sums to at most
# about 1.05 (MAX_FREQUENCY_SUM), which keeps 2 pa pb near .55.
func nrcGenotype(pa as float, pb as float, homozygous as bool, theta as float) {
    def t as float init $theta;
    def den as float init (1.0 + $t) * (1.0 + 2.0 * $t);
    if ($homozygous) {
        return (2.0 * $t + (1.0 - $t) * $pa) * (3.0 * $t + (1.0 - $t) * $pa) / $den;
    }
    return 2.0 * ($t + (1.0 - $t) * $pa) * ($t + (1.0 - $t) * $pb) / $den;
}

# --- size-bias correction ---------------------------------------------------

/**
 * One allele's size-bias corrected frequency for a single-source match.
 *
 * The suspect's and the offender's genotype `g` are added to the database:
 * `(x + c) / (2N + 4)`, where `x = p · 2N` is the allele's count (0 if
 * unlisted) and `c` the copies the two people carry (2 per heterozygote
 * allele, 4 for a homozygote, 0 for an allele not in `g`). No 5/(2N) floor.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param g {forensics.Genotype} the matching genotype
 * @param allele {string} the allele to score
 * @return {float} the corrected frequency
 * @throws {Error} kind "reporting" on a locus the database lacks or one with
 *   `n == 0`
 */
export func sizeBiasAlleleFrequency(
    fdb as forensics.FrequencyDb,
    g as forensics.Genotype,
    allele as string) {
    def lf as forensics.LocusFrequencies init locusCounts($fdb, $g.locus);
    def added as float init 0.0;
    if ($allele == $g.a) {
        $added = $added + 2.0;
    }
    if ($allele == $g.b) {
        $added = $added + 2.0;
    }
    def twoN as float init 2.0 * convert.toFloat($lf.n);
    return (count($lf, $allele) + $added) / ($twoN + 4.0);
}

/**
 * The NRC II 4.2 genotype frequency at size-bias corrected allele
 * frequencies (`sizeBiasAlleleFrequency`).
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param g {forensics.Genotype} the matching genotype
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the corrected genotype frequency
 * @throws {Error} kind "reporting" on an out-of-range theta, a locus the
 *   database lacks or one with `n == 0`
 */
export func sizeBiasGenotypeFrequency(
    fdb as forensics.FrequencyDb,
    g as forensics.Genotype,
    theta as float) {
    checkTheta($theta);
    def pa as float init sizeBiasAlleleFrequency($fdb, $g, $g.a);
    def homozygous as bool init forensics.isHomozygous($g);
    def pb as float init $pa;
    if (not $homozygous) {
        $pb = sizeBiasAlleleFrequency($fdb, $g, $g.b);
    }
    return nrcGenotype($pa, $pb, $homozygous, $theta);
}

/**
 * The random match probability with the size-bias correction: the product of
 * `sizeBiasGenotypeFrequency` over the profile's loci. Every locus must be in
 * the database with its sample size.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param p {forensics.Profile} the matching profile
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the corrected random match probability
 * @throws {Error} kind "reporting" on an empty profile, an out-of-range theta,
 *   a locus the database lacks or one with `n == 0`
 */
export func sizeBiasRandomMatchProbability(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float) {
    checkProfile($p);
    def out as float init 1.0;
    for (def locus in $p.loci) {
        $out = $out * sizeBiasGenotypeFrequency($fdb, $p.loci[$locus], $theta);
    }
    return $out;
}

/**
 * The base-10 logarithm of `sizeBiasRandomMatchProbability`, summed per locus
 * so a long profile does not underflow.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param p {forensics.Profile} the matching profile
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the corrected random match probability
 * @throws {Error} kind "reporting" as `sizeBiasRandomMatchProbability`
 */
export func log10SizeBiasRandomMatchProbability(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float) {
    checkProfile($p);
    def out as float init 0.0;
    for (def locus in $p.loci) {
        $out = $out + math.log10(sizeBiasGenotypeFrequency($fdb, $p.loci[$locus], $theta));
    }
    return $out;
}

/**
 * The single-source likelihood ratio with the size-bias correction,
 * `1 / sizeBiasRandomMatchProbability`.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param p {forensics.Profile} the matching profile
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the corrected likelihood ratio
 * @throws {Error} kind "reporting" as `sizeBiasRandomMatchProbability`, or if
 *   the product underflows to zero
 */
export func sizeBiasLikelihoodRatio(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float) {
    def rmp as float init sizeBiasRandomMatchProbability($fdb, $p, $theta);
    if ($rmp <= 0.0) {
        fail("profile '" + $p.id + "' has a zero random match probability");
    }
    return 1.0 / $rmp;
}

# --- random deviates --------------------------------------------------------

# stream starts a generator at a seed in [1, 2^31 - 2].
func stream(seed as int) {
    return Stream{state: $seed, value: 0.0, spare: 0.0, hasSpare: false};
}

# nextUniform is one Park-Miller step; the value state / (2^31 - 1) lies in
# (0, 1), never 0 or 1, so its logarithm is always finite.
func nextUniform(s as Stream) {
    def st as Stream init $s;
    $st.state = ($s.state * RAND_MULTIPLIER) % RAND_MODULUS;
    $st.value = convert.toFloat($st.state) / 2147483647.0;
    return $st;
}

# nextNormal is a standard normal deviate by the polar method (Marsaglia &
# Bray 1964): (v1, v2) uniform in the unit disc, s = v1^2 + v2^2, then
# v1 m and v2 m with m = sqrt(-2 ln s / s), two independent deviates. The
# second is kept for the next call.
func nextNormal(s as Stream) {
    def st as Stream init $s;
    if ($st.hasSpare) {
        $st.hasSpare = false;
        $st.value = $st.spare;
        return $st;
    }
    while (true) {
        $st = nextUniform($st);
        def v1 as float init 2.0 * $st.value - 1.0;
        $st = nextUniform($st);
        def v2 as float init 2.0 * $st.value - 1.0;
        def r as float init $v1 * $v1 + $v2 * $v2;
        if ($r < 1.0 and $r > 0.0) {
            def m as float init math.sqrt(-2.0 * math.ln($r) / $r);
            $st.value = $v1 * $m;
            $st.spare = $v2 * $m;
            $st.hasSpare = true;
            return $st;
        }
    }
    return $st;
}

# nextLnGamma draws ln X for X ~ Gamma(alpha, 1) (Marsaglia & Tsang 2000). For
# alpha >= 1: d = alpha - 1/3, c = 1/sqrt(9d); x normal, v = (1 + c x)^3 > 0,
# u uniform; accept d v if u < 1 - 0.0331 x^4 (the squeeze) or
# ln u < x^2/2 + d (1 - v + ln v). For alpha < 1, X = Gamma(alpha + 1) *
# U^(1/alpha). The logarithm is returned because U^(1/alpha) underflows for
# small alpha.
func nextLnGamma(s as Stream, alpha as float) {
    def st as Stream init $s;
    def shift as float init 0.0;
    def a as float init $alpha;
    if ($a < 1.0) {
        $st = nextUniform($st);
        $shift = math.ln($st.value) / $a;
        $a = $a + 1.0;
    }
    def d as float init $a - 1.0 / 3.0;
    def c as float init 1.0 / math.sqrt(9.0 * $d);
    while (true) {
        $st = nextNormal($st);
        def x as float init $st.value;
        def v as float init 1.0 + $c * $x;
        if ($v > 0.0) {
            $v = $v * $v * $v;
            $st = nextUniform($st);
            def u as float init $st.value;
            def x2 as float init $x * $x;
            if ($u < 1.0 - 0.0331 * $x2 * $x2 or
                math.ln($u) < 0.5 * $x2 + $d * (1.0 - $v + math.ln($v))) {
                $st.value = math.ln($d * $v) + $shift;
                return $st;
            }
        }
    }
    return $st;
}

# --- posterior quantiles ----------------------------------------------------

func checkPrior(prior as int) {
    if ($prior != PRIOR_UNIFORM and $prior != PRIOR_ONE_OVER_K) {
        fail("prior must be PRIOR_UNIFORM or PRIOR_ONE_OVER_K, got {$prior}");
    }
}

func checkSimulation(draws as int, seed as int) {
    if ($draws < 1 or $draws > MAX_DRAWS) {
        fail("draws must be between 1 and 1000000, got {$draws}");
    }
    if ($seed < 1 or $seed >= RAND_MODULUS) {
        fail("seed must be between 1 and 2147483646, got {$seed}");
    }
}

# locusPrior builds one locus's aggregated Dirichlet posterior: the profile's
# allele(s) with x + alpha each, and every other allele in one class with
# their counts and priors summed, plus the count the listed frequencies leave
# short of 2N (unlisted alleles, no prior).
func locusPrior(fdb as forensics.FrequencyDb, g as forensics.Genotype, prior as int) {
    def lf as forensics.LocusFrequencies init locusCounts($fdb, $g.locus);
    def homozygous as bool init forensics.isHomozygous($g);
    def twoN as float init 2.0 * convert.toFloat($lf.n);
    def k as int init len($lf.alleles);
    def listed as float init 0.0;
    for (def a in $lf.alleles) {
        $listed = $listed + $lf.alleles[$a] * $twoN;
    }
    if (not maps.has($lf.alleles, $g.a)) {
        $k = $k + 1;
    }
    if (not $homozygous and not maps.has($lf.alleles, $g.b)) {
        $k = $k + 1;
    }
    def w as float init 1.0;
    if ($prior == PRIOR_ONE_OVER_K) {
        $w = 1.0 / convert.toFloat($k);
    }
    def own as int init 2;
    def xa as float init count($lf, $g.a);
    def xb as float init 0.0;
    if ($homozygous) {
        $own = 1;
    } else {
        $xb = count($lf, $g.b);
    }
    def rest as float init $listed - $xa - $xb + convert.toFloat($k - $own) * $w;
    if ($twoN - $listed > REST_TOLERANCE * $twoN) {
        $rest = $rest + ($twoN - $listed);
    }
    if ($rest < REST_TOLERANCE) {
        $rest = 0.0;
    }
    return LocusPrior{
        homozygous: $homozygous,
        alphaA: $xa + $w,
        alphaB: $xb + $w,
        alphaRest: $rest
    };
}

# simulate draws `draws` values of log10 RMP from the posterior, sorted
# ascending. Per draw and locus the profile's allele frequencies are gamma
# deviates normalised by their sum (the Dirichlet construction), all kept in
# logarithms so a tiny frequency cannot underflow to a zero RMP.
func simulate(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float,
    prior as int,
    draws as int,
    seed as int) {
    checkProfile($p);
    checkTheta($theta);
    checkPrior($prior);
    checkSimulation($draws, $seed);
    def priors as list of LocusPrior init [];
    for (def locus in $p.loci) {
        $priors[] = locusPrior($fdb, $p.loci[$locus], $prior);
    }
    def nl as int init len($priors);
    def ln10 as float init math.ln(10.0);
    def log2 as float init math.log10(2.0);
    def st as Stream init stream($seed);
    def out as list of float init [];
    for (def d in 0..$draws) {
        def total as float init 0.0;
        for (def i in 0..$nl) {
            $st = nextLnGamma($st, $priors[$i].alphaA);
            def la as float init $st.value;
            def lb as float init $la;
            if (not $priors[$i].homozygous) {
                $st = nextLnGamma($st, $priors[$i].alphaB);
                $lb = $st.value;
            }
            def m as float init $la;
            def lr as float init 0.0;
            def hasRest as bool init $priors[$i].alphaRest > 0.0;
            if ($hasRest) {
                $st = nextLnGamma($st, $priors[$i].alphaRest);
                $lr = $st.value;
                $m = math.max($m, $lr);
            }
            if (not $priors[$i].homozygous) {
                $m = math.max($m, $lb);
            }
            def sum as float init math.exp($la - $m);
            if (not $priors[$i].homozygous) {
                $sum = $sum + math.exp($lb - $m);
            }
            if ($hasRest) {
                $sum = $sum + math.exp($lr - $m);
            }
            def ls as float init $m + math.ln($sum);
            $total = $total +
                locusLog10($la - $ls, $lb - $ls, $priors[$i].homozygous, $theta, $ln10, $log2);
        }
        $out[] = $total;
    }
    return lists.sort($out);
}

# locusLog10 is log10 of one locus's genotype frequency from the natural logs
# of its allele frequencies: exact in logarithms at theta 0, NRC II 4.2
# otherwise (where theta keeps every factor away from zero).
func locusLog10(
    lnA as float,
    lnB as float,
    homozygous as bool,
    theta as float,
    ln10 as float,
    log2 as float) {
    if ($theta == 0.0) {
        if ($homozygous) {
            return 2.0 * $lnA / $ln10;
        }
        return $log2 + ($lnA + $lnB) / $ln10;
    }
    return math.log10(nrcGenotype(math.exp($lnA), math.exp($lnB), $homozygous, $theta));
}

# quantileOf is the nearest-rank quantile of sorted values: the ceil(q M)-th
# smallest of M.
func quantileOf(sorted as list of float, q as float) {
    def m as int init len($sorted);
    def rank as int init math.ceil($q * convert.toFloat($m));
    if ($rank < 1) {
        $rank = 1;
    }
    if ($rank > $m) {
        $rank = $m;
    }
    return $sorted[$rank - 1];
}

/**
 * Posterior quantiles of the log10 random match probability under a
 * Dirichlet posterior on each locus's allele frequencies (Curran, Buckleton,
 * Triggs & Weir 2002), by seeded Monte Carlo.
 *
 * Locus by locus the frequencies are `Dir(x_i + α)`, `x_i = p_i · 2N`, with
 * `α` from `prior`; the count `2N - Σ x_i` a partial table leaves is one more
 * class with no prior. Loci are independent. Each draw evaluates NRC II 4.2 at
 * the drawn frequencies and multiplies over loci. The q-quantile is the
 * nearest-rank one, the `ceil(q · draws)`-th smallest draw.
 * @param fdb {forensics.FrequencyDb} the reference frequencies, `n > 0` at
 *   every profile locus
 * @param p {forensics.Profile} the matching profile
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @param prior {int} `PRIOR_UNIFORM` or `PRIOR_ONE_OVER_K`
 * @param quantiles {list of float} the quantiles wanted, each in (0, 1)
 * @param draws {int} the simulation size, 1 to `MAX_DRAWS`
 * @param seed {int} the generator seed, 1 to 2147483646
 * @return {list of float} log10 RMP at each quantile, in `quantiles` order
 * @throws {Error} kind "reporting" on an empty profile or quantile list, a
 *   quantile outside (0, 1), an out-of-range theta, prior, draw count or seed,
 *   a locus the database lacks or one with `n == 0`
 */
export func log10RmpQuantiles(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float,
    prior as int,
    quantiles as list of float,
    draws as int,
    seed as int) {
    if (len($quantiles) == 0) {
        fail("no quantiles requested");
    }
    for (def q in $quantiles) {
        if ($q <= 0.0 or $q >= 1.0) {
            fail("a quantile must be in (0, 1), got {$q}");
        }
    }
    def sorted as list of float init simulate($fdb, $p, $theta, $prior, $draws, $seed);
    def out as list of float init [];
    for (def q in $quantiles) {
        $out[] = quantileOf($sorted, $q);
    }
    return $out;
}

/**
 * A one-sided lower bound on the single-source likelihood ratio: `1 / RMP_level`,
 * `RMP_level` the `level` posterior quantile of the RMP from
 * `log10RmpQuantiles`. A `level` of 0.99 gives the bound the LR exceeds with
 * posterior probability 0.99.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param p {forensics.Profile} the matching profile
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @param prior {int} `PRIOR_UNIFORM` or `PRIOR_ONE_OVER_K`
 * @param level {float} the posterior probability, in (0, 1)
 * @param draws {int} the simulation size, 1 to `MAX_DRAWS`
 * @param seed {int} the generator seed, 1 to 2147483646
 * @return {float} the likelihood-ratio bound
 * @throws {Error} kind "reporting" as `log10RmpQuantiles`, or if the bound
 *   exceeds 1e300 (use `log10RmpQuantiles` then)
 */
export func likelihoodRatioLowerBound(
    fdb as forensics.FrequencyDb,
    p as forensics.Profile,
    theta as float,
    prior as int,
    level as float,
    draws as int,
    seed as int) {
    def q as list of float init log10RmpQuantiles(
        $fdb,
        $p,
        $theta,
        $prior,
        [$level],
        $draws,
        $seed);
    def exponent as float init -$q[0];
    if ($exponent > 300.0) {
        fail("the bound is 10^{$exponent}, past float range; use log10RmpQuantiles");
    }
    return math.pow(10.0, $exponent);
}

# --- verbal equivalents -----------------------------------------------------

# The guideline's six bands above LR = 1, by log10 upper bound (the last has
# none). Labels and the first phrasing are from the table on p. 17.
func bandLabel(band as int) {
    def labels as list of string init [
        "weak support",
        "moderate support",
        "moderately strong support",
        "strong support",
        "very strong support",
        "extremely strong support"
    ];
    return $labels[$band];
}

func bandStatement(band as int, first as bool) {
    def x as string init "the first proposition";
    def y as string init "the alternative";
    if (not $first) {
        $x = "the alternative";
        $y = "the first proposition";
    }
    if ($band == 0) {
        return "The forensic findings provide weak support for " + $x + " relative to " + $y + ".";
    }
    return "The forensic findings provide " + bandLabel($band) + " for " + $x + " rather than " +
        $y + ".";
}

# band picks the band of a magnitude given as log10 |log10 LR| > 0: (0, 1],
# (1, 2], (2, 3], (3, 4], (4, 6], (6, inf).
func bandOfLog(m as float) {
    def tops as list of float init [1.0, 2.0, 3.0, 4.0, 6.0];
    for (def i in 0..len($tops)) {
        if ($m <= $tops[$i]) {
            return $i;
        }
    }
    return 5;
}

# bandOf is bandOfLog decided on the LR itself, so a boundary LR such as 1000
# or 0.001 lands where the scale puts it and not where log10 rounding would.
func bandOf(lr as float) {
    if ($lr > 1.0) {
        def tops as list of float init [10.0, 100.0, 1000.0, 10000.0, 1000000.0];
        for (def i in 0..len($tops)) {
            if ($lr <= $tops[$i]) {
                return $i;
            }
        }
        return 5;
    }
    def bottoms as list of float init [0.1, 0.01, 0.001, 0.0001, 0.000001];
    for (def i in 0..len($bottoms)) {
        if ($lr >= $bottoms[$i]) {
            return $i;
        }
    }
    return 5;
}

func neutral() {
    return VerbalEquivalent{
        label: "no support",
        supports: SUPPORTS_NEITHER,
        statement: "The forensic findings do not support one proposition over the other.",
        lower: 1.0,
        lowerInclusive: true,
        upper: 1.0,
        upperInclusive: true,
        hasUpper: true
    };
}

# describe builds the band's record; bounds come from the exponent table so
# they are the scale's round numbers.
func describe(band as int, first as bool) {
    def exps as list of float init [0.0, 1.0, 2.0, 3.0, 4.0, 6.0];
    def lo as float init math.pow(10.0, $exps[$band]);
    def hasHi as bool init $band < 5;
    def hi as float init 0.0;
    if ($hasHi) {
        $hi = math.pow(10.0, $exps[$band + 1]);
    }
    if ($first) {
        return VerbalEquivalent{
            label: bandLabel($band),
            supports: SUPPORTS_FIRST,
            statement: bandStatement($band, true),
            lower: $lo,
            lowerInclusive: false,
            upper: $hi,
            upperInclusive: true,
            hasUpper: $hasHi
        };
    }
    # Mirror image on the LR scale: [1/hi, 1/lo), down to 0 for the top band.
    def mirrorLo as float init 0.0;
    if ($hasHi) {
        $mirrorLo = math.pow(10.0, -$exps[$band + 1]);
    }
    return VerbalEquivalent{
        label: bandLabel($band),
        supports: SUPPORTS_ALTERNATIVE,
        statement: bandStatement($band, false),
        lower: $mirrorLo,
        lowerInclusive: $hasHi,
        upper: math.pow(10.0, -$exps[$band]),
        upperInclusive: false,
        hasUpper: true
    };
}

/**
 * A likelihood ratio's verbal equivalent on the ENFSI Guideline for
 * Evaluative Reporting in Forensic Science (2015) scale.
 *
 * Above 1: weak (to 10), moderate (to 100), moderately strong (to 1000),
 * strong (to 10,000), very strong (to 1,000,000) and extremely strong support
 * for the first proposition. Below 1 the guideline reads `1/LR` on the same
 * scale as support for the alternative. Exactly 1 supports neither.
 *
 * The guideline's table writes ranges ("2 - 10", "10 - 100", ...) without
 * saying which end is inclusive, and calls them indicative. Here each band
 * includes its upper bound - `(10, 100]`, mirrored as `[0.01, 0.1)` - as the
 * guideline's speaker-recognition example scale (p. 64) writes its own bands.
 * The table's first band starts at 2; LRs in (1, 2) are put in it.
 * @param lr {float} the likelihood ratio, greater than 0
 * @return {VerbalEquivalent} the band, its wording and its bounds
 * @throws {Error} kind "reporting" if `lr` is not positive
 */
export func verbalEquivalent(lr as float) {
    if ($lr <= 0.0) {
        fail("a likelihood ratio must be positive for the verbal scale, got {$lr}");
    }
    if ($lr == 1.0) {
        return neutral();
    }
    return describe(bandOf($lr), $lr > 1.0);
}

/**
 * `verbalEquivalent` for a likelihood ratio given as its base-10 logarithm,
 * for LRs past float range (`mixture.log10MixtureLR`). Band edges are the
 * integers 1, 2, 3, 4 and 6 (and their negatives), compared exactly.
 * @param log10Lr {float} log10 of the likelihood ratio
 * @return {VerbalEquivalent} the band, its wording and its bounds
 */
export func verbalEquivalentLog10(log10Lr as float) {
    if ($log10Lr == 0.0) {
        return neutral();
    }
    if ($log10Lr > 0.0) {
        return describe(bandOfLog($log10Lr), true);
    }
    # The mirrored band [-1, 0) is (0, 1] in magnitude: the same closed edge.
    return describe(bandOfLog(-$log10Lr), false);
}
