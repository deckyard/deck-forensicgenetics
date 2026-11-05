# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * The semi-continuous (qualitative) mixture likelihood ratio: the drop-out /
 * drop-in model of Curran et al. (2005) and Gill et al. (2007), as implemented
 * in LRmix (Haned et al. 2012; Gill & Haned 2013).
 *
 * The evidence at a locus is the set of alleles a core `Mixture` records -
 * presence only, no peak heights. A `Hypothesis` names its known (conditioned)
 * contributors with a drop-out probability each, and a number of unknown
 * contributors sharing one drop-out probability. Per locus,
 *
 * ```
 *   P(E | H) = Σ over genotypes G of the unknowns of
 *              P(E | knowns, G; d, c) · P(G | typed profiles; θ)
 *
 *   LR = Π over loci of P(E | Hp) / P(E | Hd)
 * ```
 *
 * Drop-out is independent per allele copy: `d` for a heterozygote's allele,
 * `d²` for a homozygote's (the Gill et al. 2007 approximation, not the
 * two-parameter form of Balding & Buckleton 2009); an allele several
 * contributors carry is lost only if every copy drops out. An observed allele
 * no contributor carries is a drop-in, `c · p_a` each; no drop-in is `1 - c`.
 *
 * Unknown genotypes follow the Balding-Nichols sampling formula (Balding &
 * Nichols 1994; for mixtures Curran et al. 1999), conditioned on the known
 * contributors of both hypotheses (`conditioningProfiles`); `θ = 0` is
 * Hardy-Weinberg. One known suspect against one unknown with `d = 0`, `c = 0`
 * gives `forensicgenetics.matchLikelihoodRatio` at the same θ.
 *
 * Genotypes range over the alleles in the stain or a typed profile plus one
 * lumped class holding the remaining frequency (exact: those alleles are never
 * observed). Unknowns share `d`, so combinations are enumerated as multisets
 * weighted by their orderings. `MAX_UNKNOWNS` and `MAX_LOCUS_ALLELES` bound
 * the `(k(k+1)/2)^u` cost.
 * @module mixture
 * @example
 * import "@deckyard/forensicgenetics/" as forensics;
 * import "@deckyard/forensicgenetics/mixture.j" as mixture;
 *
 * def hp as mixture.Hypothesis init mixture.hypothesis("Hp", [$suspect], [0.1], 1, 0.1);
 * def hd as mixture.Hypothesis init mixture.hypothesis("Hd", [], [], 2, 0.1);
 * def lr as float init mixture.mixtureLR($db, $stain, $hp, $hd, 0.05, 0.01);
 */
use math;
use maps;
use lists;
use convert;

import "./forensicgenetics.j" as forensics;

/**
 * The most unknown contributors a hypothesis may name.
 *
 * Three unknowns at `MAX_LOCUS_ALLELES` alleles take about 4 s per locus.
 */
export def const MAX_UNKNOWNS as int init 3;

/**
 * The most distinct alleles one locus may enumerate: those in the stain plus
 * those of typed profiles, excluding the lumped class.
 */
export def const MAX_LOCUS_ALLELES as int init 10;

/**
 * One hypothesis about who contributed to a mixture.
 *
 * Build one with `hypothesis`, which validates it.
 * @field name {string} a label, used in error messages
 * @field known {list of forensics.Profile} the known contributors
 * @field knownDropout {list of float} the drop-out probability of each known
 *   contributor, in `known` order
 * @field unknowns {int} the number of unknown contributors, 0 to `MAX_UNKNOWNS`
 * @field unknownDropout {float} the drop-out probability shared by every
 *   unknown contributor
 */
export def struct Hypothesis {
    name as string,
    known as list of forensics.Profile,
    knownDropout as list of float,
    unknowns as int,
    unknownDropout as float
};

# Frame is one locus prepared for one hypothesis: the allele classes the
# unknowns' genotypes range over, and everything about them that does not
# depend on those genotypes.
def struct Frame {
    labels as list of string,
    freq as list of float,
    observed as list of bool,
    typedCount as list of float,
    typedTotal as float,
    knownDrop as list of float,
    knownCarry as list of bool,
    gA as list of int,
    gB as list of int,
    hw as list of float,
    zeros as list of int,
    unknownPow as list of float
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "mixture", message: $msg, file: "", line: 0, col: 0};
}

# --- hypotheses -------------------------------------------------------------

# checkProbability rejects a probability outside [0, 1].
func checkProbability(what as string, p as float) {
    if ($p < 0.0 or $p > 1.0) {
        fail($what + " must be in [0, 1], got {$p}");
    }
}

# checkHypothesis validates a hypothesis, however it was built: a struct
# literal skips `hypothesis`, so the likelihood functions check again.
func checkHypothesis(h as Hypothesis) {
    def nk as int init len($h.known);
    def nd as int init len($h.knownDropout);
    if ($nk != $nd) {
        fail("hypothesis '" + $h.name + "' has {$nk} known contributors but {$nd}" +
            " drop-out probabilities");
    }
    if ($h.unknowns < 0 or $h.unknowns > MAX_UNKNOWNS) {
        fail("hypothesis '" + $h.name + "': unknown contributors must be 0 to" +
            " {MAX_UNKNOWNS}, got {$h.unknowns} - the genotype sum grows as" +
            " (alleles²)^unknowns");
    }
    if ($nk + $h.unknowns == 0) {
        fail("hypothesis '" + $h.name + "' names no contributor");
    }
    def seen as map of string to bool init {};
    for (def i in 0..$nk) {
        checkProbability(
            "hypothesis '" + $h.name + "': a drop-out probability",
            $h.knownDropout[$i]);
        def id as string init $h.known[$i].id;
        if (maps.has($seen, $id)) {
            fail("hypothesis '" + $h.name + "' names '" + $id + "' twice");
        }
        $seen[$id] = true;
    }
    checkProbability(
        "hypothesis '" + $h.name + "': the unknowns' drop-out probability",
        $h.unknownDropout);
}

/**
 * A validated hypothesis.
 *
 * Pass empty lists for `known` and `knownDropout` when there is no known
 * contributor. A drop-out probability of 0 means the contributor's alleles
 * are always seen.
 * @param name {string} a label, used in error messages
 * @param known {list of forensics.Profile} the known contributors
 * @param knownDropout {list of float} one drop-out probability in [0, 1] per
 *   known contributor, in `known` order
 * @param unknowns {int} the number of unknown contributors, 0 to `MAX_UNKNOWNS`
 * @param unknownDropout {float} the drop-out probability, in [0, 1], shared by
 *   every unknown contributor
 * @return {Hypothesis} the hypothesis
 * @throws {Error} kind "mixture" on a length mismatch, a probability outside
 *   [0, 1], an unknown count outside 0 to `MAX_UNKNOWNS`, no contributor at
 *   all, or one profile id named twice
 */
export func hypothesis(
    name as string,
    known as list of forensics.Profile,
    knownDropout as list of float,
    unknowns as int,
    unknownDropout as float) {
    def h as Hypothesis init Hypothesis{
        name: $name,
        known: $known,
        knownDropout: $knownDropout,
        unknowns: $unknowns,
        unknownDropout: $unknownDropout
    };
    checkHypothesis($h);
    return $h;
}

/**
 * The typed profiles the θ correction conditions on: every known contributor
 * of either hypothesis, once each, Hp's first.
 *
 * A suspect named only under Hp is still typed under Hd; conditioning on him
 * there is what makes the single-source case reproduce NRC II 4.2.
 * @param hp {Hypothesis} the prosecution hypothesis
 * @param hd {Hypothesis} the defence hypothesis
 * @return {list of forensics.Profile} the typed profiles, deduplicated by id
 * @throws {Error} kind "mixture" if one id carries different genotypes in the
 *   two hypotheses
 */
export func conditioningProfiles(hp as Hypothesis, hd as Hypothesis) {
    def out as list of forensics.Profile init [];
    def at as map of string to int init {};
    for (def i in 0..len($hp.known)) {
        if (not maps.has($at, $hp.known[$i].id)) {
            $at[$hp.known[$i].id] = len($out);
            $out[] = $hp.known[$i];
        }
    }
    for (def i in 0..len($hd.known)) {
        def id as string init $hd.known[$i].id;
        if (not maps.has($at, $id)) {
            $at[$id] = len($out);
            $out[] = $hd.known[$i];
        } elseif ($out[$at[$id]].loci != $hd.known[$i].loci) {
            fail("profile '" + $id + "' has different genotypes under '" + $hp.name +
                "' and '" + $hd.name + "'");
        }
    }
    return $out;
}

# --- mixtures ---------------------------------------------------------------

/**
 * The fewest contributors that explain a mixture with no drop-in: the most
 * alleles seen at any locus, halved and rounded up.
 * @param m {forensics.Mixture} the mixture
 * @return {int} `ceil(max alleles at a locus / 2)`
 * @throws {Error} kind "mixture" on a mixture with no loci
 */
export func minimumContributors(m as forensics.Mixture) {
    if (len($m.loci) == 0) {
        fail("mixture '" + $m.id + "' has no loci");
    }
    def most as int init 0;
    for (def locus in $m.loci) {
        def seen as int init len(distinctLabels($m.loci[$locus]));
        if ($seen > $most) {
            $most = $seen;
        }
    }
    return ($most + 1) // 2;
}

/**
 * A copy of a mixture keeping only the named loci, in `keep` order.
 *
 * The likelihood functions use every locus of the mixture; narrow it here to
 * leave one out on purpose.
 * @param m {forensics.Mixture} the mixture to narrow
 * @param keep {list of string} the locus names to retain
 * @return {forensics.Mixture} the restricted mixture
 * @throws {Error} kind "mixture" if the mixture has no such locus
 */
export func restrictMixture(m as forensics.Mixture, keep as list of string) {
    def out as forensics.Mixture init forensics.Mixture{id: $m.id, loci: {}};
    for (def locus in $keep) {
        if (not maps.has($m.loci, $locus)) {
            fail("mixture '" + $m.id + "' has no locus '" + $locus + "'");
        }
        $out.loci[$locus] = $m.loci[$locus];
    }
    return $out;
}

# distinctLabels drops repeats, order preserved.
func distinctLabels(xs as list of string) {
    def seen as map of string to bool init {};
    def out as list of string init [];
    for (def x in $xs) {
        if (not maps.has($seen, $x)) {
            $seen[$x] = true;
            $out[] = $x;
        }
    }
    return $out;
}

# --- one locus --------------------------------------------------------------

func checkTheta(theta as float) {
    if ($theta < 0.0 or $theta >= 1.0) {
        fail("theta must be in [0, 1)");
    }
}

func checkDropIn(dropIn as float) {
    if ($dropIn < 0.0 or $dropIn >= 1.0) {
        fail("the drop-in probability must be in [0, 1), got {$dropIn}");
    }
}

# genotypeAtLocus fetches a typed person's genotype, failing in this module's
# kind rather than the core's, with the locus and the role named.
func genotypeAtLocus(p as forensics.Profile, locus as string, role as string) {
    if (not maps.has($p.loci, $locus)) {
        fail($role + " '" + $p.id + "' is not typed at locus '" + $locus + "'; narrow the" +
            " mixture with restrictMixture to leave the locus out");
    }
    return $p.loci[$locus];
}

# addLabel appends an allele to the label list if it is new.
func addLabel(labels as list of string, a as string) {
    if (not lists.contains($labels, $a)) {
        $labels[] = $a;
    }
    return $labels;
}

# buildFrame prepares one locus for one hypothesis. Labels are the observed
# alleles, then the typed profiles' other alleles, then (if any frequency mass
# remains) one lumped class for every other allele.
func buildFrame(
    fdb as forensics.FrequencyDb,
    observed as list of string,
    locus as string,
    h as Hypothesis,
    typed as list of forensics.Profile) {
    if (not forensics.hasLocus($fdb, $locus)) {
        fail("frequency database '" + $fdb.name + "' has no locus '" + $locus + "'");
    }
    def labels as list of string init distinctLabels($observed);
    def nObserved as int init len($labels);
    def typedG as list of forensics.Genotype init [];
    for (def i in 0..len($typed)) {
        def g as forensics.Genotype init genotypeAtLocus($typed[$i], $locus, "typed profile");
        $typedG[] = $g;
        $labels = addLabel($labels, $g.a);
        $labels = addLabel($labels, $g.b);
    }
    def nk as int init len($labels);
    if ($nk > MAX_LOCUS_ALLELES) {
        fail("locus '" + $locus + "' has {$nk} distinct alleles among the stain and the" +
            " typed profiles; at most {MAX_LOCUS_ALLELES} can be enumerated");
    }
    def f as Frame init Frame{
        labels: $labels,
        freq: [],
        observed: [],
        typedCount: [],
        typedTotal: 2.0 * convert.toFloat(len($typed)),
        knownDrop: [],
        knownCarry: [],
        gA: [],
        gB: [],
        hw: [],
        zeros: [],
        unknownPow: [1.0]
    };
    def mass as float init 0.0;
    for (def i in 0..$nk) {
        def p as float init forensics.alleleFrequency($fdb, $locus, $labels[$i]);
        $mass = $mass + $p;
        $f.freq[] = $p;
        $f.observed[] = $i < $nObserved;
        $f.typedCount[] = 0.0;
        $f.knownDrop[] = 1.0;
        $f.knownCarry[] = false;
    }
    # Floored frequencies of the enumerated alleles past 1 leave no proper
    # distribution for the unknowns; rescaling would push a floored allele below
    # its floor. Same policy as pedigree.j.
    if ($mass > 1.0 + 0.000000001) {
        fail("locus '" + $locus + "': the floored frequencies of the alleles in the stain" +
            " and the typed profiles sum to {$mass}, past 1; narrow the mixture with" +
            " restrictMixture or use a larger table");
    }
    if ($mass < 1.0 - 0.000000001) {
        $f.labels[] = "";
        $f.freq[] = 1.0 - $mass;
        $f.observed[] = false;
        $f.typedCount[] = 0.0;
        $f.knownDrop[] = 1.0;
        $f.knownCarry[] = false;
    }
    for (def g in $typedG) {
        $f.typedCount[lists.indexOf($labels, $g.a)] = $f.typedCount[lists.indexOf($labels, $g.a)] +
            1.0;
        $f.typedCount[lists.indexOf($labels, $g.b)] = $f.typedCount[lists.indexOf($labels, $g.b)] +
            1.0;
    }
    for (def i in 0..len($h.known)) {
        def g as forensics.Genotype init genotypeAtLocus($h.known[$i], $locus, "known contributor");
        def d as float init $h.knownDropout[$i];
        def ia as int init lists.indexOf($labels, $g.a);
        def ib as int init lists.indexOf($labels, $g.b);
        # one factor of d per allele copy: a homozygote's allele drops with d²
        $f.knownDrop[$ia] = $f.knownDrop[$ia] * $d;
        $f.knownDrop[$ib] = $f.knownDrop[$ib] * $d;
        $f.knownCarry[$ia] = true;
        $f.knownCarry[$ib] = true;
    }
    def nc as int init len($f.labels);
    for (def i in 0..$nc) {
        $f.zeros[] = 0;
        for (def j in $i..$nc) {
            $f.gA[] = $i;
            $f.gB[] = $j;
            $f.hw[] = hardyWeinberg($f.freq[$i], $f.freq[$j], $i == $j);
        }
    }
    # du^k for the k copies the unknowns can carry of one allele
    for (def k in 1..2 * $h.unknowns + 1) {
        $f.unknownPow[] = $f.unknownPow[$k - 1] * $h.unknownDropout;
    }
    return $f;
}

# hardyWeinberg is p² for a homozygote, 2pq for a heterozygote.
func hardyWeinberg(p as float, q as float, homozygous as bool) {
    if ($homozygous) {
        return $p * $q;
    }
    return 2.0 * $p * $q;
}

# checkTyped confirms every known contributor of `h` is among the typed
# profiles; the sampling formula is only right if their alleles are counted.
func checkTyped(h as Hypothesis, typed as list of forensics.Profile) {
    def ids as list of string init [];
    for (def i in 0..len($typed)) {
        $ids[] = $typed[$i].id;
    }
    for (def i in 0..len($h.known)) {
        if (not lists.contains($ids, $h.known[$i].id)) {
            fail("known contributor '" + $h.known[$i].id + "' of '" + $h.name +
                "' is not among the typed profiles");
        }
    }
}

# evidenceGivenGenotypes is P(E | knowns, the unknowns' genotypes `combo`)
# for the drop-out / drop-in model. `combo` holds genotype indices into the
# frame's gA / gB tables. An allele is lost only if every copy drops out.
func evidenceGivenGenotypes(f as Frame, combo as list of int, dropIn as float) {
    def nc as int init len($f.labels);
    def copies as list of int init $f.zeros;
    for (def g in $combo) {
        $copies[$f.gA[$g]] = $copies[$f.gA[$g]] + 1;
        $copies[$f.gB[$g]] = $copies[$f.gB[$g]] + 1;
    }
    def out as float init 1.0;
    def droppedIn as bool init false;
    for (def i in 0..$nc) {
        if ($f.knownCarry[$i] or $copies[$i] > 0) {
            def allGone as float init $f.knownDrop[$i] * $f.unknownPow[$copies[$i]];
            if ($f.observed[$i]) {
                $out = $out * (1.0 - $allGone);
            } else {
                $out = $out * $allGone;
            }
        } elseif ($f.observed[$i]) {
            $out = $out * $dropIn * $f.freq[$i];
            $droppedIn = true;
        }
    }
    if (not $droppedIn) {
        $out = $out * (1.0 - $dropIn);
    }
    return $out;
}

# genotypesProbability is P(the unknowns carry `combo`, in that order | the
# typed profiles) by the Balding-Nichols sampling formula: the next allele is
# `a` with probability (n_a θ + (1-θ) p_a) / (1 + (n-1) θ) after `n` alleles,
# `n_a` of them `a`. A heterozygote counts both orders. At θ = 0 this is the
# product of the precomputed Hardy-Weinberg genotype frequencies.
func genotypesProbability(f as Frame, combo as list of int, theta as float) {
    def out as float init 1.0;
    if ($theta == 0.0) {
        for (def g in $combo) {
            $out = $out * $f.hw[$g];
        }
        return $out;
    }
    # Counts of alleles drawn so far are recovered by scanning `drawn` (at most
    # 2 · MAX_UNKNOWNS entries) rather than copying a per-label count list per
    # combination.
    def drawn as list of int init [];
    def n as float init $f.typedTotal;
    for (def g in $combo) {
        $drawn[] = $f.gA[$g];
        $drawn[] = $f.gB[$g];
    }
    for (def i in 0..len($drawn)) {
        def a as int init $drawn[$i];
        def na as float init $f.typedCount[$a];
        for (def j in 0..$i) {
            if ($drawn[$j] == $a) {
                $na = $na + 1.0;
            }
        }
        $out = $out * ($na * $theta + (1.0 - $theta) * $f.freq[$a]) /
            (1.0 + ($n - 1.0) * $theta);
        $n = $n + 1.0;
    }
    for (def g in $combo) {
        if ($f.gA[$g] != $f.gB[$g]) {
            $out = $out * 2.0;
        }
    }
    return $out;
}

# orderings is the number of distinct orderings of a non-decreasing tuple:
# u! / Π (run length)!.
func orderings(combo as list of int) {
    def u as int init len($combo);
    if ($u < 2) {
        return 1;
    }
    def out as int init math.factorial($u);
    def run as int init 1;
    for (def i in 1..$u) {
        if ($combo[$i] == $combo[$i - 1]) {
            $run = $run + 1;
        } else {
            $out = $out // math.factorial($run);
            $run = 1;
        }
    }
    return $out // math.factorial($run);
}

# nextCombo advances a non-decreasing tuple of genotype indices in [0, g);
# returns an empty list once every tuple has been visited.
func nextCombo(combo as list of int, g as int) {
    def i as int init len($combo) - 1;
    while ($i >= 0 and $combo[$i] == $g - 1) {
        $i = $i - 1;
    }
    if ($i < 0) {
        return [];
    }
    def v as int init $combo[$i] + 1;
    for (def j in $i..len($combo)) {
        $combo[$j] = $v;
    }
    return $combo;
}

/**
 * P(E | H) at one locus: the probability of the alleles the stain shows there,
 * summed over every genotype combination of the unknowns.
 *
 * The lower-level form of `locusLR`, with the θ conditioning set explicit.
 * `typed` must include every known contributor of `h` (normally it is
 * `conditioningProfiles(hp, hd)`). An empty observed list is valid evidence:
 * everything dropped out.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param locus {string} the locus to score
 * @param h {Hypothesis} the hypothesis
 * @param typed {list of forensics.Profile} the profiles θ conditions on
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the probability of the evidence at the locus
 * @throws {Error} kind "mixture" on an invalid hypothesis, `dropIn` or
 *   `theta`; a locus the mixture or the table lacks; a typed profile not typed
 *   at the locus; a known contributor missing from `typed`; more than
 *   `MAX_LOCUS_ALLELES` alleles; or floored frequencies summing past 1
 */
export func evidenceProbability(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    locus as string,
    h as Hypothesis,
    typed as list of forensics.Profile,
    dropIn as float,
    theta as float) {
    checkHypothesis($h);
    checkDropIn($dropIn);
    checkTheta($theta);
    checkTyped($h, $typed);
    if (not maps.has($m.loci, $locus)) {
        fail("mixture '" + $m.id + "' has no locus '" + $locus + "'");
    }
    def f as Frame init buildFrame($fdb, $m.loci[$locus], $locus, $h, $typed);
    def nGenotypes as int init len($f.gA);
    def combo as list of int init [];
    for (def k in 0..$h.unknowns) {
        $combo[] = 0;
    }
    def total as float init 0.0;
    repeat {
        def e as float init evidenceGivenGenotypes($f, $combo, $dropIn);
        if ($e > 0.0) {
            $total = $total +
                convert.toFloat(orderings($combo)) * $e *
                genotypesProbability($f, $combo, $theta);
        }
        if (len($combo) == 0) {
            break;
        }
        $combo = nextCombo($combo, $nGenotypes);
    } until (len($combo) == 0);
    return $total;
}

/**
 * The likelihood ratio at one locus, `P(E | Hp) / P(E | Hd)`, θ conditioned on
 * `conditioningProfiles(hp, hd)`.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param locus {string} the locus to score
 * @param hp {Hypothesis} the prosecution hypothesis (numerator)
 * @param hd {Hypothesis} the defence hypothesis (denominator)
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the locus LR; 0 when Hp cannot explain the evidence
 * @throws {Error} kind "mixture" as `evidenceProbability`, or if Hd gives
 *   the evidence probability zero
 */
export func locusLR(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    locus as string,
    hp as Hypothesis,
    hd as Hypothesis,
    dropIn as float,
    theta as float) {
    def typed as list of forensics.Profile init conditioningProfiles($hp, $hd);
    def den as float init evidenceProbability($fdb, $m, $locus, $hd, $typed, $dropIn, $theta);
    if ($den <= 0.0) {
        fail("locus '" + $locus + "': '" + $hd.name + "' gives the evidence probability" +
            " zero (add an unknown, drop-out or drop-in), so the ratio is undefined");
    }
    def num as float init evidenceProbability($fdb, $m, $locus, $hp, $typed, $dropIn, $theta);
    return $num / $den;
}

/**
 * The likelihood ratio at every locus of the mixture.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param hp {Hypothesis} the prosecution hypothesis
 * @param hd {Hypothesis} the defence hypothesis
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {map of string to float} locus -> LR, in mixture order
 * @throws {Error} kind "mixture" on a mixture with no loci, or as `locusLR`
 */
export func locusLRs(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    hp as Hypothesis,
    hd as Hypothesis,
    dropIn as float,
    theta as float) {
    if (len($m.loci) == 0) {
        fail("mixture '" + $m.id + "' has no loci");
    }
    def out as map of string to float init {};
    for (def locus in $m.loci) {
        $out[$locus] = locusLR($fdb, $m, $locus, $hp, $hd, $dropIn, $theta);
    }
    return $out;
}

/**
 * The mixture likelihood ratio: `locusLR` multiplied over every locus of the
 * mixture. Every locus must be in the table and typed in every known profile;
 * narrow with `restrictMixture`. On a long panel use `log10MixtureLR`.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param hp {Hypothesis} the prosecution hypothesis
 * @param hd {Hypothesis} the defence hypothesis
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "mixture" as `locusLRs`
 */
export func mixtureLR(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    hp as Hypothesis,
    hd as Hypothesis,
    dropIn as float,
    theta as float) {
    def per as map of string to float init locusLRs($fdb, $m, $hp, $hd, $dropIn, $theta);
    def out as float init 1.0;
    for (def locus in $per) {
        $out = $out * $per[$locus];
    }
    return $out;
}

/**
 * The base-10 logarithm of the mixture likelihood ratio, summed per locus so
 * it does not underflow or overflow.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param hp {Hypothesis} the prosecution hypothesis
 * @param hd {Hypothesis} the defence hypothesis
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "mixture" as `locusLRs`, or if a locus LR is zero
 *   (Hp cannot explain the evidence there; the log is undefined)
 */
export func log10MixtureLR(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    hp as Hypothesis,
    hd as Hypothesis,
    dropIn as float,
    theta as float) {
    def per as map of string to float init locusLRs($fdb, $m, $hp, $hd, $dropIn, $theta);
    def out as float init 0.0;
    for (def locus in $per) {
        if ($per[$locus] <= 0.0) {
            fail("locus '" + $locus + "': '" + $hp.name + "' cannot explain the evidence" +
                " (LR 0), so log10 LR is undefined");
        }
        $out = $out + math.log10($per[$locus]);
    }
    return $out;
}
