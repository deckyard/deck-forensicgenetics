# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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
 * the `(k(k+1)/2)^u` cost per locus, `MAX_COMBINATIONS` the total of a call.
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
 * Three unknowns at `MAX_LOCUS_ALLELES` alleles enumerate 50116 genotype
 * combinations per hypothesis per locus (`C(68, 3)`, 66 genotypes over 11
 * classes), about 2.3 s; `locusLR` evaluates both hypotheses. `MAX_COMBINATIONS`
 * caps the total over a whole mixture.
 */
export def const MAX_UNKNOWNS as int init 3;

/**
 * The most distinct alleles one locus may enumerate: those in the stain plus
 * those of typed profiles, excluding the lumped class.
 */
export def const MAX_LOCUS_ALLELES as int init 10;

/**
 * The most unknowns' genotype combinations `locusLRs`, `mixtureLR` and
 * `log10MixtureLR` enumerate for one mixture, summed over every locus and both
 * hypotheses; checked before any locus is evaluated.
 *
 * A combination costs about 35-45 µs, so the cap is about 11 s. A 24-locus
 * three-person stain with six alleles per locus (Hp suspect plus two
 * unknowns, Hd three unknowns) needs 107184 (3.7 s). Past the cap, score the
 * loci in groups with `restrictMixture` and multiply.
 */
export def const MAX_COMBINATIONS as int init 250000;

# FLOAT_MAX is the largest finite float. The interpreter raises a runtime error
# on a product past it, so `mixtureLR` checks before each multiplication.
def const FLOAT_MAX as float init 1.7976931348623157e308;

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
# depend on those genotypes. `seen + flip · g` is the factor for class i when
# all its copies drop out with probability g: `1 - g` observed, `g` not.
# `draw[i · width + k]` is the sampling-formula
# numerator of the unknowns drawing k copies of class i, `drawDenom` the
# denominator of all 2u draws (`width` = 2u + 1).
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
    unknowns as int,
    unknownPow as list of float,
    width as int,
    draw as list of float,
    drawDenom as float,
    seen as list of float,
    flip as list of float
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

# frequencyOf is the core's floored frequency. A table with no sample size
# (n = 0) has no floor, so an allele it does not list throws here, in this
# module's kind.
func frequencyOf(fdb as forensics.FrequencyDb, locus as string, allele as string) {
    if (forensics.sampleSize($fdb, $locus) == 0 and
        (not maps.has($fdb.loci[$locus].alleles, $allele) or
        $fdb.loci[$locus].alleles[$allele] <= 0.0)) {
        fail("locus '" + $locus + "' allele '" + $allele + "': no frequency in '" +
            $fdb.name + "', and no 5/2N floor (sample size 0)");
    }
    return forensics.alleleFrequency($fdb, $locus, $allele);
}

# buildFrame prepares one locus for one hypothesis. Labels are the observed
# alleles, then the typed profiles' other alleles, then (if any frequency mass
# remains) one lumped class for every other allele.
func buildFrame(
    fdb as forensics.FrequencyDb,
    observed as list of string,
    locus as string,
    h as Hypothesis,
    typed as list of forensics.Profile,
    theta as float) {
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
        unknowns: $h.unknowns,
        unknownPow: [1.0],
        width: 2 * $h.unknowns + 1,
        draw: [],
        drawDenom: 1.0,
        seen: [],
        flip: []
    };
    def mass as float init 0.0;
    for (def i in 0..$nk) {
        def p as float init frequencyOf($fdb, $locus, $labels[$i]);
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
        if ($f.observed[$i]) {
            $f.seen[] = 1.0;
            $f.flip[] = -1.0;
        } else {
            $f.seen[] = 0.0;
            $f.flip[] = 1.0;
        }
        for (def j in $i..$nc) {
            $f.gA[] = $i;
            $f.gB[] = $j;
        }
    }
    # du^k for the k copies the unknowns can carry of one allele
    for (def k in 1..$f.width) {
        $f.unknownPow[] = $f.unknownPow[$k - 1] * $h.unknownDropout;
    }
    $f.draw = drawTable($f, $theta);
    for (def k in 0..2 * $h.unknowns) {
        $f.drawDenom = $f.drawDenom *
            (1.0 + ($f.typedTotal + convert.toFloat($k) - 1.0) * $theta);
    }
    return $f;
}

# drawTable holds the sampling-formula numerators. Balding-Nichols: after n
# alleles, n_a of them `a`, the next is `a` with probability
# (n_a θ + (1-θ) p_a) / (1 + (n-1) θ). The product over the unknowns' 2u draws
# depends only on how many copies of each class they draw, so the numerators
# are tabled per (class, copies) and the denominators multiply once
# (`drawDenom`). At θ = 0 the table holds p_a^k.
func drawTable(f as Frame, theta as float) {
    def out as list of float init [];
    for (def i in 0..len($f.labels)) {
        def acc as float init 1.0;
        $out[] = $acc;
        for (def k in 1..$f.width) {
            $acc = $acc *
                (($f.typedCount[$i] + convert.toFloat($k - 1)) * $theta +
                (1.0 - $theta) * $f.freq[$i]);
            $out[] = $acc;
        }
    }
    return $out;
}

# checkTyped confirms every known contributor of `h` is among the typed
# profiles with the same genotype at `locus`; the sampling formula is only
# right if the alleles it counts are the ones the hypothesis uses.
func checkTyped(h as Hypothesis, typed as list of forensics.Profile, locus as string) {
    for (def i in 0..len($h.known)) {
        def id as string init $h.known[$i].id;
        def at as int init -1;
        for (def k in 0..len($typed)) {
            if ($typed[$k].id == $id) {
                $at = $k;
                break;
            }
        }
        if ($at < 0) {
            fail("known contributor '" + $id + "' of '" + $h.name +
                "' is not among the typed profiles");
        }
        if (maps.has($typed[$at].loci, $locus) and maps.has($h.known[$i].loci, $locus) and
            $typed[$at].loci[$locus] != $h.known[$i].loci[$locus]) {
            fail("profile '" + $id + "' has different genotypes at locus '" + $locus +
                "' in typed and in '" + $h.name + "'");
        }
    }
}

# multisets is the number of non-decreasing u-tuples over g genotypes,
# C(g + u - 1, u): the combinations `evaluate` visits.
func multisets(g as int, u as int) {
    def out as int init 1;
    for (def k in 1..$u + 1) {
        # out is C(g + k - 2, k - 1) here, so the division is exact
        $out = $out * ($g + $k - 1) // $k;
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

# evaluate sums P(E | knowns, G; d, c) · P(G | typed; θ) over the unknowns'
# genotypes G, enumerated as non-decreasing tuples of genotype indices weighted
# by their orderings. `copies` (per class) and `hetFactor` (2 per heterozygote,
# for both allele orders) follow the tuple as it advances, so a combination
# costs one pass over the classes and allocates nothing.
func evaluate(f as Frame, dropIn as float) {
    def u as int init $f.unknowns;
    def nc as int init len($f.labels);
    def last as int init len($f.gA) - 1;
    def w as int init $f.width;
    def combo as list of int init [];
    def copies as list of int init [];
    for (def i in 0..$nc) {
        $copies[] = 0;
    }
    for (def k in 0..$u) {
        $combo[] = 0;
    }
    # genotype 0 is class 0 twice
    $copies[0] = 2 * $u;
    def hetFactor as float init 1.0;
    def total as float init 0.0;
    while (true) {
        # An allele is lost only if every copy drops out; an observed allele
        # nobody carries is a drop-in.
        def e as float init 1.0;
        def draws as float init 1.0;
        def droppedIn as bool init false;
        for (def i in 0..$nc) {
            def c as int init $copies[$i];
            if ($f.knownCarry[$i] or $c > 0) {
                $e = $e *
                    ($f.seen[$i] + $f.flip[$i] * $f.knownDrop[$i] * $f.unknownPow[$c]);
                $draws = $draws * $f.draw[$i * $w + $c];
            } elseif ($f.observed[$i]) {
                $e = $e * $dropIn * $f.freq[$i];
                $droppedIn = true;
            }
        }
        if (not $droppedIn) {
            $e = $e * (1.0 - $dropIn);
        }
        if ($e > 0.0) {
            $total = $total +
                convert.toFloat(orderings($combo)) * $hetFactor * $e * $draws /
                $f.drawDenom;
        }
        # next non-decreasing tuple: bump the rightmost entry below `last` and
        # set every entry after it to the same value
        def i as int init $u - 1;
        while ($i >= 0 and $combo[$i] == $last) {
            $i = $i - 1;
        }
        if ($i < 0) {
            break;
        }
        def v as int init $combo[$i] + 1;
        for (def j in $i..$u) {
            def old as int init $combo[$j];
            $copies[$f.gA[$old]] = $copies[$f.gA[$old]] - 1;
            $copies[$f.gB[$old]] = $copies[$f.gB[$old]] - 1;
            if ($f.gA[$old] != $f.gB[$old]) {
                $hetFactor = $hetFactor / 2.0;
            }
            $combo[$j] = $v;
            $copies[$f.gA[$v]] = $copies[$f.gA[$v]] + 1;
            $copies[$f.gB[$v]] = $copies[$f.gB[$v]] + 1;
            if ($f.gA[$v] != $f.gB[$v]) {
                $hetFactor = $hetFactor * 2.0;
            }
        }
    }
    return $total;
}

# frameFor checks one hypothesis at one locus and prepares its frame.
func frameFor(
    fdb as forensics.FrequencyDb,
    m as forensics.Mixture,
    locus as string,
    h as Hypothesis,
    typed as list of forensics.Profile,
    theta as float) {
    checkTyped($h, $typed, $locus);
    if (not maps.has($m.loci, $locus)) {
        fail("mixture '" + $m.id + "' has no locus '" + $locus + "'");
    }
    return buildFrame($fdb, $m.loci[$locus], $locus, $h, $typed, $theta);
}

# checkDefence rejects a zero denominator before the numerator is computed.
func checkDefence(locus as string, hd as Hypothesis, den as float) {
    if ($den <= 0.0) {
        fail("locus '" + $locus + "': '" + $hd.name + "' gives the evidence probability" +
            " zero (add an unknown, drop-out or drop-in), so the ratio is undefined");
    }
}

/**
 * P(E | H) at one locus: the probability of the alleles the stain shows there,
 * summed over every genotype combination of the unknowns.
 *
 * The lower-level form of `locusLR`, with the θ conditioning set explicit.
 * `typed` must include every known contributor of `h`, with the same genotype
 * at the locus (normally it is `conditioningProfiles(hp, hd)`). An empty
 * observed list, which only a `Mixture` struct literal can hold, means
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
 *   at the locus; a known contributor missing from `typed` or typed there with
 *   a different genotype; more than `MAX_LOCUS_ALLELES` alleles; an allele the
 *   table does not list at a locus with sample size 0; or floored frequencies
 *   summing past 1
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
    return evaluate(frameFor($fdb, $m, $locus, $h, $typed, $theta), $dropIn);
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
    checkDefence($locus, $hd, $den);
    def num as float init evidenceProbability($fdb, $m, $locus, $hp, $typed, $dropIn, $theta);
    return $num / $den;
}

/**
 * The likelihood ratio at every locus of the mixture.
 *
 * Every locus is checked and its genotype combinations counted before any is
 * evaluated.
 * @param fdb {forensics.FrequencyDb} the reference frequencies
 * @param m {forensics.Mixture} the stain
 * @param hp {Hypothesis} the prosecution hypothesis
 * @param hd {Hypothesis} the defence hypothesis
 * @param dropIn {float} the drop-in probability `c`, in [0, 1)
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {map of string to float} locus -> LR, in mixture order
 * @throws {Error} kind "mixture" on a mixture with no loci, more than
 *   `MAX_COMBINATIONS` genotype combinations in total, or as `locusLR`
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
    def typed as list of forensics.Profile init conditioningProfiles($hp, $hd);
    checkHypothesis($hp);
    checkHypothesis($hd);
    checkDropIn($dropIn);
    checkTheta($theta);
    # Frames alternate Hd, Hp per locus.
    def frames as list of Frame init [];
    def combos as int init 0;
    for (def locus in $m.loci) {
        $frames[] = frameFor($fdb, $m, $locus, $hd, $typed, $theta);
        $frames[] = frameFor($fdb, $m, $locus, $hp, $typed, $theta);
        def pair as int init len($frames) - 2;
        $combos = $combos + multisets(len($frames[$pair].gA), $hd.unknowns) +
            multisets(len($frames[$pair + 1].gA), $hp.unknowns);
    }
    if ($combos > MAX_COMBINATIONS) {
        fail("mixture '" + $m.id + "' needs {$combos} genotype combinations over its loci" +
            " and both hypotheses, past MAX_COMBINATIONS ({MAX_COMBINATIONS}); score fewer" +
            " loci per call (restrictMixture) or name fewer unknowns");
    }
    def out as map of string to float init {};
    def at as int init 0;
    for (def locus in $m.loci) {
        def den as float init evaluate($frames[$at], $dropIn);
        checkDefence($locus, $hd, $den);
        $out[$locus] = evaluate($frames[$at + 1], $dropIn) / $den;
        $at = $at + 2;
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
 * @throws {Error} kind "mixture" as `locusLRs`, or if the product exceeds the
 *   float range (use `log10MixtureLR`)
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
        # `out <= 1` cannot overflow; for `out > 1` the quotient is finite
        if ($out > 1.0 and $per[$locus] > FLOAT_MAX / $out) {
            fail("the mixture LR of '" + $m.id + "' exceeds the float range (above" +
                " 1.8e308); use log10MixtureLR");
        }
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
