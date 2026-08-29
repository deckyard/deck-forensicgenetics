# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <jennifer@mplx.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Lineage markers: mitochondrial DNA and Y-STR haplotypes.
 *
 * mtDNA and the Y chromosome are haploid, non-recombining and uniparentally
 * inherited, so the entire marker set is inherited as one **haplotype**. There
 * is no Hardy-Weinberg equilibrium to invoke and no product rule across loci:
 * the match probability is a direct **counting** estimate from a reference
 * database - the haplotype was seen `k` times among `N` profiles.
 *
 * That counting estimate is unstable for rare haplotypes (a haplotype seen
 * once, or never, in a database of a few thousand), so the conservative
 * figure to quote is the **Clopper-Pearson exact-binomial upper bound** rather
 * than the raw `k/N`. `lineageMatchProbability` returns that bound.
 *
 * The module is deliberately **database-independent**: it supplies the
 * estimators and the haplotype types, and the caller feeds in `k` and `N` from
 * whichever database they queried (EMPOP, YHRD, an in-house set). No reference
 * data is bundled, and neither is the rCRS - `align` takes the reference
 * region as an argument.
 * @module lineage
 * @see https://empop.online/ - the EMPOP mtDNA database
 * @see https://yhrd.org/ - the Y-STR Haplotype Reference Database
 * @example
 * import "@jennifer/forensicgenetics/lineage.j" as lineage;
 *
 * # A haplotype seen 2 times in a database of 5000.
 * def p as float init lineage.lineageMatchProbability(2, 5000);   # ~0.00133
 * def lr as float init lineage.lineageLikelihoodRatio(2, 5000);   # ~751
 *
 * def h as lineage.MtHaplotype init lineage.mtHaplotype("Q1", ["16024-16365"],
 *     "16189C 16223T 16278T");
 */
use math;
use maps;
use lists;
use stats;
use strings;
use convert;

/**
 * The default confidence level for a conservative lineage estimate: 0.95,
 * i.e. a one-sided 95% upper bound.
 */
export def const DEFAULT_CONFIDENCE as float init 0.95;

# The largest reference position a difference designation may name. The human
# mitochondrial genome is 16569 bases, so this is far above anything real; the
# bound exists because positions are sorted through a fixed-width key, which a
# wider number would order wrongly.
def const MAX_POSITION as int init 999999;

/**
 * One difference from the reference sequence, in rCRS nomenclature.
 *
 * A substitution at 16189 to C is `{position: 16189, insertion: 0, base: "C",
 * kind: "sub"}` and prints `16189C`. The first base inserted after 315 is
 * `{position: 315, insertion: 1, base: "C", kind: "ins"}` and prints `315.1C`.
 * A deletion of 249 is `{position: 249, insertion: 0, base: "", kind: "del"}`
 * and prints `249d`.
 * @field position {int} the reference position
 * @field insertion {int} 0 for a substitution or deletion; 1, 2, ... for the
 *   n-th base inserted after `position`
 * @field base {string} the observed base; "" for a deletion
 * @field kind {string} "sub", "ins" or "del"
 */
export def struct Difference {
    position as int,
    insertion as int,
    base as string,
    kind as string
};

/**
 * An mtDNA haplotype: the sequenced ranges plus the differences from the
 * reference within them.
 *
 * The ranges matter for comparison. Two haplotypes can only be compared
 * like-for-like over regions both were sequenced across, so a haplotype that
 * records only HVR-I is not "equal" to one covering the whole control region
 * merely because their listed differences coincide.
 * @field id {string} a caller-chosen label
 * @field ranges {list of string} the sequenced regions, e.g. ["16024-16365", "73-340"]
 * @field differences {list of Difference} the differences from the reference,
 *   in ascending position order
 */
export def struct MtHaplotype {
    id as string,
    ranges as list of string,
    differences as list of Difference
};

/**
 * A Y-STR haplotype: a repeat count per locus.
 *
 * A multi-copy locus (DYS385a/b, and DYS389I/II when reported together) holds
 * its alleles comma-separated - `"11,14"`. Copy counts must agree between two
 * haplotypes being compared at that locus.
 * @field id {string} a caller-chosen label
 * @field loci {map of string to string} locus name -> allele designation(s)
 */
export def struct YHaplotype {
    id as string,
    loci as map of string to string
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "lineage", message: $msg, file: "", line: 0, col: 0};
}

# --- counting estimators ----------------------------------------------------

# checkCount validates a (k, N) observation pair.
func checkCount(k as int, n as int) {
    if ($n <= 0) {
        fail("database size N must be positive");
    }
    if ($k < 0 or $k > $n) {
        fail("observed count k must be in [0, N]");
    }
}

/**
 * The raw counting estimate of a haplotype's frequency: `k / N`.
 *
 * The point estimate, not the figure to quote in a report - a haplotype seen
 * zero times gets 0.0 here, which no reference database can actually support.
 * Use `lineageMatchProbability` for the conservative estimate.
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @return {float} the point estimate k/N
 * @throws {Error} kind "lineage" on a non-positive N or a k outside [0, N]
 */
export func haplotypeFrequency(k as int, n as int) {
    checkCount($k, $n);
    return convert.toFloat($k) / convert.toFloat($n);
}

/**
 * The Clopper-Pearson exact-binomial **one-sided upper bound** on a haplotype
 * frequency.
 *
 * The conservative counting estimate: the largest population frequency that
 * would still make `k` sightings in `N` profiles unsurprising at the given
 * confidence. It is well defined at `k == 0`, where it reduces to the familiar
 * `1 - alpha^(1/N)` (the "rule of three" at 95%), and it is exactly 1.0 at
 * `k == N`.
 *
 * Computed as the upper end of a two-sided Clopper-Pearson interval at level
 * `2c - 1`, which is the standard identity between a two-sided interval and a
 * one-sided bound.
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @param confidence {float} the one-sided confidence level, in (0.5, 1)
 * @return {float} the upper bound on the frequency
 * @throws {Error} kind "lineage" on a bad count pair or a confidence outside (0.5, 1)
 */
export func haplotypeUpperBound(k as int, n as int, confidence as float) {
    checkCount($k, $n);
    if ($confidence <= 0.5 or $confidence >= 1.0) {
        fail("confidence must be in (0.5, 1)");
    }
    if ($k == $n) {
        return 1.0;
    }
    def twoSided as float init 2.0 * $confidence - 1.0;
    def ci as stats.Interval init stats.proportionCi($k, $n, $twoSided, "clopper-pearson");
    return $ci.upper;
}

/**
 * The conservative lineage match probability: the 95% upper bound on `k / N`.
 *
 * This is the number a lineage-marker report quotes. It is `haplotypeUpperBound`
 * at `DEFAULT_CONFIDENCE`.
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @return {float} the conservative match probability
 * @throws {Error} kind "lineage" on a non-positive N or a k outside [0, N]
 */
export func lineageMatchProbability(k as int, n as int) {
    return haplotypeUpperBound($k, $n, DEFAULT_CONFIDENCE);
}

/**
 * The lineage likelihood ratio: the reciprocal of `lineageMatchProbability`.
 *
 * The weight of evidence for "the sample came from this maternal (mtDNA) or
 * paternal (Y) lineage" against "it came from an unrelated lineage".
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @return {float} the likelihood ratio
 * @throws {Error} kind "lineage" on a bad count pair, or if the bound is 1.0
 *   (`k == N`), which carries no weight
 */
export func lineageLikelihoodRatio(k as int, n as int) {
    def p as float init lineageMatchProbability($k, $n);
    if ($p <= 0.0) {
        fail("a match probability of zero has no likelihood ratio");
    }
    return 1.0 / $p;
}

# --- rCRS difference nomenclature -------------------------------------------

# isDigit / isBase classify one character of a difference designation. The base
# set is the IUPAC nucleotide alphabet, so an ambiguity code from a
# heteroplasmic position still parses.
func isDigit(ch as string) {
    return strings.indexOf("0123456789", $ch) >= 0;
}

func isBase(ch as string) {
    return strings.indexOf("ACGTURYSWKMBDHVN", $ch) >= 0;
}

/**
 * Parse one difference designation in rCRS nomenclature.
 *
 * Accepts `16189C`, the redundant `A263G` form (the reference base is dropped),
 * the insertion form `315.1C`, and the deletion forms `249d`, `249D`, `249DEL`
 * and `249-`.
 * @param text {string} the designation
 * @return {Difference} the parsed difference
 * @throws {Error} kind "lineage" if the designation does not parse
 */
export func parseDifference(text as string) {
    def t as string init strings.upper(strings.trim($text));
    if (len($t) == 0) {
        fail("empty difference designation");
    }
    def ch as list of string init strings.chars($t);
    def i as int init 0;
    # An explicit reference base is redundant (the position determines it), so
    # "A263G" and "263G" are the same difference. Drop the leading base.
    if (isBase($ch[0]) and len($ch) > 1 and isDigit($ch[1])) {
        $i = 1;
    }
    def digits as string init "";
    while ($i < len($ch) and isDigit($ch[$i])) {
        $digits = $digits + $ch[$i];
        $i = $i + 1;
    }
    if (len($digits) == 0) {
        fail("difference '" + $text + "' has no position");
    }
    if (len($digits) > 6) {
        fail("difference '" + $text + "': position exceeds {MAX_POSITION}");
    }
    def position as int init convert.toInt($digits);
    if ($position < 1 or $position > MAX_POSITION) {
        fail("difference '" + $text + "': position must be in [1, {MAX_POSITION}]");
    }
    def insertion as int init 0;
    if ($i < len($ch) and $ch[$i] == ".") {
        $i = $i + 1;
        def more as string init "";
        while ($i < len($ch) and isDigit($ch[$i])) {
            $more = $more + $ch[$i];
            $i = $i + 1;
        }
        if (len($more) == 0) {
            fail("difference '" + $text + "' has a '.' with no insertion index");
        }
        if (len($more) > 4) {
            fail("difference '" + $text + "': insertion index is implausibly large");
        }
        $insertion = convert.toInt($more);
        if ($insertion < 1) {
            fail("difference '" + $text + "': an insertion index starts at 1");
        }
    }
    def rest as string init strings.join($ch[$i..], "");
    if (len($rest) == 0) {
        fail("difference '" + $text + "' names neither a base nor a deletion");
    }
    if ($rest == "D" or $rest == "DEL" or $rest == "-") {
        if ($insertion > 0) {
            fail("difference '" + $text + "': an inserted position cannot also be deleted");
        }
        return Difference{position: $position, insertion: 0, base: "", kind: "del"};
    }
    if (len($rest) != 1 or not isBase($rest)) {
        fail("difference '" + $text + "': '" + $rest + "' is not a single IUPAC base");
    }
    def kind as string init "sub";
    if ($insertion > 0) {
        $kind = "ins";
    }
    return Difference{position: $position, insertion: $insertion, base: $rest, kind: $kind};
}

/**
 * Render a difference back into rCRS nomenclature.
 * @param d {Difference} the difference
 * @return {string} its designation, e.g. "16189C", "315.1C" or "249d"
 * @throws {Error} kind "lineage" on an unknown `kind`
 */
export func formatDifference(d as Difference) {
    match ($d.kind) {
        when "sub" { return "{$d.position}{$d.base}"; }
        when "ins" { return "{$d.position}.{$d.insertion}{$d.base}"; }
        when "del" { return "{$d.position}d"; }
        else { fail("unknown difference kind '" + $d.kind + "'"); }
    }
}

/**
 * An mtDNA haplotype from a space-separated list of differences.
 *
 * Tabs and newlines separate as well as spaces, so a designation string copied
 * out of a report reads directly.
 * @param id {string} a label for the haplotype
 * @param ranges {list of string} the sequenced regions, e.g. ["16024-16365"]
 * @param text {string} the differences, e.g. "16189C 16223T 315.1C"
 * @return {MtHaplotype} the haplotype, differences sorted by position
 * @throws {Error} kind "lineage" if a designation does not parse
 */
export func mtHaplotype(id as string, ranges as list of string, text as string) {
    def flat as string init strings.replace(strings.replace($text, "\t", " "), "\n", " ");
    def out as list of Difference init [];
    for (def token in strings.split($flat, " ")) {
        if (len(strings.trim($token)) > 0) {
            $out[] = parseDifference($token);
        }
    }
    return MtHaplotype{id: $id, ranges: $ranges, differences: sortDifferences($out)};
}

# sortKey renders a difference as a zero-padded sortable string: position,
# then insertion index, so 315 < 315.1 < 316 and 3150 does not sort between
# them.
func sortKey(d as Difference) {
    return pad(convert.toString($d.position), 6) + pad(convert.toString($d.insertion), 4);
}

func pad(s as string, width as int) {
    def out as string init $s;
    while (len($out) < $width) {
        $out = "0" + $out;
    }
    return $out;
}

# sortDifferences orders differences by position then insertion index.
func sortDifferences(ds as list of Difference) {
    return lists.sortBy($ds, sortKey);
}

/**
 * The differences of a haplotype as designation strings.
 * @param h {MtHaplotype} the haplotype
 * @return {list of string} the designations in position order
 */
export func mtDesignations(h as MtHaplotype) {
    def out as list of string init [];
    for (def d in $h.differences) {
        $out[] = formatDifference($d);
    }
    return $out;
}

/**
 * Whether two mtDNA haplotypes carry exactly the same differences.
 *
 * This compares difference lists only. It does **not** check that the two were
 * sequenced over the same regions - test that separately with `mtRangesAgree`,
 * because comparing an HVR-I-only haplotype against a whole-control-region one
 * is not a like-for-like match.
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {bool} true if the difference sets are identical
 */
export func mtMatches(a as MtHaplotype, b as MtHaplotype) {
    return mtDistance($a, $b) == 0;
}

/**
 * The number of differences that separate two mtDNA haplotypes.
 *
 * The size of the symmetric difference of the two designation sets: a position
 * one carries and the other does not counts once, and a position both carry
 * with different bases counts twice (each haplotype has a designation the other
 * lacks).
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {int} the number of differing designations
 */
export func mtDistance(a as MtHaplotype, b as MtHaplotype) {
    def da as list of string init mtDesignations($a);
    def db as list of string init mtDesignations($b);
    def inA as map of string to bool init setOf($da);
    def inB as map of string to bool init setOf($db);
    def n as int init 0;
    for (def d in $da) {
        if (not maps.has($inB, $d)) {
            $n = $n + 1;
        }
    }
    for (def d in $db) {
        if (not maps.has($inA, $d)) {
            $n = $n + 1;
        }
    }
    return $n;
}

# setOf indexes a list for membership testing, so a comparison is a hash lookup
# rather than a scan.
func setOf(xs as list of string) {
    def out as map of string to bool init {};
    for (def x in $xs) {
        $out[$x] = true;
    }
    return $out;
}

/**
 * Whether two haplotypes declare the same sequenced ranges, in any order.
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {bool} true if each declares exactly the ranges the other does
 */
export func mtRangesAgree(a as MtHaplotype, b as MtHaplotype) {
    if (len($a.ranges) != len($b.ranges)) {
        return false;
    }
    for (def r in $a.ranges) {
        if (not lists.contains($b.ranges, $r)) {
            return false;
        }
    }
    return true;
}

# --- alignment to the reference ---------------------------------------------

# Needleman-Wunsch scores. A gap costs more than a mismatch, so an isolated
# difference is called a substitution rather than an insertion plus a deletion.
def const MATCH as int init 1;
def const MISMATCH as int init -1;
def const GAP as int init -2;

# The traceback move stored for each cell. Only the direction is kept - one
# byte per cell rather than a full score matrix - which is what lets a
# control-region alignment run in megabytes instead of gigabytes.
def const DIR_DIAG as int init 0;
def const DIR_UP as int init 1;
def const DIR_LEFT as int init 2;

/**
 * The longest region `align` will accept, in bases.
 *
 * Alignment is O(nm) in time, so the cost grows with the square of the region:
 * a 1000-base region is a second or so, and the whole 16.5 kb mitochondrial
 * genome would be well over half an hour. The cap turns that into an immediate,
 * catchable refusal rather than a process that appears to hang. It is far above
 * the control region (~1.1 kb), which is what this nomenclature is for.
 */
export def const MAX_ALIGN_LENGTH as int init 5000;

# clean strips whitespace and upper-cases a sequence, so a FASTA block pasted
# with line breaks aligns as one string.
func clean(seq as string) {
    # Collected into a list and joined once. Appending to a string would
    # reallocate and copy the whole prefix on every base, which is quadratic in
    # the region length.
    def out as list of string init [];
    for (def ch in strings.chars(strings.upper($seq))) {
        if (isBase($ch)) {
            $out[] = $ch;
        } elseif (strings.indexOf(" \t\r\n", $ch) < 0) {
            fail("sequence contains a non-nucleotide character '" + $ch + "'");
        }
    }
    return strings.join($out, "");
}

# A cell outside the computed band. Low enough that no real path can reach it,
# and far enough from the int64 edge that adding GAP to it cannot overflow.
def const NEG as int init -1000000000;

# The first band half-width tried. Two aligned mtDNA sequences differ by a
# handful of positions, so the optimal path hugs the diagonal and a narrow band
# holds it; the width doubles until the path no longer touches the edge.
def const INITIAL_BAND as int init 16;

# Band is one banded pass: the two gapped strings, and whether the band was wide
# enough for the result to be the true global optimum.
def struct Band {
    exact as bool,
    top as string,
    bottom as string
};

# alignPair returns the two gapped strings as [gappedRef, gappedSample].
#
# Full Needleman-Wunsch is O(nm), which for a control region is a million cells
# of interpreted arithmetic. Instead the matrix is computed only within a band
# around the diagonal, which is where the optimal path lies when the sequences
# are near-identical - as an mtDNA haplotype and its reference are. If the
# traceback touches the band edge the band may have constrained it, so the width
# doubles and the pass repeats; at full width the band is the whole matrix, so
# the result is always the true optimum, never an approximation.
func alignPair(reference as string, sample as string) {
    def r as list of string init strings.chars($reference);
    def s as list of string init strings.chars($sample);
    def rb as bytes init convert.bytesFromString($reference, "utf-8");
    def sb as bytes init convert.bytesFromString($sample, "utf-8");
    def limit as int init math.max(len($r), len($s));
    def w as int init INITIAL_BAND;
    if ($w > $limit) {
        $w = $limit;
    }
    def band as Band init alignBand($rb, $sb, $r, $s, $w);
    while (not $band.exact and $w < $limit) {
        $w = $w * 2;
        if ($w > $limit) {
            $w = $limit;
        }
        $band = alignBand($rb, $sb, $r, $s, $w);
    }
    return [$band.top, $band.bottom];
}

# alignBand fills the cells within `w` of the diagonal and traces back through
# them. Two score rows are live at a time and each row keeps one direction byte
# per band cell, so both time and memory are O(w * m) rather than O(nm).
#
# The innermost loop does as little as possible: bases are compared as bytes
# rather than one-character strings, the row's own base is hoisted out, and the
# score to the left is carried in a variable instead of being re-indexed.
func alignBand(rb as bytes, sb as bytes, r as list of string, s as list of string, w as int) {
    def m as int init len($r);
    def n as int init len($s);
    if (math.abs($m - $n) > $w) {
        # The end cell lies outside the band, so this width cannot even express
        # the alignment. Ask for a wider one.
        return Band{exact: false, top: "", bottom: ""};
    }
    def prev as list of int init [];
    for (def j in 0..($n + 2)) {
        $prev[] = NEG;
    }
    def cur as list of int init $prev;
    def dirs as list of bytes init [];
    def los as list of int init [];
    def firstHi as int init math.min($n, $w);
    def firstRow as bytes;
    $firstRow[] = DIR_DIAG;
    $prev[0] = 0;
    for (def j in 1..($firstHi + 1)) {
        $prev[$j] = $j * GAP;
        $firstRow[] = DIR_LEFT;
    }
    $prev[$firstHi + 1] = NEG;
    $dirs[] = $firstRow;
    $los[] = 0;
    for (def i in 1..($m + 1)) {
        def lo as int init math.max(0, $i - $w);
        def hi as int init math.min($n, $i + $w);
        def base as int init $rb[$i - 1];
        def moves as bytes;
        def left as int init NEG;
        if ($lo > 0) {
            # Guard the cell just left of the band: the next row reads it as its
            # diagonal, and it must not see a value two rows stale.
            $cur[$lo - 1] = NEG;
        } else {
            $cur[0] = $i * GAP;
            $moves[] = DIR_UP;
            $left = $i * GAP;
        }
        for (def j in math.max(1, $lo)..($hi + 1)) {
            def up as int init $prev[$j];
            def diag as int init $prev[$j - 1];
            def best as int init $diag + MISMATCH;
            if ($base == $sb[$j - 1]) {
                $best = $diag + MATCH;
            }
            def move as int init DIR_DIAG;
            # Ties keep the earlier branch, so this is the preference order:
            # diagonal, then deletion, then insertion. A run of equal-scoring
            # cells reads as substitutions rather than as an indel pair.
            if ($up + GAP > $best) {
                $best = $up + GAP;
                $move = DIR_UP;
            }
            if ($left + GAP > $best) {
                $best = $left + GAP;
                $move = DIR_LEFT;
            }
            $moves[] = $move;
            $cur[$j] = $best;
            $left = $best;
        }
        $cur[$hi + 1] = NEG;
        $dirs[] = $moves;
        $los[] = $lo;
        def spare as list of int init $prev;
        $prev = $cur;
        $cur = $spare;
    }
    return traceback($dirs, $los, $r, $s, $w);
}

# traceback walks the stored moves back to the origin, reporting whether the
# path ever touched the band edge - which is what tells `alignPair` the width
# may have constrained the result.
func traceback(
    dirs as list of bytes,
    los as list of int,
    r as list of string,
    s as list of string,
    w as int) {
    def i as int init len($r);
    def j as int init len($s);
    def exact as bool init true;
    def gappedRef as list of string init [];
    def gappedSample as list of string init [];
    while ($i > 0 or $j > 0) {
        if ($i - $j >= $w or $j - $i >= $w) {
            $exact = false;
        }
        def move as int init $dirs[$i][$j - $los[$i]];
        if ($i > 0 and $j > 0 and $move == DIR_DIAG) {
            $gappedRef[] = $r[$i - 1];
            $gappedSample[] = $s[$j - 1];
            $i = $i - 1;
            $j = $j - 1;
        } elseif ($i > 0 and $move == DIR_UP) {
            $gappedRef[] = $r[$i - 1];
            $gappedSample[] = "-";
            $i = $i - 1;
        } else {
            $gappedRef[] = "-";
            $gappedSample[] = $s[$j - 1];
            $j = $j - 1;
        }
    }
    return Band{
        exact: $exact,
        top: strings.join(lists.reverse($gappedRef), ""),
        bottom: strings.join(lists.reverse($gappedSample), "")
    };
}

# shiftGapsRight slides every gap in `a` as far 3' as it can go without
# breaking the alignment, which is what puts an insertion into a homopolymer
# run at the run's 3' end - the convention that makes an extra C in the
# 311-315 stretch read as 315.1C rather than 311.1C.
func shiftGapsRight(a as string, b as string) {
    def x as list of string init strings.chars($a);
    def y as list of string init strings.chars($b);
    def moved as bool init true;
    while ($moved) {
        $moved = false;
        for (def i in 0..(len($x) - 1)) {
            if ($x[$i] == "-" and $x[$i + 1] != "-" and $x[$i + 1] == $y[$i]) {
                $x[$i] = $x[$i + 1];
                $x[$i + 1] = "-";
                $moved = true;
            }
        }
    }
    return strings.join($x, "");
}

/**
 * Align a sequenced region against its reference and render the result as rCRS
 * differences.
 *
 * This is the substantive step that turns a raw sequence into standard
 * nomenclature: a global Needleman-Wunsch alignment (match +1, mismatch -1,
 * gap -2), with indels slid to the 3' end of any homopolymer run so an extra
 * base in the 311-315 C-stretch reads as `315.1C`.
 *
 * No reference is bundled - pass the rCRS (NC_012920.1) region you sequenced,
 * with `refStart` its first position, so the returned differences are numbered
 * in rCRS coordinates.
 *
 * **Caveat**: this renders a plain optimal alignment. EMPOP's phylogenetic
 * alignment rules, which shift indels to the position preferred by the mtDNA
 * phylogeny rather than merely 3', are not applied - so a haplotype with a
 * complex indel may differ from an EMPOP-conformant designation and should be
 * checked before it is quoted.
 * @param reference {string} the reference sequence for the region
 * @param sample {string} the sequenced sample over the same region
 * @param refStart {int} the rCRS position of the reference's first base
 * @return {list of Difference} the differences, in position order
 * @throws {Error} kind "lineage" on an empty sequence, a non-positive
 *   `refStart`, or a non-nucleotide character
 */
export func align(reference as string, sample as string, refStart as int) {
    def r as string init clean($reference);
    def s as string init clean($sample);
    if (len($r) == 0 or len($s) == 0) {
        fail("both the reference and the sample must be non-empty");
    }
    if (len($r) > MAX_ALIGN_LENGTH or len($s) > MAX_ALIGN_LENGTH) {
        fail("sequence longer than {MAX_ALIGN_LENGTH} bases: alignment is O(nm), so a" +
            " whole-genome region would take far longer than anyone expects it to." +
            " Align the sequenced region (the control region is ~1.1 kb), not the genome");
    }
    if ($refStart < 1) {
        fail("refStart must be a positive reference position");
    }
    def pair as list of string init alignPair($r, $s);
    def gr as string init shiftGapsRight($pair[0], $pair[1]);
    def gs as string init shiftGapsRight($pair[1], $gr);
    return renderDifferences(strings.chars($gr), strings.chars($gs), $refStart);
}

# renderDifferences walks the aligned columns and emits one Difference per
# column that is not a plain match.
func renderDifferences(gr as list of string, gs as list of string, refStart as int) {
    def out as list of Difference init [];
    def pos as int init $refStart - 1;
    def ins as int init 0;
    for (def i in 0..len($gr)) {
        def a as string init $gr[$i];
        def b as string init $gs[$i];
        if ($a == "-") {
            # A gap in the reference is an insertion, numbered after the last
            # reference position seen.
            $ins = $ins + 1;
            $out[] = Difference{position: $pos, insertion: $ins, base: $b, kind: "ins"};
        } else {
            $pos = $pos + 1;
            $ins = 0;
            if ($b == "-") {
                $out[] = Difference{position: $pos, insertion: 0, base: "", kind: "del"};
            } elseif ($a != $b) {
                $out[] = Difference{position: $pos, insertion: 0, base: $b, kind: "sub"};
            }
        }
    }
    return $out;
}

/**
 * An mtDNA haplotype built by aligning a sequence against its reference.
 *
 * The convenience pairing of `align` with `mtHaplotype`: the sequenced range is
 * recorded from `refStart` and the reference length, so the result carries the
 * coverage its differences are meaningful over.
 * @param id {string} a label for the haplotype
 * @param reference {string} the reference sequence for the region
 * @param sample {string} the sequenced sample
 * @param refStart {int} the rCRS position of the reference's first base
 * @return {MtHaplotype} the haplotype
 * @throws {Error} kind "lineage" as for `align`
 */
export func mtHaplotypeFrom(id as string, reference as string, sample as string, refStart as int) {
    def r as string init clean($reference);
    def stop as int init $refStart + len($r) - 1;
    def range as string init "{$refStart}-{$stop}";
    return MtHaplotype{
        id: $id,
        ranges: [$range],
        differences: align($reference, $sample, $refStart)
    };
}

# --- Y-STR haplotypes -------------------------------------------------------

/**
 * A Y-STR haplotype.
 * @param id {string} a label for the haplotype
 * @param loci {map of string to string} locus name -> allele designation(s),
 *   multi-copy loci comma-separated ("11,14")
 * @return {YHaplotype} the haplotype
 * @throws {Error} kind "lineage" on an empty locus name or allele designation
 */
export func yHaplotype(id as string, loci as map of string to string) {
    for (def locus in $loci) {
        if (len($locus) == 0) {
            fail("Y haplotype '" + $id + "' has an empty locus name");
        }
        def parts as list of string init alleleList($loci[$locus]);
        if (len($parts) == 0) {
            fail("Y haplotype '" + $id + "' locus '" + $locus + "' has no alleles");
        }
        for (def a in $parts) {
            repeatCount($a);
        }
    }
    return YHaplotype{id: $id, loci: $loci};
}

# alleleList splits a possibly multi-copy designation into its alleles. Both
# "11,14" and "11/14" are accepted, since databases differ.
func alleleList(text as string) {
    def out as list of string init [];
    for (def part in strings.split(strings.replace($text, "/", ","), ",")) {
        def t as string init strings.trim($part);
        if (len($t) > 0) {
            $out[] = $t;
        }
    }
    return $out;
}

/**
 * The loci a Y haplotype records, in insertion order.
 * @param h {YHaplotype} the haplotype
 * @return {list of string} the locus names
 */
export func yLoci(h as YHaplotype) {
    return maps.keys($h.loci);
}

/**
 * The loci two Y haplotypes have in common, in the first one's order.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @return {list of string} the shared locus names
 */
export func ySharedLoci(a as YHaplotype, b as YHaplotype) {
    def out as list of string init [];
    for (def locus in $a.loci) {
        if (maps.has($b.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * Whether two Y haplotypes agree exactly at every shared locus.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @return {bool} true if every shared locus carries the same allele set
 * @throws {Error} kind "lineage" if the two share no locus
 */
export func yMatches(a as YHaplotype, b as YHaplotype) {
    def shared as list of string init ySharedLoci($a, $b);
    if (len($shared) == 0) {
        fail("Y haplotypes '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    for (def locus in $shared) {
        if (sortedAlleles($a.loci[$locus]) != sortedAlleles($b.loci[$locus])) {
            return false;
        }
    }
    return true;
}

# sortedAlleles normalises a designation to an ascending numeric list, so
# "14,11" and "11,14" compare equal. A non-numeric designation is rejected here
# rather than being left to `convert`, so the caller sees a lineage error naming
# the offending allele instead of a bare toFloat failure.
func sortedAlleles(text as string) {
    def vals as list of float init [];
    for (def a in alleleList($text)) {
        $vals[] = repeatCount($a);
    }
    return lists.sort($vals);
}

# repeatCount parses a Y-STR allele designation as a repeat count.
func repeatCount(allele as string) {
    def dots as int init 0;
    def digits as int init 0;
    for (def ch in strings.chars($allele)) {
        if ($ch == ".") {
            $dots = $dots + 1;
        } elseif (isDigit($ch)) {
            $digits = $digits + 1;
        } else {
            fail("allele '" + $allele + "' is not a repeat count, so it cannot be compared");
        }
    }
    if ($digits == 0 or $dots > 1) {
        fail("allele '" + $allele + "' is not a repeat count, so it cannot be compared");
    }
    return convert.toFloat($allele);
}

/**
 * The stepwise distance between two Y haplotypes.
 *
 * The sum, over shared loci, of the absolute repeat-count differences - the
 * distance a single-step mutation model counts in. Two haplotypes at distance
 * 1 differ by one repeat at one locus.
 *
 * A multi-copy locus is compared after sorting both sides, so `"11,14"` against
 * `"14,11"` is distance 0. Differing copy counts at a shared locus are an error
 * rather than a guess.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @return {float} the total stepwise distance
 * @throws {Error} kind "lineage" if the two share no locus, if a shared locus
 *   has different copy counts, or if an allele is not numeric
 */
export func yStepDistance(a as YHaplotype, b as YHaplotype) {
    def shared as list of string init ySharedLoci($a, $b);
    if (len($shared) == 0) {
        fail("Y haplotypes '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    def total as float init 0.0;
    for (def locus in $shared) {
        def xs as list of float init sortedAlleles($a.loci[$locus]);
        def ys as list of float init sortedAlleles($b.loci[$locus]);
        def na as int init len($xs);
        def nb as int init len($ys);
        if ($na != $nb) {
            fail("locus '" + $locus + "' has {$na} alleles in '" + $a.id + "' but {$nb} in '" +
                $b.id + "'");
        }
        for (def i in 0..len($xs)) {
            $total = $total + math.abs($xs[$i] -$ys[$i]);
        }
    }
    return $total;
}

/**
 * Whether two Y haplotypes are within `steps` repeat units of each other.
 *
 * `yMatchesWithin(a, b, 1.0)` is the single-step comparison used to keep a
 * paternal lineage together across one mutation; `yMatchesWithin(a, b, 0.0)` is
 * `yMatches` over numeric loci.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @param steps {float} the largest total stepwise distance still counted as a match
 * @return {bool} true if the stepwise distance is at most `steps`
 * @throws {Error} kind "lineage" as for `yStepDistance`, or on a negative `steps`
 */
export func yMatchesWithin(a as YHaplotype, b as YHaplotype, steps as float) {
    if ($steps < 0.0) {
        fail("steps must not be negative");
    }
    return yStepDistance($a, $b) <= $steps;
}
