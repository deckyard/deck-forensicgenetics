# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * mtDNA and Y-STR haplotypes and their counting estimates.
 *
 * A lineage haplotype is inherited as one unit, so there is no product rule:
 * the match probability is a count, `k` sightings among `N` database profiles.
 * `lineageMatchProbability` returns the Clopper-Pearson one-sided 95% upper
 * bound on `k/N`, the figure to quote.
 *
 * The caller supplies `k` and `N` from the database searched (EMPOP, YHRD, an
 * in-house set). No reference data is bundled, and no rCRS: `align` takes the
 * reference region as an argument. This module imports nothing from the rest
 * of the deck.
 * @module lineage
 * @see https://empop.online/ - the EMPOP mtDNA database
 * @see https://yhrd.org/ - the Y-STR Haplotype Reference Database
 * @example
 * import "@deckyard/forensicgenetics/lineage.j" as lineage;
 *
 * # A haplotype seen 2 times in a database of 5000.
 * def p as float init lineage.lineageMatchProbability(2, 5000);   # ~0.00126
 * def lr as float init lineage.lineageLikelihoodRatio(2, 5000);   # ~795
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
use regex;
use decimal;

/**
 * The one-sided confidence level `lineageMatchProbability` uses.
 */
export def const DEFAULT_CONFIDENCE as float init 0.95;

# The largest reference position a difference designation may name. The human
# mitochondrial genome is 16569 bases, so this is far above anything real; the
# bound exists because positions are sorted through a fixed-width key, which a
# wider number would order wrongly.
def const MAX_POSITION as int init 999999;

/**
 * One difference from the reference, in rCRS nomenclature.
 *
 * `16189C` is `{position: 16189, insertion: 0, base: "C", kind: "sub"}`,
 * `315.1C` is `{position: 315, insertion: 1, base: "C", kind: "ins"}`, and
 * `249d` is `{position: 249, insertion: 0, base: "", kind: "del"}`.
 * @field position {int} the 1-based reference position; 0 for an insertion
 *   before the first base
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
 * An mtDNA haplotype: the sequenced ranges and the differences from the
 * reference within them. `mtMatches` ignores the ranges; compare them with
 * `mtRangesAgree`.
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
 * A multi-copy locus (DYS385a/b, DYS389I/II reported together) holds its
 * alleles comma- or slash-separated, `"11,14"` or `"11/14"`.
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
 * The point estimate of a haplotype's frequency, `k / N`.
 *
 * Not the figure to report: it is 0.0 at `k == 0`. Quote
 * `lineageMatchProbability`.
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
 * The Clopper-Pearson exact-binomial one-sided upper bound on a haplotype
 * frequency.
 *
 * At `k == 0` it is `1 - alpha^(1/N)`, with `alpha = 1 - confidence` (the
 * "rule of three" at 95%); at `k == N` it is exactly 1.0. Computed as the upper
 * end of the two-sided Clopper-Pearson interval at level `2c - 1`.
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
 * The lineage match probability to report: `haplotypeUpperBound` at
 * `DEFAULT_CONFIDENCE` (0.95).
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @return {float} the conservative match probability
 * @throws {Error} kind "lineage" on a non-positive N or a k outside [0, N]
 */
export func lineageMatchProbability(k as int, n as int) {
    return haplotypeUpperBound($k, $n, DEFAULT_CONFIDENCE);
}

/**
 * The lineage likelihood ratio, `1 / lineageMatchProbability`.
 *
 * "The sample came from this maternal (mtDNA) or paternal (Y) lineage" against
 * "it came from an unrelated lineage". It is 1.0 at `k == N`.
 * @param k {int} how many times the haplotype was seen in the database
 * @param n {int} the database size
 * @return {float} the likelihood ratio
 * @throws {Error} kind "lineage" on a non-positive N or a k outside [0, N]
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
# heteroplasmic position still parses. A designation ending in D is read as a
# deletion before isBase is consulted, so the code D (A/G/T) never reaches it.
func isDigit(ch as string) {
    return strings.indexOf("0123456789", $ch) >= 0;
}

func isBase(ch as string) {
    return strings.indexOf("ACGTURYSWKMBDHVN", $ch) >= 0;
}

# A whole designation: optional reference base, position, optional insertion
# index, then one IUPAC base or a deletion mark. The ranges are checked in code.
def const DESIGNATION as string init '^[ACGTURYSWKMBDHVN]?([0-9]{1,6})' +
    '(?:\.([0-9]{1,4}))?([ACGTURYSWKMBDHVN]|DEL|-)$';

# quickDifference parses a well-formed, upper-cased designation with one
# anchored match. Its result has kind "" when the designation needs the full
# scan in `parseDifference`: malformed, out of range, or an error to report.
func quickDifference(t as string) {
    def none as Difference init Difference{position: 0, insertion: 0, base: "", kind: ""};
    def m as regex.Match init regex.find(DESIGNATION, $t);
    if ($m.start != 0) {
        return $none;
    }
    def position as int init convert.toInt($m.groups[0]);
    def insertion as int init 0;
    if (len($m.groups[1]) > 0) {
        $insertion = convert.toInt($m.groups[1]);
        if ($insertion < 1) {
            return $none;
        }
    }
    def rest as string init $m.groups[2];
    if ($position < 1 and $insertion == 0) {
        return $none;
    }
    if ($rest == "D" or $rest == "DEL" or $rest == "-") {
        if ($insertion > 0) {
            return $none;
        }
        return Difference{position: $position, insertion: 0, base: "", kind: "del"};
    }
    def kind as string init "sub";
    if ($insertion > 0) {
        $kind = "ins";
    }
    return Difference{position: $position, insertion: $insertion, base: $rest, kind: $kind};
}

/**
 * Parse one difference designation in rCRS nomenclature.
 *
 * Accepts `16189C`, `A263G` (the reference base is dropped), `315.1C`, and the
 * deletions `249d`, `249D`, `249DEL` and `249-`. Case-insensitive; surrounding
 * whitespace is ignored. The base may be any IUPAC code except `D`, which is
 * always read as a deletion (`249D`), never as the ambiguity code A/G/T. The
 * position must be in [1, 999999] and an insertion index in [1, 9999].
 * Position 0 is accepted for an insertion only: `0.1A` is a base inserted
 * before position 1, the form `align` emits for an insertion before the first
 * reference base.
 * @param text {string} the designation
 * @return {Difference} the parsed difference
 * @throws {Error} kind "lineage" if the designation does not parse
 */
export func parseDifference(text as string) {
    def t as string init strings.upper(strings.trim($text));
    if (len($t) == 0) {
        fail("empty difference designation");
    }
    # Most designations are well formed; anything the fast path does not take
    # goes through the scan below, which owns every error message.
    def quick as Difference init quickDifference($t);
    if (len($quick.kind) > 0) {
        return $quick;
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
    if ($position > MAX_POSITION) {
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
    if ($position == 0 and $insertion == 0) {
        fail("difference '" + $text + "': position must be in [1, {MAX_POSITION}]; 0 is" +
            " only valid for an insertion before the first base (0.1A)");
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
 * A difference in rCRS nomenclature.
 * @param d {Difference} the difference to render
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
 * An mtDNA haplotype from whitespace-separated designations.
 * @param id {string} a label for the haplotype
 * @param ranges {list of string} the sequenced regions, e.g. ["16024-16365"]
 * @param text {string} the differences, e.g. "16189C 16223T 315.1C"
 * @return {MtHaplotype} the haplotype, sorted by position then insertion index
 * @throws {Error} kind "lineage" if a designation does not parse, or if two
 *   designations name the same position and insertion index (a repeat such as
 *   `263G 263G` or a conflict such as `263G 263A`)
 */
export func mtHaplotype(id as string, ranges as list of string, text as string) {
    def flat as string init strings.replace(strings.replace($text, "\t", " "), "\n", " ");
    def out as list of Difference init [];
    for (def token in strings.split($flat, " ")) {
        if (len(strings.trim($token)) > 0) {
            $out[] = parseDifference($token);
        }
    }
    def sorted as list of Difference init sortDifferences($out);
    # After sorting, a repeated or conflicting designation (263G 263G, 263G
    # 263A, 263G 263d) sits next to its twin. Either would count twice in
    # mtDistance, so it is an input error.
    for (def i in 0..len($sorted)) {
        if ($i > 0 and sortKey($sorted[$i]) == sortKey($sorted[$i - 1])) {
            fail("haplotype '" + $id + "': '" + formatDifference($sorted[$i - 1]) + "' and '" +
                formatDifference($sorted[$i]) + "' name the same position");
        }
    }
    return MtHaplotype{id: $id, ranges: $ranges, differences: $sorted};
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

# sortDifferences orders differences by position then insertion index, stable.
# It sorts fixed-width string keys with lists.sort, not the structs with
# lists.sortBy: a func value anywhere in a call tree makes the interpreter
# deep-copy compound arguments up the chain (a caller's whole mtDNA tree in
# mitoleaf).
func sortDifferences(ds as list of Difference) {
    def keys as list of string init [];
    for (def i in 0..len($ds)) {
        $keys[] = sortKey($ds[$i]) + pad(convert.toString($i), 7);
    }
    $keys = lists.sort($keys);
    def out as list of Difference init [];
    for (def k in $keys) {
        $out[] = $ds[convert.toInt($k[10..len($k)])];
    }
    return $out;
}

/**
 * The differences of a haplotype as designation strings.
 * @param h {MtHaplotype} the haplotype to render
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
 * Whether two mtDNA haplotypes carry the same differences.
 *
 * Compares differences only, not sequenced ranges: check `mtRangesAgree` too.
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {bool} true if the difference sets are identical
 */
export func mtMatches(a as MtHaplotype, b as MtHaplotype) {
    return mtDistance($a, $b) == 0;
}

/**
 * The size of the symmetric difference of two haplotypes' designation sets.
 *
 * A designation only one carries counts once; a position both carry with
 * different bases counts twice.
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {int} the number of differing designations
 */
export func mtDistance(a as MtHaplotype, b as MtHaplotype) {
    def n as int init mergeDistance($a, $b);
    if ($n >= 0) {
        return $n;
    }
    return setDistance($a, $b);
}

# mergeDistance is mtDistance by one merge over the two difference lists,
# ordered by (position, insertion) as mtHaplotype leaves them. Designations
# differ whenever (position, insertion) does; at an equal key the same kind
# and base give the same designation, anything else is rendered and compared.
# It returns -1 when either list is not strictly ascending or holds an unknown
# kind (a struct literal); setDistance handles those. Keys are compared inline:
# a helper call per comparison cost more than the comparison.
func mergeDistance(a as MtHaplotype, b as MtHaplotype) {
    if (not mergeable($a) or not mergeable($b)) {
        return -1;
    }
    def na as int init len($a.differences);
    def nb as int init len($b.differences);
    def i as int init 0;
    def j as int init 0;
    def n as int init 0;
    while ($i < $na and $j < $nb) {
        def c as int init $a.differences[$i].position - $b.differences[$j].position;
        if ($c == 0) {
            $c = $a.differences[$i].insertion - $b.differences[$j].insertion;
        }
        if ($c < 0) {
            $n = $n + 1;
            $i = $i + 1;
        } elseif ($c > 0) {
            $n = $n + 1;
            $j = $j + 1;
        } else {
            if (not sameDifference($a.differences[$i], $b.differences[$j])) {
                $n = $n + 2;
            }
            $i = $i + 1;
            $j = $j + 1;
        }
    }
    return $n + ($na - $i) + ($nb - $j);
}

# mergeable reports whether every difference has a known kind and a larger
# (position, insertion) than the one before.
func mergeable(h as MtHaplotype) {
    def lastPos as int init -1;
    def lastIns as int init 0;
    for (def i in 0..len($h.differences)) {
        def k as string init $h.differences[$i].kind;
        def p as int init $h.differences[$i].position;
        def q as int init $h.differences[$i].insertion;
        if ($k != "sub" and $k != "ins" and $k != "del") {
            return false;
        }
        if ($p < $lastPos or ($p == $lastPos and $q <= $lastIns)) {
            return false;
        }
        $lastPos = $p;
        $lastIns = $q;
    }
    return true;
}

# sameDifference reports whether two differences at one key render alike.
func sameDifference(x as Difference, y as Difference) {
    if ($x.kind == $y.kind and $x.base == $y.base) {
        return true;
    }
    return formatDifference($x) == formatDifference($y);
}

# setDistance is the symmetric difference of the designation sets.
func setDistance(a as MtHaplotype, b as MtHaplotype) {
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
 * Duplicates are ignored: the ranges compare as sets.
 * @param a {MtHaplotype} one haplotype
 * @param b {MtHaplotype} the other
 * @return {bool} true if each declares exactly the ranges the other does
 */
export func mtRangesAgree(a as MtHaplotype, b as MtHaplotype) {
    # Compared as sets: a range listed twice is still one range, and a count
    # check plus membership would let ["x", "x"] agree with ["x", "y"].
    return lists.sort(maps.keys(setOf($a.ranges))) == lists.sort(maps.keys(setOf($b.ranges)));
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
 * The longest sequence `align` accepts, in bases.
 *
 * Banded alignment of near-identical sequences is roughly linear (about 0.7 s
 * at this cap with a dozen differences). The cap makes a whole-genome call
 * fail at once instead of appearing to hang. The control region is about
 * 1.1 kb.
 */
export def const MAX_ALIGN_LENGTH as int init 5000;

/**
 * The most matrix cells one banded pass of `align` may fill.
 *
 * Divergent sequences widen the band towards the full O(nm) matrix, at about
 * 1.5 us per cell; the doubling passes before the last add at most as much
 * again, so the cap holds a call to about 4 s. It admits a full-width
 * alignment of two 1100-base sequences (a whole control region, however
 * divergent), and a band of half-width 124 at `MAX_ALIGN_LENGTH`. A control
 * region against the rCRS with a few dozen differences needs a band of 16 to
 * 32 (under 75000 cells).
 */
export def const MAX_ALIGN_CELLS as int init 1250000;

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
# holds it; the width doubles until the band's score proves no path outside it
# can do as well (see `alignBand`).
def const INITIAL_BAND as int init 16;

# Band is one banded pass: the two gapped strings, their score, and whether the
# band was wide enough for the result to be the true global optimum.
def struct Band {
    exact as bool,
    score as int,
    top as string,
    bottom as string
};

# alignPair returns the two gapped strings as [gappedRef, gappedSample].
#
# Full Needleman-Wunsch is O(nm), which for a control region is a million cells
# of interpreted arithmetic. Instead the matrix is computed only within a band
# around the diagonal, which is where the optimal path lies when the sequences
# are near-identical - as an mtDNA haplotype and its reference are. A band is
# accepted only when its score beats every path that leaves it (`bandBound`);
# otherwise the width doubles and the pass repeats. At full width the band is
# the whole matrix, so the result is always the true optimum, never an
# approximation.
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
    checkBandCells(len($r), len($s), $w);
    def band as Band init alignBand($rb, $sb, $r, $s, $w);
    while (not $band.exact and $w < $limit) {
        $w = $w * 2;
        # A band at least as wide as a row costs as much as the full matrix,
        # so go straight to full width rather than pay for a pass that might
        # still not prove the optimum.
        if ($w > $limit or 2 * $w + 1 >= len($s) + 1) {
            $w = $limit;
        }
        checkBandCells(len($r), len($s), $w);
        $band = alignBand($rb, $sb, $r, $s, $w);
    }
    return [$band.top, $band.bottom];
}

# checkBandCells refuses a pass of half-width `w` that would fill more than
# MAX_ALIGN_CELLS cells: m rows of 2w + 1 cells, at most the full matrix.
func checkBandCells(m as int, n as int, w as int) {
    def cells as int init math.min($m * (2 * $w + 1), $m * ($n + 1));
    if ($cells > MAX_ALIGN_CELLS) {
        fail("sequences too divergent to align: proving the optimum needs a band of" +
            " {$cells} cells, over MAX_ALIGN_CELLS ({MAX_ALIGN_CELLS}). Align a shorter" +
            " region, or check that the sample belongs to this reference region");
    }
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
        return Band{exact: false, score: NEG, top: "", bottom: ""};
    }
    # Two score rows are live at a time, held in one list and alternated by
    # index: $p is the row just finished, $q the one being filled. Rotating
    # them by assignment instead - prev/cur/spare - deep-copies a full-width
    # row on every pass, which is O(n) work per row where the band needs only
    # O(w), and quadratic across the whole alignment.
    # The two live score rows share one flat list: row k starts at k * width.
    # $p is the row just finished, $q the one being filled.
    def width as int init $n + 2;
    def rows as list of int init [];
    for (def j in 0..(2 * $width)) {
        $rows[] = NEG;
    }
    def p as int init 0;
    def dirs as list of bytes init [];
    def los as list of int init [];
    def firstHi as int init math.min($n, $w);
    def firstRow as bytes;
    $firstRow[] = DIR_DIAG;
    $rows[0] = 0;
    for (def j in 1..($firstHi + 1)) {
        $rows[$j] = $j * GAP;
        $firstRow[] = DIR_LEFT;
    }
    $rows[$firstHi + 1] = NEG;
    $dirs[] = $firstRow;
    $los[] = 0;
    for (def i in 1..($m + 1)) {
        def q as int init 1 - $p;
        def po as int init $p * $width;
        def qo as int init $q * $width;
        def lo as int init math.max(0, $i - $w);
        def hi as int init math.min($n, $i + $w);
        def base as int init $rb[$i - 1];
        def moves as bytes;
        def left as int init NEG;
        if ($lo > 0) {
            # Guard the cell just left of the band: the next row reads it as its
            # diagonal, and it must not see a value two rows stale.
            $rows[$qo + $lo - 1] = NEG;
        } else {
            $rows[$qo] = $i * GAP;
            $moves[] = DIR_UP;
            $left = $i * GAP;
        }
        def from as int init math.max(1, $lo);
        # Each cell's `up` is the next cell's `diag`, so it is read once.
        def diag as int init $rows[$po + $from - 1];
        for (def j in $from..($hi + 1)) {
            def up as int init $rows[$po + $j];
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
            $rows[$qo + $j] = $best;
            $left = $best;
            $diag = $up;
        }
        $rows[$qo + $hi + 1] = NEG;
        $dirs[] = $moves;
        $los[] = $lo;
        $p = $q;
    }
    def score as int init $rows[$p * $width + $n];
    # A traceback that merely stays off the band edge proves nothing: a path
    # far outside the band can still score higher. The band is exact when its
    # score beats the best any path leaving it could reach; a tie widens too,
    # so the traceback equals the full-width one.
    def exact as bool init $w >= math.max($m, $n) or $score > bandBound($m, $n, $w);
    return traceback($dirs, $los, $r, $s, $exact, $score);
}

# bandBound is the highest score any alignment path leaving the band of
# half-width `w` can reach (the Fickett / Ukkonen bound). Such a path visits
# diagonal i - j = w + 1 or -(w + 1) and ends on diagonal m - n, so it has at
# least g = 2(w + 1) - |m - n| gap moves, and with them at most (m + n - g) / 2
# diagonal moves, each worth at most MATCH. Its score is therefore at most
# floor((m + n - g) / 2) * MATCH + g * GAP.
func bandBound(m as int, n as int, w as int) {
    def g as int init 2 * ($w + 1) - math.abs($m - $n);
    return (($m + $n - $g) // 2) * MATCH + $g * GAP;
}

# traceback walks the stored moves back to the origin and packs the result.
func traceback(
    dirs as list of bytes,
    los as list of int,
    r as list of string,
    s as list of string,
    exact as bool,
    score as int) {
    def i as int init len($r);
    def j as int init len($s);
    def gappedRef as list of string init [];
    def gappedSample as list of string init [];
    while ($i > 0 or $j > 0) {
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
        score: $score,
        top: strings.join(lists.reverse($gappedRef), ""),
        bottom: strings.join(lists.reverse($gappedSample), "")
    };
}

# shiftGapsRight slides every gap run in `a` as far 3' as it can go without
# lowering the score, which is what puts an insertion into a homopolymer or
# tandem-repeat run at the run's 3' end - the convention that makes an extra C
# in the 311-315 stretch read as 315.1C rather than 311.1C, and an extra CA in
# a CA repeat land after the last unit.
#
# A maximal run of k gaps starting at column c moves one column right when the
# base after it equals `b[c]`: that base moves to column c, where it matches,
# and the run now covers c+1 .. c+k. A run that reaches the next run merges
# with it. One left-to-right pass does every move: the output is built in a
# list, the pending run is only a count, and each column is written once, so
# the cost is linear in the alignment length (a swap-until-stable loop was
# O(gaps * distance moved): 22 s for 5000 bases against one).
func shiftGapsRight(a as string, b as string) {
    def x as list of string init strings.chars($a);
    def y as list of string init strings.chars($b);
    def out as list of string init [];
    def run as int init 0;
    for (def c in $x) {
        if ($c == "-") {
            $run = $run + 1;
        } elseif ($run > 0 and $c == $y[len($out)]) {
            # The run's first column is len($out): the base takes it and the
            # run moves one column right.
            $out[] = $c;
        } else {
            for (def k in 0..$run) {
                $out[] = "-";
            }
            $run = 0;
            $out[] = $c;
        }
    }
    for (def k in 0..$run) {
        $out[] = "-";
    }
    return strings.join($out, "");
}

/**
 * The rCRS differences of a sequenced region against its reference.
 *
 * A global Needleman-Wunsch alignment (match +1, mismatch -1, gap -2; ties
 * prefer the diagonal), with each indel run slid to the 3' end of any
 * homopolymer or repeat run, so an extra base in the 311-315 C-stretch reads
 * as `315.1C`. Whitespace in either sequence is ignored and bases are
 * upper-cased. The alignment is computed in a band around the diagonal that
 * widens until no path outside it can score as well, so the result is the
 * full-matrix optimum. An insertion before the reference's first base is
 * reported at position `refStart - 1`, which is `0.1A` when `refStart` is 1
 * (`parseDifference` reads it back).
 *
 * Pass the rCRS (NC_012920.1) region you sequenced, with `refStart` the 1-based
 * position of its first base. EMPOP's phylogenetic alignment rules are not
 * applied, so a complex indel may differ from the EMPOP designation.
 * @param reference {string} the reference sequence for the region
 * @param sample {string} the sequenced sample over the same region
 * @param refStart {int} the rCRS position of the reference's first base
 * @return {list of Difference} the differences, in position order
 * @throws {Error} kind "lineage" on an empty sequence, a non-positive
 *   `refStart`, a `refStart` that puts the reference's last base past position
 *   999999, a non-IUPAC character, a sequence longer than
 *   `MAX_ALIGN_LENGTH`, or sequences so divergent that a band pass would
 *   exceed `MAX_ALIGN_CELLS`
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
    # The last reference base must have a designation parseDifference accepts.
    if ($refStart > MAX_POSITION - len($r) + 1) {
        fail("refStart {$refStart} puts the reference's last base past position" +
            " {MAX_POSITION}, the largest a designation may name");
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
 * An mtDNA haplotype from `align`, with the range
 * `"<refStart>-<refStart + reference length - 1>"`.
 * @param id {string} a label for the haplotype
 * @param reference {string} the reference sequence for the region
 * @param sample {string} the sequenced sample
 * @param refStart {int} the rCRS position of the reference's first base
 * @return {MtHaplotype} the haplotype
 * @throws {Error} kind "lineage" as for `align`
 */
export func mtHaplotypeFrom(id as string, reference as string, sample as string, refStart as int) {
    # align validates refStart before the range end is computed from it
    def diffs as list of Difference init align($reference, $sample, $refStart);
    def stop as int init $refStart + len(clean($reference)) - 1;
    return MtHaplotype{id: $id, ranges: ["{$refStart}-{$stop}"], differences: $diffs};
}

# --- Y-STR haplotypes -------------------------------------------------------

/**
 * A Y-STR haplotype.
 * @param id {string} a label for the haplotype
 * @param loci {map of string to string} locus name -> allele designation(s),
 *   multi-copy loci as "11,14" or "11/14"
 * @return {YHaplotype} the haplotype
 * @throws {Error} kind "lineage" on an empty locus name, a locus with no allele,
 *   or an allele that is not a repeat count
 */
export func yHaplotype(id as string, loci as map of string to string) {
    for (def locus in $loci) {
        if (len($locus) == 0) {
            fail("Y haplotype '" + $id + "' has an empty locus name");
        }
        if (not multiCopy($loci[$locus]) and len(strings.trim($loci[$locus])) > 0) {
            repeatCount(strings.trim($loci[$locus]));
            continue;
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
 * @param h {YHaplotype} the haplotype to read
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
 * Whether two Y haplotypes agree exactly at every shared locus. Multi-copy
 * loci compare as sorted lists, so "11,14" equals "14,11".
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @return {bool} true if every shared locus carries the same alleles
 * @throws {Error} kind "lineage" if the two share no locus, or an allele is
 *   not a repeat count
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
    if (not multiCopy($text)) {
        def t as string init strings.trim($text);
        if (len($t) == 0) {
            return [];
        }
        return [repeatCount($t)];
    }
    def vals as list of float init [];
    for (def a in alleleList($text)) {
        $vals[] = repeatCount($a);
    }
    return lists.sort($vals);
}

# MAX_ALLELE_DIGITS bounds the digits on each side of a repeat count's dot;
# real counts have two or three, and a longer string would overflow toFloat.
def const MAX_ALLELE_DIGITS as int init 6;

# REPEAT_COUNT is a repeat count: digits, optionally one '.' with digits on
# both sides ("11", "17.2"; not "11." or ".5").
def const REPEAT_COUNT as string init '^[0-9]{1,6}(\.[0-9]{1,6})?$';

# repeatCount parses a Y-STR allele designation as a repeat count.
func repeatCount(allele as string) {
    if (not regex.matches(REPEAT_COUNT, $allele)) {
        fail("allele '" + $allele + "' is not a repeat count (at most {MAX_ALLELE_DIGITS}" +
            " digits, optionally '.' and at most {MAX_ALLELE_DIGITS} more), so it cannot" +
            " be compared");
    }
    return convert.toFloat($allele);
}

# multiCopy reports whether a designation lists more than one allele.
func multiCopy(text as string) {
    return strings.contains($text, ",") or strings.contains($text, "/");
}

# byRepeatCount sorts allele labels by repeat count, stable (insertion sort
# over the parsed values: a locus has a handful of copies). No lists.sortBy:
# see sortDifferences.
func byRepeatCount(xs as list of string) {
    def vals as list of float init [];
    def out as list of string init [];
    for (def x in $xs) {
        def v as float init repeatCount($x);
        def j as int init len($out);
        $vals[] = $v;
        $out[] = $x;
        while ($j > 0 and $vals[$j - 1] > $v) {
            $vals[$j] = $vals[$j - 1];
            $out[$j] = $out[$j - 1];
            $j = $j - 1;
        }
        $vals[$j] = $v;
        $out[$j] = $x;
    }
    return $out;
}

# pairSteps is |x - y| for one pair of alleles, exact in decimal; -2.0 when
# it is not a whole number of repeats.
func pairSteps(x as string, y as string) {
    if ($x == $y) {
        return 0.0;
    }
    def d as decimal.Decimal init decimal.abs(decimal.sub($x, $y));
    if (not decimal.isInteger($d)) {
        return -2.0;
    }
    return decimal.toFloat($d);
}

# locusSteps is the stepwise distance at one locus, exact in decimal: the
# alleles of each side are sorted by repeat count and paired, and each pair
# contributes |x - y|. It returns -1.0 when the copy counts differ and -2.0
# when a pair differs by a non-integer (17.2 against 18): a microvariant is
# reached from another allele only through alleles with the same fractional
# part, so such a pair has no stepwise distance.
func locusSteps(x as string, y as string) {
    # Single-copy on both sides (most loci): no split, no sort.
    if (not multiCopy($x) and not multiCopy($y)) {
        def tx as string init strings.trim($x);
        def ty as string init strings.trim($y);
        if (len($tx) > 0 and len($ty) > 0) {
            def vx as float init repeatCount($tx);
            def vy as float init repeatCount($ty);
            # whole counts of at most 6 digits are exact in float
            if (not strings.contains($tx, ".") and not strings.contains($ty, ".")) {
                return math.abs($vx - $vy);
            }
            return pairSteps($tx, $ty);
        }
    }
    def xs as list of string init byRepeatCount(alleleList($x));
    def ys as list of string init byRepeatCount(alleleList($y));
    if (len($xs) != len($ys)) {
        return -1.0;
    }
    # Each pair is a whole number of repeats, exact in float, so the sum is.
    def total as float init 0.0;
    for (def i in 0..len($xs)) {
        def d as float init pairSteps($xs[$i], $ys[$i]);
        if ($d < 0.0) {
            return -2.0;
        }
        $total = $total + $d;
    }
    return $total;
}

/**
 * The stepwise distance between two Y haplotypes: the sum over shared loci of
 * the absolute repeat-count differences, computed exactly in decimal.
 * Multi-copy loci are sorted before comparing, so `"11,14"` against `"14,11"`
 * is 0. Two alleles with different fractional parts (`"17.2"` against `"18"`)
 * are not stepwise neighbours and throw; `"17.2"` against `"18.2"` is 1 step.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @return {float} the total stepwise distance, a whole number
 * @throws {Error} kind "lineage" if the two share no locus, if a shared locus
 *   has different copy counts, if paired alleles differ in their fractional
 *   part, or if an allele is not a repeat count
 */
export func yStepDistance(a as YHaplotype, b as YHaplotype) {
    def shared as list of string init ySharedLoci($a, $b);
    if (len($shared) == 0) {
        fail("Y haplotypes '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    def total as float init 0.0;
    for (def locus in $shared) {
        def x as string init $a.loci[$locus];
        def y as string init $b.loci[$locus];
        def d as float init locusSteps($x, $y);
        if ($d == -1.0) {
            def na as int init len(alleleList($x));
            def nb as int init len(alleleList($y));
            fail("locus '" + $locus + "' has {$na} alleles in '" + $a.id + "' but {$nb} in '" +
                $b.id + "'");
        }
        if ($d == -2.0) {
            fail("locus '" + $locus + "': '" + $x + "' in '" + $a.id + "' and '" + $y + "' in '" +
                $b.id + "' differ in their fractional part, so no step count joins them");
        }
        $total = $total + $d;
    }
    return $total;
}

/**
 * Whether `yStepDistance` is at most `steps`.
 *
 * `steps = 1.0` keeps a paternal lineage together across one single-step
 * mutation; `steps = 0.0` is `yMatches`. A shared locus where `yStepDistance`
 * throws - different copy counts, or alleles that differ in their fractional
 * part - makes the result false, as it does for `yMatches`.
 * @param a {YHaplotype} one haplotype
 * @param b {YHaplotype} the other
 * @param steps {float} the largest distance still counted as a match
 * @return {bool} true if the stepwise distance is at most `steps`
 * @throws {Error} kind "lineage" if the two share no locus, if an allele is
 *   not a repeat count, or on a negative `steps`
 */
export func yMatchesWithin(a as YHaplotype, b as YHaplotype, steps as float) {
    if ($steps < 0.0) {
        fail("steps must not be negative");
    }
    def shared as list of string init ySharedLoci($a, $b);
    if (len($shared) == 0) {
        fail("Y haplotypes '" + $a.id + "' and '" + $b.id + "' share no locus");
    }
    def total as float init 0.0;
    for (def locus in $shared) {
        def d as float init locusSteps($a.loci[$locus], $b.loci[$locus]);
        if ($d < 0.0) {
            return false;
        }
        $total = $total + $d;
    }
    return $total <= $steps;
}
