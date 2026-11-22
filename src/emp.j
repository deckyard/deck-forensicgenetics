# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.27.0

/**
 * EMPOP's EMP files: mtDNA haplotypes as text, read and written.
 *
 * An EMP file is the tab-delimited format EMPOP takes for submission and
 * search: `#` starts a comment, a `#!` line names the sequenced range(s) for
 * every haplotype below it (`#! 16024-576`, `#! 16024-16365 73-340`), and each
 * haplotype is one line of sample ID, haplogroup (`?` if unknown), count (0 to
 * 9999; 0 excludes it from analysis) and its differences from the rCRS, one
 * per cell (`16189C`, `315.1C`, `249DEL`).
 *
 * The module takes and returns text and does no I/O: read the file with
 * `fs.readString` and pass it in. The haplotypes are `lineage.MtHaplotype`
 * values, so `lineage`'s comparisons and counting estimates apply to them.
 * @module emp
 * @see https://empop-use.readthedocs.io/en/latest/tools/ - the EMP file format
 * @see https://empop.online/ - the EMPOP mtDNA database
 * @example
 * import "@deckyard/forensicgenetics/emp.j" as emp;
 * import "@deckyard/forensicgenetics/lineage.j" as lineage;
 *
 * def db as emp.EmpFile init emp.readEmp(fs.readString("population.emp"));
 * def h as lineage.MtHaplotype init lineage.mtHaplotype("evidence", ["16024-576"],
 *     "16189C 16223T 263G 315.1C");
 * def k as int init emp.empCount($db, $h);
 * def p as float init lineage.lineageMatchProbability($k, emp.empSampleSize($db));
 */
use lists;
use maps;
use regex;
use strings;
use convert;

import "./lineage.j" as lineage;

/**
 * The length of the rCRS (NC_012920.1), the largest position a range or a
 * difference in an EMP file may name.
 */
export def const GENOME_LENGTH as int init 16569;

/** The largest count a haplotype line may carry. */
export def const MAX_COUNT as int init 9999;

/**
 * The most haplotypes `readEmp` takes from one file. EMPOP's whole database
 * holds a few tens of thousands; the cap stops a runaway input early.
 */
export def const MAX_ENTRIES as int init 200000;

/**
 * One haplotype line of an EMP file.
 * @field haplotype {lineage.MtHaplotype} the sample ID as `id`, the ranges of
 *   the `#!` line above it, and its differences from the rCRS
 * @field haplogroup {string} the haplogroup as written; "?" if unknown
 * @field count {int} how often the haplotype was seen, 0 to `MAX_COUNT`; 0
 *   excludes it from analysis
 * @field line {int} the 1-based line it was read from; 0 for an entry built
 *   with `empEntry`
 */
export def struct EmpEntry {
    haplotype as lineage.MtHaplotype,
    haplogroup as string,
    count as int,
    line as int
};

/**
 * An EMP file.
 * @field comments {list of string} the comment lines above the first `#!`
 *   line (title, contact, population), without the `#` and one following
 *   space
 * @field entries {list of EmpEntry} the haplotype lines, in file order
 */
export def struct EmpFile {
    comments as list of string,
    entries as list of EmpEntry
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "emp", message: $msg, file: "", line: 0, col: 0};
}

/**
 * Throw a reached cap: the generic kind "limit".
 *
 * @param msg {string} what is over the cap
 * @throws {Error} kind "limit", always
 */
func failLimit(msg as string) {
    throw Error{kind: "limit", message: $msg, file: "", line: 0, col: 0};
}

# --- checks -----------------------------------------------------------------

# EMPOP allows letters, digits, "-", "_" and "/" in a sample ID: no blanks,
# no umlauts, nothing a later tool might split on.
def const ID_PATTERN as string init '^[A-Za-z0-9_/-]+$';
def const COUNT_PATTERN as string init '^[0-9]{1,4}$';
def const RANGE_PATTERN as string init '^([0-9]{1,5})-([0-9]{1,5})$';
def const SPACE_PATTERN as string init '\s';

# A stretch of the genome a range covers, lo <= hi. A range that runs past
# 16569 back to 1 (16024-576, across the origin of the circular genome) is two
# segments.
def struct Segment {
    lo as int,
    hi as int
};

# rangeProblem says what is wrong with one sequence range, or "" if nothing.
func rangeProblem(r as string) {
    def m as regex.Match init regex.find(RANGE_PATTERN, $r);
    if ($m.start != 0) {
        return "sequence range '" + $r + "' is not start-end, as in 16024-576";
    }
    def lo as int init convert.toInt($m.groups[0]);
    def hi as int init convert.toInt($m.groups[1]);
    if ($lo < 1 or $lo > GENOME_LENGTH or $hi < 1 or $hi > GENOME_LENGTH) {
        return "sequence range '" + $r + "' leaves 1-{GENOME_LENGTH}";
    }
    return "";
}

# segmentsOf turns valid ranges into segments, in range order.
func segmentsOf(ranges as list of string) {
    def out as list of Segment init [];
    for (def r in $ranges) {
        def m as regex.Match init regex.find(RANGE_PATTERN, $r);
        def lo as int init convert.toInt($m.groups[0]);
        def hi as int init convert.toInt($m.groups[1]);
        if ($lo <= $hi) {
            $out[] = Segment{lo: $lo, hi: $hi};
        } else {
            $out[] = Segment{lo: $lo, hi: GENOME_LENGTH};
            $out[] = Segment{lo: 1, hi: $hi};
        }
    }
    return $out;
}

# segmentIndex is the first segment holding a position, or -1.
func segmentIndex(segs as list of Segment, position as int) {
    for (def i in 0..len($segs)) {
        if ($position >= $segs[$i].lo and $position <= $segs[$i].hi) {
            return $i;
        }
    }
    return -1;
}

/**
 * A difference as an EMP file writes it: `16189C`, `315.1C`, `249DEL`.
 *
 * @param d {lineage.Difference} the difference
 * @return {string} its designation
 * @throws {Error} kind "emp" on an unknown `kind`
 */
export func empDesignation(d as lineage.Difference) {
    match ($d.kind) {
        when "sub" { return "{$d.position}{$d.base}"; }
        when "ins" { return "{$d.position}.{$d.insertion}{$d.base}"; }
        when "del" { return "{$d.position}DEL"; }
        else { fail("unknown difference kind '" + $d.kind + "'"); }
    }
}

# haplotypeProblem checks what EMPOP checks beyond the notation: the ranges,
# every position within 1-16569 and within a range, and insertions numbered
# from .1 without a gap. The differences must be sorted (mtHaplotype sorts).
func haplotypeProblem(h as lineage.MtHaplotype) {
    if (len($h.ranges) == 0) {
        return "no sequence range";
    }
    for (def r in $h.ranges) {
        def problem as string init rangeProblem($r);
        if (len($problem) > 0) {
            return $problem;
        }
    }
    def segs as list of Segment init segmentsOf($h.ranges);
    def ds as list of lineage.Difference init $h.differences;
    for (def i in 0..len($ds)) {
        def d as lineage.Difference init $ds[$i];
        if ($d.position < 1 or $d.position > GENOME_LENGTH) {
            return "'" + empDesignation($d) + "' leaves 1-{GENOME_LENGTH}";
        }
        if (segmentIndex($segs, $d.position) < 0) {
            return "'" + empDesignation($d) + "' lies outside the sequence range " +
                strings.join($h.ranges, " ");
        }
        if ($d.kind == "ins") {
            def expected as int init 1;
            if ($i > 0 and $ds[$i - 1].kind == "ins" and $ds[$i - 1].position == $d.position) {
                $expected = $ds[$i - 1].insertion + 1;
            }
            if ($d.insertion != $expected) {
                return "'" + empDesignation($d) + "' needs " + "{$d.position}.{$expected}" +
                    " before it";
            }
        }
    }
    return "";
}

# entryProblem checks a whole entry; "" if it can be written.
func entryProblem(e as EmpEntry) {
    if (not regex.matches(ID_PATTERN, $e.haplotype.id)) {
        return "sample ID '" + $e.haplotype.id +
            "' may hold only letters, digits, '-', '_' and '/'";
    }
    if (len($e.haplogroup) == 0 or regex.matches(SPACE_PATTERN, $e.haplogroup)) {
        return "haplogroup '" + $e.haplogroup + "' must be one word; '?' if unknown";
    }
    if ($e.count < 0 or $e.count > MAX_COUNT) {
        return "count {$e.count} is outside 0-{MAX_COUNT}";
    }
    return haplotypeProblem($e.haplotype);
}

# canonical says whether differences are as mtHaplotype leaves them: strictly
# ascending by (position, insertion), each of a known kind, a deletion without
# a base. Two canonical lists name the same designations exactly when they are
# equal, so they compare without the call into lineage, whose argument copies
# were most of the time of counting.
func canonical(ds as list of lineage.Difference) {
    for (def i in 0..len($ds)) {
        def kind as string init $ds[$i].kind;
        if ($kind != "sub" and $kind != "ins" and $kind != "del") {
            return false;
        }
        if ($kind == "del" and len($ds[$i].base) > 0) {
            return false;
        }
        if ($i > 0) {
            def p as int init $ds[$i - 1].position;
            if ($ds[$i].position < $p or
                ($ds[$i].position == $p and $ds[$i].insertion <= $ds[$i - 1].insertion)) {
                return false;
            }
        }
    }
    return true;
}

# rangeKey is a set of ranges as one string: sorted, each once, as
# lineage.mtRangesAgree compares them.
func rangeKey(ranges as list of string) {
    def set as map of string to bool init {};
    for (def r in $ranges) {
        $set[$r] = true;
    }
    return strings.join(lists.sort(maps.keys($set)), " ");
}

# normalised rebuilds a haplotype through lineage.mtHaplotype, which sorts
# the differences and refuses a repeated or conflicting one; the "lineage"
# error comes back as this module's, after `where`.
func normalised(h as lineage.MtHaplotype, where as string) {
    def designations as list of string init [];
    for (def d in $h.differences) {
        $designations[] = empDesignation($d);
    }
    def out as lineage.MtHaplotype init $h;
    try {
        $out = lineage.mtHaplotype($h.id, $h.ranges, strings.join($designations, " "));
    } catch (e) {
        if ($e.kind == "lineage") {
            fail($where + $e.message);
        }
        throw $e;
    }
    return $out;
}

# --- reading ----------------------------------------------------------------

# commentText is a comment line without its "#" and one following space.
func commentText(raw as string) {
    def t as string init $raw[1..];
    if (strings.startsWith($t, " ")) {
        $t = $t[1..];
    }
    return $t;
}

# readEntry reads one haplotype line under the given ranges.
func readEntry(raw as string, n as int, ranges as list of string) {
    def at as string init "line {$n}: ";
    def cells as list of string init strings.split($raw, "\t");
    if (len($cells) < 3) {
        fail($at + "expected sample ID, haplogroup and count, separated by tabs");
    }
    def id as string init strings.trim($cells[0]);
    def countText as string init strings.trim($cells[2]);
    if (not regex.matches(COUNT_PATTERN, $countText)) {
        fail($at + "count '" + $countText + "' is not a whole number in 0-{MAX_COUNT}");
    }
    def designations as list of string init [];
    for (def cell in $cells[3..]) {
        def t as string init strings.trim($cell);
        if (len($t) == 0) {
            continue;
        }
        if (regex.matches(SPACE_PATTERN, $t)) {
            fail($at + "'" + $t + "' holds more than one variant; put each in its own cell");
        }
        $designations[] = $t;
    }
    def h as lineage.MtHaplotype init lineage.MtHaplotype{
        id: $id,
        ranges: $ranges,
        differences: []
    };
    try {
        $h = lineage.mtHaplotype($id, $ranges, strings.join($designations, " "));
    } catch (e) {
        if ($e.kind == "lineage") {
            fail($at + $e.message);
        }
        throw $e;
    }
    def entry as EmpEntry init EmpEntry{
        haplotype: $h,
        haplogroup: strings.trim($cells[1]),
        count: convert.toInt($countText),
        line: $n
    };
    def problem as string init entryProblem($entry);
    if (len($problem) > 0) {
        fail($at + $problem);
    }
    return $entry;
}

/**
 * Read an EMP file.
 *
 * Accepts LF or CRLF line ends, a leading byte-order mark, blank lines and
 * empty trailing cells (a haplotype equal to the rCRS has no difference).
 * Comment lines below the first `#!` line, among them a haplotype commented
 * out to exclude it, are skipped. Differences are checked as EMPOP's EMPCheck
 * checks them: the notation, a position within 1-16569 and within the
 * sequence range, no position named twice, and insertions numbered `.1`,
 * `.2`, ... without a gap. The reference base itself given as a difference
 * cannot be caught without the rCRS, which is not bundled.
 *
 * @param text {string} the file's content
 * @return {EmpFile} its header comments and haplotypes
 * @throws {Error} kind "emp" naming the line, for a haplotype before the first
 *   `#!` line, a malformed range, sample ID, haplogroup, count or difference,
 *   a repeated sample ID, or a difference EMPCheck would refuse; kind "limit"
 *   past `MAX_ENTRIES` haplotypes
 */
export func readEmp(text as string) {
    def lines as list of string init strings.split($text, "\n");
    def comments as list of string init [];
    def entries as list of EmpEntry init [];
    def ranges as list of string init [];
    def seen as map of string to int init {};
    for (def i in 0..len($lines)) {
        def raw as string init strings.trimRight($lines[$i]);
        def n as int init $i + 1;
        if ($i == 0 and strings.startsWith($raw, "﻿")) {
            $raw = $raw[1..];
        }
        if (len(strings.trim($raw)) == 0) {
            continue;
        }
        if (strings.startsWith($raw, "#!")) {
            def rest as string init strings.trim($raw[2..]);
            if (len($rest) == 0) {
                fail("line {$n}: the '#!' line names no sequence range");
            }
            $ranges = regex.split('\s+', $rest);
            for (def r in $ranges) {
                def problem as string init rangeProblem($r);
                if (len($problem) > 0) {
                    fail("line {$n}: " + $problem);
                }
            }
            continue;
        }
        if (strings.startsWith($raw, "#")) {
            if (len($ranges) == 0) {
                $comments[] = commentText($raw);
            }
            continue;
        }
        if (len($ranges) == 0) {
            fail("line {$n}: a haplotype before the first '#!' sequence range line");
        }
        if (len($entries) >= MAX_ENTRIES) {
            failLimit("line {$n}: more than {MAX_ENTRIES} haplotypes");
        }
        def entry as EmpEntry init readEntry($raw, $n, $ranges);
        def id as string init $entry.haplotype.id;
        if (maps.has($seen, $id)) {
            fail("line {$n}: sample ID '" + $id + "' is already on line {$seen[$id]}");
        }
        $seen[$id] = $n;
        $entries[] = $entry;
    }
    return EmpFile{comments: $comments, entries: $entries};
}

# --- building and writing ---------------------------------------------------

/**
 * An EMP entry from a haplotype, checked as `readEmp` checks a line.
 *
 * @param h {lineage.MtHaplotype} the haplotype; its `id` is the sample ID and
 *   its `ranges` the sequence ranges
 * @param haplogroup {string} the haplogroup, "?" if unknown
 * @param count {int} how often it was seen, 0 to `MAX_COUNT`
 * @return {EmpEntry} the entry, its differences sorted, `line` 0
 * @throws {Error} kind "emp" for anything `readEmp` would refuse
 */
export func empEntry(h as lineage.MtHaplotype, haplogroup as string, count as int) {
    def entry as EmpEntry init EmpEntry{
        haplotype: normalised($h, "sample '" + $h.id + "': "),
        haplogroup: $haplogroup,
        count: $count,
        line: 0
    };
    def problem as string init entryProblem($entry);
    if (len($problem) > 0) {
        fail("sample '" + $h.id + "': " + $problem);
    }
    return $entry;
}

# entryLine renders one haplotype line: the differences in range order, so a
# control region 16024-576 lists 16xxx before 73; an entry with no difference
# ends in an empty cell, as EMPOP's template writes it.
func entryLine(e as EmpEntry) {
    def cells as list of string init [$e.haplotype.id, $e.haplogroup, convert.toString($e.count)];
    def segs as list of Segment init segmentsOf($e.haplotype.ranges);
    for (def s in 0..len($segs)) {
        for (def d in $e.haplotype.differences) {
            if (segmentIndex($segs, $d.position) == $s) {
                $cells[] = empDesignation($d);
            }
        }
    }
    if (len($e.haplotype.differences) == 0) {
        $cells[] = "";
    }
    return strings.join($cells, "\t");
}

/**
 * Write an EMP file.
 *
 * The comments come first, each as `# ` and its text; a `#!` line precedes
 * the first entry and every entry whose ranges differ from the one before.
 * Differences are written in range order, deletions as `DEL`. Lines end in LF.
 * `readEmp` reads the result back to the same comments and entries (`line`
 * aside).
 *
 * @param f {EmpFile} the comments and entries
 * @return {string} the file's content
 * @throws {Error} kind "emp" naming the entry (1-based) for anything
 *   `readEmp` would refuse, a repeated sample ID, or a comment holding a line
 *   break
 */
export func writeEmp(f as EmpFile) {
    def out as list of string init [];
    for (def c in $f.comments) {
        if (strings.contains($c, "\n") or strings.contains($c, "\r")) {
            fail("comment '" + strings.trim($c) + "' holds a line break");
        }
        if (len($c) == 0) {
            $out[] = "#";
        } else {
            $out[] = "# " + $c;
        }
    }
    def last as list of string init [];
    def seen as map of string to int init {};
    for (def i in 0..len($f.entries)) {
        def at as string init "entry {$i + 1}: ";
        def e as EmpEntry init $f.entries[$i];
        if (not canonical($e.haplotype.differences)) {
            $e.haplotype = normalised($e.haplotype, $at);
        }
        def problem as string init entryProblem($e);
        if (len($problem) > 0) {
            fail($at + $problem);
        }
        if (maps.has($seen, $e.haplotype.id)) {
            fail($at + "sample ID '" + $e.haplotype.id + "' is already entry " +
                "{$seen[$e.haplotype.id]}");
        }
        $seen[$e.haplotype.id] = $i + 1;
        if ($i == 0 or $e.haplotype.ranges != $last) {
            $out[] = "#! " + strings.join($e.haplotype.ranges, " ");
            $last = $e.haplotype.ranges;
        }
        $out[] = entryLine($e);
    }
    return strings.join($out, "\n") + "\n";
}

# --- counting ---------------------------------------------------------------

/**
 * The database size `N` for `lineage`'s estimates: the sum of the counts.
 * An entry with count 0 is excluded and adds nothing.
 *
 * @param f {EmpFile} the file
 * @return {int} the summed counts
 */
export func empSampleSize(f as EmpFile) {
    def n as int init 0;
    for (def e in $f.entries) {
        $n = $n + $e.count;
    }
    return $n;
}

/**
 * How often a haplotype was seen, `k` for `lineage`'s estimates: the summed
 * counts of the entries with the same differences and the same sequence
 * ranges (`lineage.mtMatches` and `lineage.mtRangesAgree`).
 *
 * A haplotype typed over other ranges than the file's is never counted: trim
 * both to the shared range first.
 *
 * @param f {EmpFile} the file
 * @param h {lineage.MtHaplotype} the haplotype searched for
 * @return {int} the summed counts of its matches
 */
export func empCount(f as EmpFile, h as lineage.MtHaplotype) {
    def k as int init 0;
    def fast as bool init canonical($h.differences);
    def key as string init rangeKey($h.ranges);
    for (def e in $f.entries) {
        def same as bool init false;
        if ($fast and canonical($e.haplotype.differences)) {
            $same = $e.haplotype.differences == $h.differences;
        } else {
            $same = lineage.mtMatches($e.haplotype, $h);
        }
        if ($same and ($e.haplotype.ranges == $h.ranges or rangeKey($e.haplotype.ranges) == $key)) {
            $k = $k + $e.count;
        }
    }
    return $k;
}
