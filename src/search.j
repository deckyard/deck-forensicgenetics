# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Database searching: one query profile against a list of profiles.
 *
 * A direct search compares genotypes locus by locus. At every locus both
 * profiles are typed at it records the alleles shared identical by state
 * (IBS 0, 1 or 2: the size of the multiset intersection of the two genotypes)
 * and whether the genotypes are identical, consistent with drop-out, or a
 * mismatch. Each candidate is then classed:
 *
 * - `"exact"`: at least `minLoci` loci compared, every one identical;
 * - `"partial"`: at least `minLoci` compared, not exact, and at most
 *   `maxMismatches` mismatching loci;
 * - `"mismatch"`: at least `minLoci` compared, more mismatches than allowed;
 * - `"insufficient"`: fewer than `minLoci` loci compared.
 *
 * With `allowDropout`, a locus where the query is homozygous `a/a` and the
 * candidate is a heterozygote carrying `a` is drop-out-consistent rather than
 * a mismatch: the query (the evidence) may have lost the candidate's other
 * allele. It is directional - an allele the query shows and the candidate
 * lacks is always a mismatch. Loci only one profile carries are not compared,
 * which is what lets a partial profile be searched at all.
 *
 * A familial search ranks candidates by the kinship likelihood ratio of a
 * chosen relationship against unrelated, `kinship.kinshipLR`, so the theta = 0
 * convention and the "never drop a covered locus" rule of `kinship.j` apply
 * unchanged.
 *
 * Locus names are compared as written: normalise the query and the database
 * with `panels.normalizeProfile` first when they came from different sources.
 * @module search
 * @example
 * import "@deckyard/forensicgenetics/search.j" as search;
 * import "@deckyard/forensicgenetics/kinship.j" as kinship;
 *
 * def opts as search.SearchOptions init search.searchOptions(10, 1, true);
 * def ranked as list of search.Comparison init search.search($stain, $database, $opts);
 * def leads as list of search.Comparison init search.hits($ranked);
 *
 * def kin as list of search.FamilialHit init search.familialSearch(
 *     $stain, $database, kinship.FULL_SIBLINGS, $db, 15);
 */
use maps;
use lists;

import "./forensicgenetics.j" as forensics;
import "./kinship.j" as kinship;

/**
 * Every compared locus identical.
 */
export def const VERDICT_EXACT as string init "exact";

/**
 * Not exact, but no more mismatching loci than allowed.
 */
export def const VERDICT_PARTIAL as string init "partial";

/**
 * More mismatching loci than allowed.
 */
export def const VERDICT_MISMATCH as string init "mismatch";

/**
 * Fewer loci compared than `minLoci`; no conclusion.
 */
export def const VERDICT_INSUFFICIENT as string init "insufficient";

/**
 * How a direct search classes a candidate.
 * @field minLoci {int} the fewest loci that must be compared for any verdict
 *   other than "insufficient"; at least 1
 * @field maxMismatches {int} the most mismatching loci a "partial" candidate
 *   may have; at least 0
 * @field allowDropout {bool} whether a homozygous query locus whose allele the
 *   candidate carries counts as drop-out rather than a mismatch
 */
export def struct SearchOptions {
    minLoci as int,
    maxMismatches as int,
    allowDropout as bool
};

/**
 * One candidate of a direct search, compared with the query.
 * @field id {string} the candidate profile's id
 * @field index {int} the candidate's position in the database list
 * @field lociCompared {int} loci both profiles are typed at
 * @field matchingLoci {int} compared loci with identical genotypes
 * @field dropoutLoci {int} compared loci counted as drop-out (`allowDropout`)
 * @field mismatchedLoci {list of string} the mismatching loci, in query order
 * @field ibs {map of string to int} compared locus -> alleles shared (0, 1, 2)
 * @field ibs0 {int} compared loci sharing no allele
 * @field ibs1 {int} compared loci sharing one allele
 * @field ibs2 {int} compared loci sharing both
 * @field sharedAlleles {int} the sum of `ibs` over the compared loci
 * @field verdict {string} "exact", "partial", "mismatch" or "insufficient"
 */
export def struct Comparison {
    id as string,
    index as int,
    lociCompared as int,
    matchingLoci as int,
    dropoutLoci as int,
    mismatchedLoci as list of string,
    ibs as map of string to int,
    ibs0 as int,
    ibs1 as int,
    ibs2 as int,
    sharedAlleles as int,
    verdict as string
};

/**
 * One candidate of a familial search.
 * @field id {string} the candidate profile's id
 * @field index {int} the candidate's position in the database list
 * @field lociCompared {int} loci both profiles are typed at
 * @field scored {bool} false when fewer than `minLoci` loci were compared;
 *   `lr` is then 0.0 and means nothing
 * @field lr {float} `kinship.kinshipLR` of the relationship against unrelated
 */
export def struct FamilialHit {
    id as string,
    index as int,
    lociCompared as int,
    scored as bool,
    lr as float
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "search", message: $msg, file: "", line: 0, col: 0};
}

# --- options ----------------------------------------------------------------

/**
 * Search options, validated.
 *
 * There is no default: the right `minLoci` and `maxMismatches` depend on the
 * kit, the database's policy and the quality of the query.
 * @param minLoci {int} the fewest loci compared for a verdict; at least 1
 * @param maxMismatches {int} mismatching loci a "partial" may have; at least 0
 * @param allowDropout {bool} whether a homozygous query locus whose allele
 *   the candidate carries counts as drop-out
 * @return {SearchOptions} the options
 * @throws {Error} kind "search" if `minLoci` is below 1 or `maxMismatches`
 *   is negative
 */
export func searchOptions(minLoci as int, maxMismatches as int, allowDropout as bool) {
    if ($minLoci < 1) {
        fail("minLoci must be at least 1");
    }
    if ($maxMismatches < 0) {
        fail("maxMismatches must not be negative");
    }
    return SearchOptions{
        minLoci: $minLoci,
        maxMismatches: $maxMismatches,
        allowDropout: $allowDropout
    };
}

# checkOptions repeats searchOptions' checks for an options literal that
# bypassed it.
func checkOptions(opts as SearchOptions) {
    searchOptions($opts.minLoci, $opts.maxMismatches, $opts.allowDropout);
}

# --- one locus --------------------------------------------------------------

/**
 * The number of alleles two genotypes share identical by state: the size of
 * the multiset intersection of the two allele pairs.
 *
 * `10/11` and `10/11` share 2; `10/10` and `10/11` share 1 (the second `10` of
 * the homozygote has no partner); `10/10` and `10/10` share 2; `10/11` and
 * `12/13` share 0.
 * @param ga {Genotype} one genotype
 * @param gb {Genotype} the other, at the same locus
 * @return {int} 0, 1 or 2
 * @throws {Error} kind "search" if the genotypes are at different loci
 */
export func allelesShared(ga as forensics.Genotype, gb as forensics.Genotype) {
    if ($ga.locus != $gb.locus) {
        fail("cannot compare a genotype at '" + $ga.locus + "' with one at '" + $gb.locus + "'");
    }
    return ibsCount($ga, $gb);
}

# ibsCount is allelesShared without the locus check, for the inner loop where
# both genotypes come from the same key.
func ibsCount(ga as forensics.Genotype, gb as forensics.Genotype) {
    # Match ga.a first, then ga.b against whatever of gb is left.
    if ($ga.a == $gb.a) {
        if ($ga.b == $gb.b) {
            return 2;
        }
        return 1;
    }
    if ($ga.a == $gb.b) {
        if ($ga.b == $gb.a) {
            return 2;
        }
        return 1;
    }
    if ($ga.b == $gb.a or $ga.b == $gb.b) {
        return 1;
    }
    return 0;
}

# sameGenotype compares two genotypes as unordered allele pairs, so a Genotype
# literal written out of canonical order still matches.
func sameGenotype(x as forensics.Genotype, y as forensics.Genotype) {
    return ($x.a == $y.a and $x.b == $y.b) or ($x.a == $y.b and $x.b == $y.a);
}

# dropoutConsistent: the query is homozygous and the candidate, not identical,
# carries the query's allele - the candidate's other allele may have dropped
# out of the query.
func dropoutConsistent(q as forensics.Genotype, c as forensics.Genotype) {
    return $q.a == $q.b and $c.a != $c.b and ($c.a == $q.a or $c.b == $q.a);
}

# --- direct search ----------------------------------------------------------

/**
 * One candidate compared with the query: loci compared, identical, drop-out
 * and mismatching loci, allele sharing, and the verdict. `index` is 0; `search`
 * sets it to the candidate's position.
 * @param query {Profile} the profile searched for (the evidence)
 * @param candidate {Profile} one database profile
 * @param opts {SearchOptions} how to class the candidate
 * @return {Comparison} the comparison
 * @throws {Error} kind "search" on invalid options or an empty query
 */
export func compareProfiles(
    query as forensics.Profile,
    candidate as forensics.Profile,
    opts as SearchOptions) {
    checkOptions($opts);
    if (len($query.loci) == 0) {
        fail("query profile '" + $query.id + "' has no genotypes");
    }
    return comparison($query, $candidate, $opts, 0);
}

# comparison does the work of compareProfiles once the options and query are
# known good, so search validates them once rather than per candidate. Both
# profiles are read-only parameters, so neither is copied.
func comparison(
    query as forensics.Profile,
    candidate as forensics.Profile,
    opts as SearchOptions,
    index as int) {
    def out as Comparison init Comparison{
        id: $candidate.id,
        index: $index,
        lociCompared: 0,
        matchingLoci: 0,
        dropoutLoci: 0,
        mismatchedLoci: [],
        ibs: {},
        ibs0: 0,
        ibs1: 0,
        ibs2: 0,
        sharedAlleles: 0,
        verdict: VERDICT_INSUFFICIENT
    };
    for (def locus in $query.loci) {
        if (not maps.has($candidate.loci, $locus)) {
            continue;
        }
        $out.lociCompared = $out.lociCompared + 1;
        def n as int init ibsCount($query.loci[$locus], $candidate.loci[$locus]);
        $out.ibs[$locus] = $n;
        $out.sharedAlleles = $out.sharedAlleles + $n;
        if ($n == 0) {
            $out.ibs0 = $out.ibs0 + 1;
        } elseif ($n == 1) {
            $out.ibs1 = $out.ibs1 + 1;
        } else {
            $out.ibs2 = $out.ibs2 + 1;
        }
        if (sameGenotype($query.loci[$locus], $candidate.loci[$locus])) {
            $out.matchingLoci = $out.matchingLoci + 1;
        } elseif ($opts.allowDropout and
            dropoutConsistent($query.loci[$locus], $candidate.loci[$locus])) {
            $out.dropoutLoci = $out.dropoutLoci + 1;
        } else {
            $out.mismatchedLoci[] = $locus;
        }
    }
    $out.verdict = verdictOf($out.lociCompared, $out.matchingLoci, len($out.mismatchedLoci), $opts);
    return $out;
}

# verdictOf classes a candidate from its counts.
func verdictOf(compared as int, matching as int, mismatches as int, opts as SearchOptions) {
    if ($compared < $opts.minLoci) {
        return VERDICT_INSUFFICIENT;
    }
    if ($matching == $compared) {
        return VERDICT_EXACT;
    }
    if ($mismatches <= $opts.maxMismatches) {
        return VERDICT_PARTIAL;
    }
    return VERDICT_MISMATCH;
}

# verdictRank orders verdicts best first.
func verdictRank(verdict as string) {
    if ($verdict == VERDICT_EXACT) {
        return 0;
    }
    if ($verdict == VERDICT_PARTIAL) {
        return 1;
    }
    if ($verdict == VERDICT_MISMATCH) {
        return 2;
    }
    return 3;
}

# A RankKey carries only what ranking reads, so the sort passes move small
# structs rather than whole Comparisons with their per-locus maps.
def struct RankKey {
    pos as int,
    verdict as int,
    mismatches as int,
    matching as int,
    shared as int
};

func bySharedDesc(k as RankKey) {
    return 0 - $k.shared;
}

func byMatchingDesc(k as RankKey) {
    return 0 - $k.matching;
}

func byMismatches(k as RankKey) {
    return $k.mismatches;
}

func byVerdict(k as RankKey) {
    return $k.verdict;
}

/**
 * The query compared with every profile of a database, ranked.
 *
 * Order: verdict (exact, partial, mismatch, insufficient), then fewer
 * mismatching loci, then more identical loci, then more shared alleles, then
 * database order. Every candidate is returned, "insufficient" ones last;
 * `hits` keeps only the exact and partial ones.
 * @param query {Profile} the profile searched for (the evidence)
 * @param database {list of Profile} the profiles searched
 * @param opts {SearchOptions} how to class each candidate
 * @return {list of Comparison} one per database profile, best first
 * @throws {Error} kind "search" on invalid options or an empty query
 */
export func search(
    query as forensics.Profile,
    database as list of forensics.Profile,
    opts as SearchOptions) {
    checkOptions($opts);
    if (len($query.loci) == 0) {
        fail("query profile '" + $query.id + "' has no genotypes");
    }
    def found as list of Comparison init [];
    def keys as list of RankKey init [];
    # Indexed, not `for (def c in $database)`: that binds a copy of each
    # profile, and its whole locus map, on every pass.
    for (def i in 0..len($database)) {
        def c as Comparison init comparison($query, $database[$i], $opts, $i);
        $keys[] = RankKey{
            pos: $i,
            verdict: verdictRank($c.verdict),
            mismatches: len($c.mismatchedLoci),
            matching: $c.matchingLoci,
            shared: $c.sharedAlleles
        };
        $found[] = $c;
    }
    # lists.sortBy is stable, so sorting by the least significant key first and
    # the most significant last gives the lexicographic order, ties left in
    # database order.
    $keys = lists.sortBy($keys, bySharedDesc);
    $keys = lists.sortBy($keys, byMatchingDesc);
    $keys = lists.sortBy($keys, byMismatches);
    $keys = lists.sortBy($keys, byVerdict);
    def out as list of Comparison init [];
    for (def k in $keys) {
        $out[] = $found[$k.pos];
    }
    return $out;
}

/**
 * The exact and partial candidates of a search result, order kept.
 * @param results {list of Comparison} a `search` result
 * @return {list of Comparison} the candidates with verdict "exact" or "partial"
 */
export func hits(results as list of Comparison) {
    def out as list of Comparison init [];
    for (def i in 0..len($results)) {
        if ($results[$i].verdict == VERDICT_EXACT or $results[$i].verdict == VERDICT_PARTIAL) {
            $out[] = $results[$i];
        }
    }
    return $out;
}

# --- familial search --------------------------------------------------------

# sharedLoci lists the loci both profiles are typed at, in the query's order -
# the loci, and the order, kinship.kinshipLR multiplies over.
func sharedLoci(query as forensics.Profile, candidate as forensics.Profile) {
    def out as list of string init [];
    for (def locus in $query.loci) {
        if (maps.has($candidate.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

# allCovered reports whether the frequency table covers every listed locus.
func allCovered(fdb as forensics.FrequencyDb, loci as list of string) {
    for (def locus in $loci) {
        if (not maps.has($fdb.loci, $locus)) {
            return false;
        }
    }
    return true;
}

func byLrDesc(h as FamilialHit) {
    return 0.0 - $h.lr;
}

func byScored(h as FamilialHit) {
    if ($h.scored) {
        return 0;
    }
    return 1;
}

/**
 * Database candidates ranked by the kinship likelihood ratio for `rel`
 * against unrelated.
 *
 * Each candidate sharing at least `minLoci` typed loci with the query is
 * scored with the value of `kinship.kinshipLR(query, candidate, rel, fdb)`,
 * computed from `kinship.kinshipLRAt` memoised per locus and candidate
 * genotype; the rest are returned unscored, after the scored ones. Ties keep database order. As in
 * `kinship.kinshipLR`, a locus both profiles carry but `fdb` does not cover
 * throws rather than being dropped: narrow the profiles to the table's loci
 * first (`panels.restrictProfile`, `forensicgenetics.restrict`).
 *
 * A familial search produces leads, not identifications: in a large database
 * unrelated people reach high LRs by chance, and half siblings, grandparents
 * and avuncular relatives are indistinguishable.
 * @param query {Profile} the profile searched for
 * @param database {list of Profile} the profiles searched
 * @param rel {Relationship} the relationship to test, e.g. `kinship.FULL_SIBLINGS`
 * @param fdb {FrequencyDb} the reference frequencies
 * @param minLoci {int} the fewest shared loci for a candidate to be scored;
 *   at least 1
 * @return {list of FamilialHit} one per database profile, highest LR first
 * @throws {Error} kind "search" on an empty query or `minLoci` below 1;
 *   kind "kinship" (from `kinship.kinshipLR`) on a shared locus `fdb` does
 *   not cover
 */
export func familialSearch(
    query as forensics.Profile,
    database as list of forensics.Profile,
    rel as kinship.Relationship,
    fdb as forensics.FrequencyDb,
    minLoci as int) {
    if ($minLoci < 1) {
        fail("minLoci must be at least 1");
    }
    if (len($query.loci) == 0) {
        fail("query profile '" + $query.id + "' has no genotypes");
    }
    # The query is fixed, so a per-locus LR depends only on the candidate's
    # genotype there; a locus has a few dozen distinct genotypes against
    # thousands of candidates. kinshipLRAt is memoised per (locus, genotype)
    # and multiplied in kinshipLR's own order from 1.0, so the product is
    # bit-identical to kinship.kinshipLR (asserted in the overlay).
    def cache as map of string to float init {};
    def out as list of FamilialHit init [];
    for (def i in 0..len($database)) {
        def loci as list of string init sharedLoci($query, $database[$i]);
        def h as FamilialHit init FamilialHit{
            id: $database[$i].id,
            index: $i,
            lociCompared: len($loci),
            scored: false,
            lr: 0.0
        };
        if (len($loci) < $minLoci) {
            $out[] = $h;
            continue;
        }
        if (not allCovered($fdb, $loci)) {
            # Let kinshipLR raise its own "kinship" error for the uncovered
            # locus.
            kinship.kinshipLR($query, $database[$i], $rel, $fdb);
        }
        def lr as float init 1.0;
        for (def locus in $loci) {
            def g as forensics.Genotype init $database[$i].loci[$locus];
            def key as string init $locus + "\t" + $g.a + "\t" + $g.b;
            if (not maps.has($cache, $key)) {
                $cache[$key] = kinship.kinshipLRAt(
                    $fdb,
                    $query.loci[$locus],
                    $g,
                    $rel,
                    kinship.UNRELATED);
            }
            $lr = $lr * $cache[$key];
        }
        $h.scored = true;
        $h.lr = $lr;
        $out[] = $h;
    }
    $out = lists.sortBy($out, byLrDesc);
    return lists.sortBy($out, byScored);
}
