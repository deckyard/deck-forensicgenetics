# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Autosomal STR profiles, reference allele frequencies and the single-source
 * statistics built from them: NRC II genotype frequencies with theta, the
 * random match probability (RMP) and its likelihood ratio, and CPI / CPE for a
 * mixture.
 *
 * Alleles are strings: `"9.3"` and `"15.2"` are ordinary designations, and
 * `"15"` and `"15.0"` are different keys. `alleleValue` converts one to a
 * repeat count where arithmetic is meant (the mutation model in `kinship.j`).
 *
 * Every allele frequency is floored at the NRC II 5/(2N), using the locus's
 * own N; that is why `FrequencyDb` stores a sample size per locus.
 *
 * The module reads no published format. Build a `FrequencyDb` with `db` and
 * `withLocus`, or convert one with a reader deck (`@deckyard/strider` for
 * STRidER XML). For a donor whose population is uncertain use
 * `randomMatchProbabilityAcross`, not a table pooled with `combineFrequencies`,
 * which biases the RMP downwards.
 *
 * Sibling modules: `kinship.j` (pairwise LRs, paternity), `pedigree.j`
 * (whole-pedigree likelihoods), `lineage.j` (mtDNA and Y-STR).
 * @module forensicgenetics
 * @see https://strider.online/formulae - the STRidER formulae page these
 *   estimators follow
 * @see https://github.com/deckyard/deck-strider-frequency-converter - the reader
 *   deck that converts a STRidER export into a `FrequencyDb`
 * @example
 * import "@deckyard/forensicgenetics/" as forensics;
 *
 * def db as forensics.FrequencyDb init forensics.withLocus(
 *     forensics.db("example"), "D3S1358", 1034, {"15": 0.4, "16": 0.6});
 * def suspect as forensics.Profile init forensics.profile("S1", [
 *     forensics.genotype("D3S1358", "15", "16"),
 *     forensics.genotype("TH01", "9.3", "9.3")
 * ]);
 * def t as float init forensics.THETA_GENERAL;
 * def rmp as float init forensics.randomMatchProbability($db, $suspect, $t);
 * def lr as float init forensics.matchLikelihoodRatio($db, $suspect, $t);
 */
use math;
use maps;
use lists;
use strings;
use convert;

/**
 * The NRC II theta for a broad, well-mixed population.
 *
 * Theta (Fst) is the co-ancestry coefficient: the probability that two alleles
 * drawn from the same sub-population are identical by descent. A larger theta
 * gives a larger, more conservative match probability.
 */
export def const THETA_GENERAL as float init 0.01;

/**
 * The NRC II theta recommended for small, isolated or endogamous populations.
 */
export def const THETA_ISOLATED as float init 0.03;

/**
 * Theta zero: the plain product rule with no sub-population correction.
 *
 * Correct only for an idealised panmictic population; not for casework.
 */
export def const THETA_NONE as float init 0.0;

/**
 * The numerator of the NRC II minimum-allele-frequency floor, 5/(2N).
 */
export def const MIN_ALLELE_COUNT as float init 5.0;

/**
 * The largest allele-frequency sum `withLocus` accepts for one locus.
 *
 * Published tables round, so a locus rarely sums to exactly 1. A sum above
 * this means the table was misread.
 */
export def const MAX_FREQUENCY_SUM as float init 1.05;

/**
 * One locus of a reference frequency table.
 * @field locus {string} the locus name, e.g. "D3S1358"
 * @field n {int} individuals typed at this locus (2N alleles); 0 if unknown
 * @field alleles {map of string to float} allele label -> published frequency
 */
export def struct LocusFrequencies {
    locus as string,
    n as int,
    alleles as map of string to float
};

/**
 * A reference allele-frequency database, one entry per locus.
 * @field name {string} the population / dataset label
 * @field loci {map of string to LocusFrequencies} locus name -> its frequencies
 */
export def struct FrequencyDb {
    name as string,
    loci as map of string to LocusFrequencies
};

/**
 * An unordered pair of alleles at one autosomal locus.
 *
 * Build one with `genotype`, which orders the pair so that both writing
 * orders give the same value.
 * @field locus {string} the locus name
 * @field a {string} the lower allele in canonical order
 * @field b {string} the higher allele; equal to `a` for a homozygote
 */
export def struct Genotype {
    locus as string,
    a as string,
    b as string
};

/**
 * A single-source STR profile: one genotype per locus.
 * @field id {string} a caller-chosen label for the profile
 * @field loci {map of string to Genotype} locus name -> the genotype seen there
 */
export def struct Profile {
    id as string,
    loci as map of string to Genotype
};

/**
 * A mixed stain: the alleles seen at each locus, with no assumption about the
 * number of contributors or how the alleles pair up.
 * @field id {string} a caller-chosen label
 * @field loci {map of string to list of string} locus name -> the alleles observed
 */
export def struct Mixture {
    id as string,
    loci as map of string to list of string
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "forensicgenetics", message: $msg, file: "", line: 0, col: 0};
}

# --- allele helpers ---------------------------------------------------------

/**
 * The repeat count an allele label denotes.
 *
 * `"15"` gives 15.0 and `"9.3"` gives 9.3. `"X"`, `"Y"` and `"OL"` have no
 * repeat count; test with `isNumericAllele` first.
 * @param allele {string} the allele label
 * @return {float} the repeat count
 * @throws {Error} kind "forensicgenetics" if the label is not numeric
 */
export func alleleValue(allele as string) {
    if (not isNumericAllele($allele)) {
        fail("allele '" + $allele + "' has no numeric repeat count");
    }
    return convert.toFloat($allele);
}

/**
 * Whether an allele label is a repeat count: digits with at most one `.`.
 * @param allele {string} the allele label
 * @return {bool} true if `alleleValue` would succeed
 */
export func isNumericAllele(allele as string) {
    if (len($allele) == 0) {
        return false;
    }
    def dots as int init 0;
    def digits as int init 0;
    for (def ch in strings.chars($allele)) {
        if ($ch == ".") {
            $dots = $dots + 1;
        } elseif (strings.indexOf("0123456789", $ch) >= 0) {
            $digits = $digits + 1;
        } else {
            return false;
        }
    }
    return $digits > 0 and $dots <= 1;
}

# compareAlleles orders two allele labels: numerically when both carry a repeat
# count, lexicographically otherwise, so a panel mixing "15" with "X" still has
# a total order. Returns -1, 0 or 1.
func compareAlleles(x as string, y as string) {
    if (isNumericAllele($x) and isNumericAllele($y)) {
        def vx as float init convert.toFloat($x);
        def vy as float init convert.toFloat($y);
        if ($vx < $vy) {
            return -1;
        }
        if ($vx > $vy) {
            return 1;
        }
        return 0;
    }
    if ($x < $y) {
        return -1;
    }
    if ($x > $y) {
        return 1;
    }
    return 0;
}

# distinct returns the input with repeats removed, order preserved. Membership
# is a hash lookup rather than a scan of what has been kept so far, so a locus
# listing many alleles stays linear.
func distinct(xs as list of string) {
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

# --- frequency database -----------------------------------------------------

/**
 * An empty frequency database.
 * @param name {string} the population or dataset label
 * @return {FrequencyDb} a database with no loci
 */
export func db(name as string) {
    return FrequencyDb{name: $name, loci: {}};
}

# checkLocus validates one locus table, throwing on anything that would make it
# not a probability distribution. Factored out of withLocus so a function that
# already owns its FrequencyDb can validate an entry and write it directly,
# rather than paying withLocus's whole-database copy per locus.
func checkLocus(locus as string, n as int, freqs as map of string to float) {
    if (len($locus) == 0) {
        fail("locus name must not be empty");
    }
    if ($n < 0) {
        fail("locus '" + $locus + "': sample size must not be negative");
    }
    def total as float init 0.0;
    for (def a in $freqs) {
        if ($freqs[$a] <= 0.0 or $freqs[$a] > 1.0) {
            fail("locus '" + $locus + "' allele '" + $a + "': frequency must be in (0, 1]");
        }
        $total = $total + $freqs[$a];
    }
    # A locus may legitimately list only part of its alleles, so a sum below 1
    # is fine. A sum meaningfully above 1 is not a rounding artefact - it is a
    # misread table, and left alone it yields genotype "probabilities" above 1.
    if ($total > MAX_FREQUENCY_SUM) {
        fail("locus '" + $locus + "': allele frequencies sum to {$total}, which is not a" +
            " probability distribution - the table was very likely misread");
    }
}

/**
 * A copy of `fdb` with one locus added or replaced.
 *
 * Frequencies are stored as given; `alleleFrequency` applies the 5/(2N) floor
 * on read. A locus may list only some of its alleles, so a sum below 1 is
 * accepted.
 * @param fdb {FrequencyDb} the database to extend
 * @param locus {string} the locus name
 * @param n {int} individuals typed at this locus; 0 if unknown (no floor)
 * @param freqs {map of string to float} allele label -> published frequency
 * @return {FrequencyDb} the extended copy
 * @throws {Error} kind "forensicgenetics" on an empty locus name, a negative `n`, a
 *   frequency outside (0, 1], or frequencies summing past `MAX_FREQUENCY_SUM`
 */
export func withLocus(
    fdb as FrequencyDb,
    locus as string,
    n as int,
    freqs as map of string to float) {
    checkLocus($locus, $n, $freqs);
    def out as FrequencyDb init $fdb;
    $out.loci[$locus] = LocusFrequencies{locus: $locus, n: $n, alleles: $freqs};
    return $out;
}

# addObservedCounts accumulates one locus table's observed allele counts
# (frequency * 2N) into a running total, so pooling several tables is a sum of
# counts rather than an average of frequencies.
func addObservedCounts(counts as map of string to float, lf as LocusFrequencies) {
    def twoN as float init 2.0 * convert.toFloat($lf.n);
    def out as map of string to float init $counts;
    for (def a in $lf.alleles) {
        def seen as float init $lf.alleles[$a] * $twoN;
        if (maps.has($out, $a)) {
            $out[$a] = $out[$a] + $seen;
        } else {
            $out[$a] = $seen;
        }
    }
    return $out;
}

/**
 * Several frequency tables pooled into one, weighted by sample size.
 *
 * Each allele's pooled frequency is `sum(p_i * 2N_i) / sum(2N_i)`; the pooled
 * N is the sum of the inputs' N. Only loci every input carries are pooled;
 * call `loci` on the result to see which.
 *
 * A pooled table describes one merged population. Merging populations whose
 * frequencies differ leaves it short of heterozygotes (the Wahlund effect), so
 * `randomMatchProbability` on it is biased downwards and overstates the
 * evidence. For "the donor was German or Austrian" use
 * `randomMatchProbabilityAcross`.
 * @param dbs {list of FrequencyDb} the tables to pool; at least one
 * @param name {string} the dataset label for the pooled table
 * @return {FrequencyDb} the pooled table, over the shared loci
 * @throws {Error} kind "forensicgenetics" on an empty list, when the inputs
 *   share no locus, or when a shared locus has no sample size anywhere
 */
export func combineFrequencies(dbs as list of FrequencyDb, name as string) {
    if (len($dbs) == 0) {
        fail("combineFrequencies needs at least one database");
    }
    # $dbs is walked positionally throughout. `for (def d in $dbs)` binds a
    # fresh copy of the FrequencyDb on every pass, and copying a database
    # copies its whole locus map - so doing that once per locus, as the inner
    # loop below would, is quadratic in the number of loci. Indexing reads the
    # element in place instead.
    def shared as list of string init loci($dbs[0]);
    for (def i in 0..len($dbs)) {
        def table as map of string to LocusFrequencies init $dbs[$i].loci;
        def keep as list of string init [];
        for (def locus in $shared) {
            if (maps.has($table, $locus)) {
                $keep[] = $locus;
            }
        }
        $shared = $keep;
    }
    if (len($shared) == 0) {
        fail("the databases share no locus, so there is nothing to pool");
    }
    def out as FrequencyDb init db($name);
    for (def locus in $shared) {
        def counts as map of string to float init {};
        def alleles as float init 0.0;
        def pooledN as int init 0;
        for (def i in 0..len($dbs)) {
            def lf as LocusFrequencies init $dbs[$i].loci[$locus];
            $pooledN = $pooledN + $lf.n;
            $alleles = $alleles + 2.0 * convert.toFloat($lf.n);
            $counts = addObservedCounts($counts, $lf);
        }
        if ($alleles <= 0.0) {
            fail("locus '" + $locus + "': no input states a sample size, so it" +
                " cannot be weighted");
        }
        def pooled as map of string to float init {};
        for (def a in $counts) {
            $pooled[$a] = $counts[$a] / $alleles;
        }
        # $out is this function's own value, so writing into it directly costs
        # nothing; routing each locus through withLocus would copy every locus
        # pooled so far, which is the other half of the same quadratic.
        checkLocus($locus, $pooledN, $pooled);
        $out.loci[$locus] = LocusFrequencies{locus: $locus, n: $pooledN, alleles: $pooled};
    }
    return $out;
}

/**
 * The loci a database covers, in insertion order.
 * @param fdb {FrequencyDb} the database
 * @return {list of string} the locus names
 */
export func loci(fdb as FrequencyDb) {
    return maps.keys($fdb.loci);
}

/**
 * Whether a database covers a locus.
 * @param fdb {FrequencyDb} the database
 * @param locus {string} the locus name
 * @return {bool} true if the locus is present
 */
export func hasLocus(fdb as FrequencyDb, locus as string) {
    return maps.has($fdb.loci, $locus);
}

/**
 * The number of individuals typed at a locus.
 * @param fdb {FrequencyDb} the database
 * @param locus {string} the locus name
 * @return {int} the sample size N; 0 if unknown
 * @throws {Error} kind "forensicgenetics" if the database has no such locus
 */
export func sampleSize(fdb as FrequencyDb, locus as string) {
    return lookup($fdb, $locus).n;
}

# lookup fetches a locus with a message that names the missing locus, which is
# far more useful than the runtime's generic missing-key error.
func lookup(fdb as FrequencyDb, locus as string) {
    if (not maps.has($fdb.loci, $locus)) {
        fail("frequency database '" + $fdb.name + "' has no locus '" + $locus + "'");
    }
    return $fdb.loci[$locus];
}

/**
 * The NRC II minimum allele frequency for a locus: 5/(2N).
 *
 * No allele is scored rarer than this, so an unobserved allele cannot drive a
 * match probability to zero. A locus with `n == 0` has no floor.
 * @param fdb {FrequencyDb} the database
 * @param locus {string} the locus name
 * @return {float} the floor, or 0.0 when N is unknown
 * @throws {Error} kind "forensicgenetics" if the database has no such locus
 */
export func minFrequency(fdb as FrequencyDb, locus as string) {
    return alleleFloor(lookup($fdb, $locus).n);
}

# alleleFloor is minFrequency once the locus is already in hand. It takes the
# sample size rather than the LocusFrequencies so that a caller holding one
# does not look the locus up a second time - lookup returns by value, so every
# call deep-copies the locus's whole allele map.
func alleleFloor(n as int) {
    if ($n <= 0) {
        return 0.0;
    }
    return MIN_ALLELE_COUNT / (2.0 * convert.toFloat($n));
}

/**
 * The frequency of one allele, floored at 5/(2N) and capped at 1.0.
 *
 * An allele the table does not list is scored at the floor.
 * @param fdb {FrequencyDb} the database
 * @param locus {string} the locus name
 * @param allele {string} the allele label
 * @return {float} the floored frequency, in (0, 1]
 * @throws {Error} kind "forensicgenetics" if the database has no such locus, or
 *   if the allele is unlisted at a locus with `n == 0`
 */
export func alleleFrequency(fdb as FrequencyDb, locus as string, allele as string) {
    def lf as LocusFrequencies init lookup($fdb, $locus);
    def p as float init 0.0;
    if (maps.has($lf.alleles, $allele)) {
        $p = $lf.alleles[$allele];
    }
    def floor as float init alleleFloor($lf.n);
    if ($p < $floor) {
        $p = $floor;
    }
    if ($p > 1.0) {
        $p = 1.0;
    }
    if ($p <= 0.0) {
        fail("locus '" + $locus + "' allele '" + $allele +
            "': no frequency and no 5/2N floor (the table is missing its sample size)");
    }
    return $p;
}

# --- genotypes and profiles -------------------------------------------------

/**
 * A genotype at one locus, its allele pair in canonical order.
 *
 * `genotype("TH01", "9.3", "6")` and `genotype("TH01", "6", "9.3")` are equal.
 * Two numeric alleles order by repeat count; otherwise the order is
 * lexicographic.
 * @param locus {string} the locus name
 * @param a {string} one allele label
 * @param b {string} the other allele label
 * @return {Genotype} the canonical genotype
 * @throws {Error} kind "forensicgenetics" on an empty locus name or allele label
 */
export func genotype(locus as string, a as string, b as string) {
    if (len($locus) == 0) {
        fail("locus name must not be empty");
    }
    if (len($a) == 0 or len($b) == 0) {
        fail("locus '" + $locus + "': allele labels must not be empty");
    }
    if (compareAlleles($a, $b) > 0) {
        return Genotype{locus: $locus, a: $b, b: $a};
    }
    return Genotype{locus: $locus, a: $a, b: $b};
}

/**
 * Whether a genotype carries the same allele twice.
 * @param g {Genotype} the genotype to test
 * @return {bool} true for a homozygote
 */
export func isHomozygous(g as Genotype) {
    return $g.a == $g.b;
}

/**
 * The distinct alleles of a genotype: one label for a homozygote, two for a
 * heterozygote.
 * @param g {Genotype} the genotype to read
 * @return {list of string} the distinct allele labels, in canonical order
 */
export func alleles(g as Genotype) {
    if (isHomozygous($g)) {
        return [$g.a];
    }
    return [$g.a, $g.b];
}

/**
 * A profile from a list of genotypes.
 *
 * Builds the profile in one pass; prefer it to folding `withGenotype` over a
 * list, which copies the profile once per genotype.
 * @param id {string} a label for the profile
 * @param genotypes {list of Genotype} one genotype per locus
 * @return {Profile} the profile
 * @throws {Error} kind "forensicgenetics" if a locus appears twice
 */
export func profile(id as string, genotypes as list of Genotype) {
    def out as Profile init Profile{id: $id, loci: {}};
    for (def g in $genotypes) {
        if (maps.has($out.loci, $g.locus)) {
            fail("profile '" + $id + "' lists locus '" + $g.locus + "' twice");
        }
        $out.loci[$g.locus] = $g;
    }
    return $out;
}

/**
 * A copy of `p` with one genotype added, replacing any genotype already at
 * that locus.
 * @param p {Profile} the profile to copy
 * @param g {Genotype} the genotype to set
 * @return {Profile} the updated copy
 */
export func withGenotype(p as Profile, g as Genotype) {
    def out as Profile init $p;
    $out.loci[$g.locus] = $g;
    return $out;
}

/**
 * The loci a profile carries, in insertion order.
 * @param p {Profile} the profile to read
 * @return {list of string} the locus names
 */
export func profileLoci(p as Profile) {
    return maps.keys($p.loci);
}

/**
 * Whether a profile has a genotype at a locus.
 * @param p {Profile} the profile to test
 * @param locus {string} the locus name
 * @return {bool} true if typed there
 */
export func hasGenotype(p as Profile, locus as string) {
    return maps.has($p.loci, $locus);
}

/**
 * The genotype a profile records at a locus.
 * @param p {Profile} the profile to read
 * @param locus {string} the locus name
 * @return {Genotype} the genotype
 * @throws {Error} kind "forensicgenetics" if the profile is not typed at that locus
 */
export func genotypeAt(p as Profile, locus as string) {
    if (not maps.has($p.loci, $locus)) {
        fail("profile '" + $p.id + "' has no genotype at locus '" + $locus + "'");
    }
    return $p.loci[$locus];
}

/**
 * The loci a profile and a database have in common, in profile order.
 *
 * Pass the result to `restrict` to narrow a profile to a table's panel;
 * `randomMatchProbability` throws on an uncovered locus instead of dropping
 * it.
 * @param fdb {FrequencyDb} the database
 * @param p {Profile} the profile to compare
 * @return {list of string} the shared locus names
 */
export func commonLoci(fdb as FrequencyDb, p as Profile) {
    def out as list of string init [];
    for (def locus in $p.loci) {
        if (maps.has($fdb.loci, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * A copy of `p` keeping only the named loci.
 * @param p {Profile} the profile to narrow
 * @param keep {list of string} the locus names to retain
 * @return {Profile} the restricted profile
 * @throws {Error} kind "forensicgenetics" if the profile is not typed at a named locus
 */
export func restrict(p as Profile, keep as list of string) {
    def out as Profile init Profile{id: $p.id, loci: {}};
    for (def locus in $keep) {
        $out.loci[$locus] = genotypeAt($p, $locus);
    }
    return $out;
}

# --- match probabilities ----------------------------------------------------

# checkTheta rejects a theta outside [0, 1). At theta = 1 the Balding-Nichols
# denominator still divides, but every allele is identical by descent and the
# formula is meaningless, so the interval is half-open.
func checkTheta(theta as float) {
    if ($theta < 0.0 or $theta >= 1.0) {
        fail("theta must be in [0, 1)");
    }
}

/**
 * The single-locus genotype frequency with the NRC II theta correction.
 *
 * NRC II recommendation 4.2 (the Balding-Nichols sampling formula),
 * conditioned on one copy of the genotype already observed: the match
 * probability, not the unconditional population frequency.
 *
 * ```
 *   homozygote  aa: [2t + (1-t)pa] [3t + (1-t)pa] / [(1+t)(1+2t)]
 *   heterozygote ab: 2 [t + (1-t)pa] [t + (1-t)pb] / [(1+t)(1+2t)]
 * ```
 *
 * `pa` and `pb` are floored at 5/(2N). At `theta == 0` these reduce to the
 * Hardy-Weinberg `pa^2` and `2 pa pb`. The result is capped at 1.0, which a
 * floor above the listed frequencies on a small locus can otherwise exceed.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param g {Genotype} the genotype to score
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the genotype frequency, at most 1.0
 * @throws {Error} kind "forensicgenetics" on an out-of-range theta or a locus the
 *   database does not cover
 */
export func genotypeFrequency(fdb as FrequencyDb, g as Genotype, theta as float) {
    checkTheta($theta);
    def t as float init $theta;
    def den as float init (1.0 + $t) * (1.0 + 2.0 * $t);
    def pa as float init alleleFrequency($fdb, $g.locus, $g.a);
    if (isHomozygous($g)) {
        def x as float init 2.0 * $t + (1.0 - $t) * $pa;
        def y as float init 3.0 * $t + (1.0 - $t) * $pa;
        return atMostOne(($x * $y) / $den);
    }
    def pb as float init alleleFrequency($fdb, $g.locus, $g.b);
    def u as float init $t + (1.0 - $t) * $pa;
    def v as float init $t + (1.0 - $t) * $pb;
    return atMostOne((2.0 * $u * $v) / $den);
}

# atMostOne caps a genotype frequency at 1. A small locus whose 5/2N floor
# exceeds its own allele frequencies can push the product past 1, and a
# probability above 1 must never leave this module.
func atMostOne(p as float) {
    if ($p > 1.0) {
        return 1.0;
    }
    return $p;
}

/**
 * The random match probability (RMP) of a single-source profile.
 *
 * The product of `genotypeFrequency` over every locus of the profile: the
 * probability that an unrelated person from the reference population carries
 * it. Every locus must be in the database; narrow the profile with
 * `commonLoci` and `restrict` first if it is not. On a long panel the product
 * can underflow; use `log10RandomMatchProbability`.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param p {Profile} the profile to score
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the random match probability
 * @throws {Error} kind "forensicgenetics" on an empty profile, an out-of-range theta,
 *   or a locus the database does not cover
 */
export func randomMatchProbability(fdb as FrequencyDb, p as Profile, theta as float) {
    if (len($p.loci) == 0) {
        fail("profile '" + $p.id + "' has no genotypes");
    }
    def out as float init 1.0;
    for (def locus in $p.loci) {
        $out = $out * genotypeFrequency($fdb, $p.loci[$locus], $theta);
    }
    return $out;
}

/**
 * The base-10 logarithm of the random match probability.
 *
 * Sums per-locus logarithms, so it does not underflow where the product would.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param p {Profile} the profile to score
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} log10 of the random match probability (a negative number)
 * @throws {Error} kind "forensicgenetics" on an empty profile, an out-of-range theta,
 *   or a locus the database does not cover
 */
export func log10RandomMatchProbability(fdb as FrequencyDb, p as Profile, theta as float) {
    if (len($p.loci) == 0) {
        fail("profile '" + $p.id + "' has no genotypes");
    }
    def out as float init 0.0;
    for (def locus in $p.loci) {
        $out = $out + math.log10(genotypeFrequency($fdb, $p.loci[$locus], $theta));
    }
    return $out;
}

/**
 * The single-source likelihood ratio, `1 / RMP`.
 *
 * "The stain came from this person" against "it came from an unrelated person
 * from the reference population".
 * @param fdb {FrequencyDb} the reference frequencies
 * @param p {Profile} the profile to score
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "forensicgenetics" as for `randomMatchProbability`
 */
export func matchLikelihoodRatio(fdb as FrequencyDb, p as Profile, theta as float) {
    def rmp as float init randomMatchProbability($fdb, $p, $theta);
    if ($rmp <= 0.0) {
        fail("profile '" + $p.id + "' has a zero random match probability");
    }
    return 1.0 / $rmp;
}

/**
 * The random match probability when the donor's population is uncertain.
 *
 * The law of total probability over the candidate populations:
 * `P(G) = sum(w_i * P(G | population i))`, each term `randomMatchProbability`
 * in its own table. A table pooled with `combineFrequencies` averages allele
 * frequencies instead of genotype probabilities and gives a smaller number.
 *
 * The weights are the priors that the donor came from each population. They
 * must be non-negative and sum to 1 within 0.001. Every table must cover every
 * locus of the profile.
 * @param dbs {list of FrequencyDb} the candidate populations; at least one
 * @param weights {list of float} the prior per population, in `dbs` order
 * @param p {Profile} the profile to score
 * @param theta {float} the co-ancestry coefficient, applied within each
 *   population
 * @return {float} the weighted random match probability
 * @throws {Error} kind "forensicgenetics" on an empty list, a length mismatch,
 *   a negative weight, weights that do not sum to 1.0, or any error of
 *   `randomMatchProbability`
 */
export func randomMatchProbabilityAcross(
    dbs as list of FrequencyDb,
    weights as list of float,
    p as Profile,
    theta as float) {
    def nd as int init len($dbs);
    def nw as int init len($weights);
    if ($nd == 0) {
        fail("randomMatchProbabilityAcross needs at least one database");
    }
    if ($nd != $nw) {
        fail("got {$nd} databases but {$nw} weights");
    }
    def total as float init 0.0;
    for (def w in $weights) {
        if ($w < 0.0) {
            fail("a population weight must not be negative");
        }
        $total = $total + $w;
    }
    if (math.abs($total - 1.0) > 0.001) {
        fail("population weights sum to {$total}, not 1.0");
    }
    def out as float init 0.0;
    def i as int init 0;
    while ($i < $nd) {
        $out = $out + $weights[$i] * randomMatchProbability($dbs[$i], $p, $theta);
        $i = $i + 1;
    }
    return $out;
}

/**
 * The single-source likelihood ratio when the donor's population is uncertain:
 * `1 / randomMatchProbabilityAcross(...)`.
 * @param dbs {list of FrequencyDb} the candidate populations
 * @param weights {list of float} the prior per population, summing to 1.0
 * @param p {Profile} the profile to score
 * @param theta {float} the co-ancestry coefficient, in [0, 1)
 * @return {float} the likelihood ratio
 * @throws {Error} kind "forensicgenetics" as `randomMatchProbabilityAcross`,
 *   or if the weighted probability is zero
 */
export func matchLikelihoodRatioAcross(
    dbs as list of FrequencyDb,
    weights as list of float,
    p as Profile,
    theta as float) {
    def rmp as float init randomMatchProbabilityAcross($dbs, $weights, $p, $theta);
    if ($rmp <= 0.0) {
        fail("profile '" + $p.id + "' has a zero random match probability");
    }
    return 1.0 / $rmp;
}

# --- mixtures: inclusion / exclusion ----------------------------------------

/**
 * A mixture from the alleles observed at each locus. Repeated alleles in a
 * list are collapsed.
 * @param id {string} a label for the mixture
 * @param observed {map of string to list of string} locus -> observed alleles
 * @return {Mixture} the mixture
 * @throws {Error} kind "forensicgenetics" if a locus lists no alleles
 */
export func mixture(id as string, observed as map of string to list of string) {
    def out as Mixture init Mixture{id: $id, loci: {}};
    for (def locus in $observed) {
        def seen as list of string init distinct($observed[$locus]);
        if (len($seen) == 0) {
            fail("mixture '" + $id + "' locus '" + $locus + "' lists no alleles");
        }
        $out.loci[$locus] = $seen;
    }
    return $out;
}

/**
 * The loci a mixture covers, in insertion order.
 * @param m {Mixture} the mixture to read
 * @return {list of string} the locus names
 */
export func mixtureLoci(m as Mixture) {
    return maps.keys($m.loci);
}

/**
 * The probability of inclusion (PI) at one locus.
 *
 * `(sum of the observed alleles' floored frequencies)^2`: the probability that
 * a random person carries only alleles the mixture shows. The sum is capped at
 * 1.0 before squaring. Repeated alleles count once. No theta: PI counts who
 * cannot be excluded, it is not a co-ancestry-corrected match probability.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param locus {string} the locus name
 * @param observed {list of string} the alleles seen at that locus
 * @return {float} the probability of inclusion
 * @throws {Error} kind "forensicgenetics" on an empty allele list or a locus the
 *   database does not cover
 */
export func inclusionProbability(fdb as FrequencyDb, locus as string, observed as list of string) {
    def seen as list of string init distinct($observed);
    if (len($seen) == 0) {
        fail("locus '" + $locus + "': no observed alleles");
    }
    def total as float init 0.0;
    for (def a in $seen) {
        $total = $total + alleleFrequency($fdb, $locus, $a);
    }
    if ($total > 1.0) {
        $total = 1.0;
    }
    return $total * $total;
}

/**
 * The probability of exclusion (PE) at one locus: `1 - PI`.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param locus {string} the locus name
 * @param observed {list of string} the alleles seen at that locus
 * @return {float} the probability of exclusion
 * @throws {Error} kind "forensicgenetics" as for `inclusionProbability`
 */
export func exclusionProbability(fdb as FrequencyDb, locus as string, observed as list of string) {
    return 1.0 - inclusionProbability($fdb, $locus, $observed);
}

/**
 * The combined probability of inclusion (CPI): the product of
 * `inclusionProbability` over the mixture's loci.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param m {Mixture} the mixture to score
 * @return {float} the combined probability of inclusion
 * @throws {Error} kind "forensicgenetics" on an empty mixture or a locus the database
 *   does not cover
 */
export func cpi(fdb as FrequencyDb, m as Mixture) {
    if (len($m.loci) == 0) {
        fail("mixture '" + $m.id + "' has no loci");
    }
    def out as float init 1.0;
    for (def locus in $m.loci) {
        $out = $out * inclusionProbability($fdb, $locus, $m.loci[$locus]);
    }
    return $out;
}

/**
 * The combined probability of exclusion (CPE): `1 - CPI`, the share of the
 * population excluded by at least one locus.
 * @param fdb {FrequencyDb} the reference frequencies
 * @param m {Mixture} the mixture to score
 * @return {float} the combined probability of exclusion
 * @throws {Error} kind "forensicgenetics" as for `cpi`
 */
export func cpe(fdb as FrequencyDb, m as Mixture) {
    return 1.0 - cpi($fdb, $m);
}

/**
 * Whether every allele of `p` at each shared locus appears in the mixture: the
 * non-exclusion test that CPI quantifies.
 *
 * Loci only one of the two carries are skipped. Sharing no locus at all
 * throws rather than returning a vacuous `true`.
 * @param p {Profile} the profile under test
 * @param m {Mixture} the mixture
 * @return {bool} false as soon as one locus excludes the profile
 * @throws {Error} kind "forensicgenetics" if the profile and mixture share no locus
 */
export func isIncluded(p as Profile, m as Mixture) {
    def compared as int init 0;
    for (def locus in $m.loci) {
        if (maps.has($p.loci, $locus)) {
            $compared = $compared + 1;
            if (not coveredBy($p.loci[$locus], $m.loci[$locus])) {
                return false;
            }
        }
    }
    if ($compared == 0) {
        fail("profile '" + $p.id + "' and mixture '" + $m.id + "' share no locus, so the" +
            " profile is neither included nor excluded - there is nothing to compare");
    }
    return true;
}

# coveredBy reports whether every allele of `g` appears among `observed`.
func coveredBy(g as Genotype, observed as list of string) {
    for (def a in alleles($g)) {
        if (not lists.contains($observed, $a)) {
            return false;
        }
    }
    return true;
}
