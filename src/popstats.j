# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Population statistics for a reference table: the forensic-efficiency
 * parameters a lab reports beside a frequency database ("PowerStats"-style),
 * and the observed-sample statistics that table was estimated from.
 *
 * From allele frequencies, per locus and assuming Hardy-Weinberg proportions:
 * expected heterozygosity (Nei), match probability and power of discrimination
 * (Fisher 1951, Jones 1972), polymorphism information content (Botstein et al.
 * 1980), typical paternity index (Brenner & Morris 1990) and the power of
 * exclusion for a trio and for a duo (Jamieson & Taylor 1997). `locusStats`
 * gathers them for one locus, `statsTable` for every locus of a database, and
 * `combinedStats` multiplies them across loci, with log10 forms for the
 * products that underflow.
 *
 * From observed genotypes (a `list of forensics.Genotype`, all at one locus):
 * observed heterozygosity, allele counting, the observed match probability and
 * Brenner-Morris power of exclusion, and the Hardy-Weinberg exact test - exact
 * for a biallelic locus (Wigginton et al. 2005), by a seeded permutation test
 * for any number of alleles (Guo & Thompson 1992). `withEstimatedLocus` turns a
 * sample into a `FrequencyDb` locus, so the frequency statistics apply to it.
 *
 * From samples of several populations: Weir & Cockerham (1984) F-statistics
 * (θ / F_ST, F_IT, f / F_IS) per locus and over loci, and the inbreeding
 * coefficient of one sample. From one sample typed at two loci: a seeded
 * permutation test of genotypic linkage disequilibrium (phase unknown), for
 * one pair or every pair of a profile set.
 *
 * Frequencies are used as published, without the casework 5/(2N) floor, and
 * rescaled to sum to 1: a table that rounds to 0.998 still describes a whole
 * locus. A locus summing below `MIN_FREQUENCY_SUM` is a partial table and
 * throws, because every statistic here sums over all alleles.
 * @module popstats
 * @example
 * import "@deckyard/forensicgenetics/" as forensics;
 * import "@deckyard/forensicgenetics/popstats.j" as popstats;
 *
 * def table as list of popstats.LocusStats init popstats.statsTable($db);
 * def all as popstats.CombinedStats init popstats.combinedStats($table);
 * def pd as float init $all.powerOfDiscrimination;
 *
 * def sample as list of forensics.Genotype init [
 *     forensics.genotype("TH01", "6", "9.3"),
 *     forensics.genotype("TH01", "9.3", "9.3")
 * ];
 * def p as float init popstats.hwePermutationTest($sample, popstats.DEFAULT_PERMUTATIONS, 1);
 */
use math;
use maps;
use lists;
use convert;

import "./forensicgenetics.j" as forensics;

/**
 * The smallest allele-frequency sum a locus may have for the statistics here.
 *
 * Rounding leaves a published locus a little off 1; a sum below this means
 * alleles are missing, and a sum over all alleles cannot be formed.
 */
export def const MIN_FREQUENCY_SUM as float init 0.95;

/**
 * A permutation count for `hwePermutationTest`: standard error about 0.002
 * at p = 0.05, `sqrt(p (1 - p) / 10000)`.
 */
export def const DEFAULT_PERMUTATIONS as int init 10000;

/**
 * The largest permutation count `hwePermutationTest` and the LD tests
 * accept.
 */
export def const MAX_PERMUTATIONS as int init 1000000;

# Park-Miller "minimal standard" generator (Park & Miller 1988) with the
# multiplier 48271 of Park, Miller & Stockmeyer (1993): x <- 48271 x mod (2^31 - 1). Its
# products stay below 2^47, so int64 arithmetic never overflows, and the
# sequence depends on nothing but the seed - unlike math.rand, whose state is
# global and shared with the caller.
def const RAND_MODULUS as int init 2147483647;
def const RAND_MULTIPLIER as int init 48271;

# Two probabilities within this relative distance count as equal in a p-value
# sum, so a table exactly as likely as the observed one is not lost to
# rounding.
def const TIE_TOLERANCE as float init 0.0000001;

/**
 * The forensic-efficiency parameters of one locus, from its allele
 * frequencies under Hardy-Weinberg proportions.
 * @field locus {string} the locus name
 * @field n {int} individuals typed, as the table states; 0 if unknown
 * @field alleleCount {int} the number of alleles the table lists
 * @field expectedHeterozygosity {float} `1 - Σ p²`
 * @field matchProbability {float} `Σ P(g)²` over every genotype g
 * @field powerOfDiscrimination {float} `1 - matchProbability`
 * @field pic {float} the polymorphism information content
 * @field typicalPaternityIndex {float} `1 / (2 Σ p²)`
 * @field powerOfExclusionTrio {float} the share of random men a
 *   mother-child-man trio excludes
 * @field powerOfExclusionDuo {float} the share of random men a child-man duo
 *   excludes
 */
export def struct LocusStats {
    locus as string,
    n as int,
    alleleCount as int,
    expectedHeterozygosity as float,
    matchProbability as float,
    powerOfDiscrimination as float,
    pic as float,
    typicalPaternityIndex as float,
    powerOfExclusionTrio as float,
    powerOfExclusionDuo as float
};

/**
 * Per-locus statistics combined across loci, assuming independent loci.
 * @field loci {int} the number of loci combined
 * @field matchProbability {float} `Π MP`
 * @field log10MatchProbability {float} `Σ log10 MP`, for when the product
 *   underflows
 * @field powerOfDiscrimination {float} `1 - Π MP`
 * @field powerOfExclusionTrio {float} `1 - Π (1 - PE_trio)`
 * @field log10NonExclusionTrio {float} `Σ log10 (1 - PE_trio)`
 * @field powerOfExclusionDuo {float} `1 - Π (1 - PE_duo)`
 * @field log10NonExclusionDuo {float} `Σ log10 (1 - PE_duo)`
 * @field typicalPaternityIndex {float} `Π TPI`
 * @field log10TypicalPaternityIndex {float} `Σ log10 TPI`
 */
export def struct CombinedStats {
    loci as int,
    matchProbability as float,
    log10MatchProbability as float,
    powerOfDiscrimination as float,
    powerOfExclusionTrio as float,
    log10NonExclusionTrio as float,
    powerOfExclusionDuo as float,
    log10NonExclusionDuo as float,
    typicalPaternityIndex as float,
    log10TypicalPaternityIndex as float
};

# Sums holds the power sums S_k = Σ p^k of one locus's rescaled frequencies,
# from which every closed form below is written.
def struct Sums {
    k as int,
    s2 as float,
    s3 as float,
    s4 as float,
    s5 as float
};

# Sample is an observed genotype list reduced to integers: each allele by its
# index in first-seen order, two indices per individual.
def struct Sample {
    locus as string,
    labels as list of string,
    first as list of int,
    second as list of int
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "popstats", message: $msg, file: "", line: 0, col: 0};
}

# --- frequencies ------------------------------------------------------------

# locusFrequencies fetches one locus, with a message naming what is missing.
func locusFrequencies(fdb as forensics.FrequencyDb, locus as string) {
    if (not maps.has($fdb.loci, $locus)) {
        fail("frequency database '" + $fdb.name + "' has no locus '" + $locus + "'");
    }
    return $fdb.loci[$locus];
}

# powerSums rescales a locus's frequencies to sum to 1 and returns their power
# sums. It refuses an empty locus and one summing below MIN_FREQUENCY_SUM;
# withLocus has already refused one above MAX_FREQUENCY_SUM.
func powerSums(lf as forensics.LocusFrequencies) {
    def total as float init 0.0;
    for (def a in $lf.alleles) {
        $total = $total + $lf.alleles[$a];
    }
    if (len($lf.alleles) == 0) {
        fail("locus '" + $lf.locus + "' lists no alleles");
    }
    if ($total < MIN_FREQUENCY_SUM) {
        fail("locus '" + $lf.locus + "': allele frequencies sum to {$total}, below " +
            "MIN_FREQUENCY_SUM - a partial table has no heterozygosity or power to report");
    }
    def out as Sums init Sums{k: len($lf.alleles), s2: 0.0, s3: 0.0, s4: 0.0, s5: 0.0};
    for (def a in $lf.alleles) {
        def p as float init $lf.alleles[$a] / $total;
        def p2 as float init $p * $p;
        $out.s2 = $out.s2 + $p2;
        $out.s3 = $out.s3 + $p2 * $p;
        $out.s4 = $out.s4 + $p2 * $p2;
        $out.s5 = $out.s5 + $p2 * $p2 * $p;
    }
    return $out;
}

func sumsAt(fdb as forensics.FrequencyDb, locus as string) {
    return powerSums(locusFrequencies($fdb, $locus));
}

# matchFromSums: Σ over genotypes of P(g)² is Σ p⁴ + Σ_{i<j} (2 p_i p_j)²,
# and Σ_{i<j} p_i² p_j² = (S2² - S4) / 2, so MP = 2 S2² - S4.
func matchFromSums(s as Sums) {
    return 2.0 * $s.s2 * $s.s2 - $s.s4;
}

# picFromSums: Botstein et al. (1980), p. 320,
# 1 - Σ p² - Σ_{i<j} 2 p_i² p_j² = 1 - S2 - (S2² - S4).
func picFromSums(s as Sums) {
    return 1.0 - $s.s2 - ($s.s2 * $s.s2 - $s.s4);
}

# trioFromSums: Jamieson & Taylor (1997), the exclusion probability when the
# other parent's genotype is known.
func trioFromSums(s as Sums) {
    return 1.0 - 2.0 * $s.s2 + $s.s3 + 2.0 * $s.s4 - 3.0 * $s.s5 - 2.0 * $s.s2 * $s.s2 +
        3.0 * $s.s2 * $s.s3;
}

# duoFromSums: Jamieson & Taylor (1997), the exclusion probability when the
# other parent is untyped.
func duoFromSums(s as Sums) {
    return 1.0 - 4.0 * $s.s2 + 2.0 * $s.s2 * $s.s2 + 4.0 * $s.s3 - 3.0 * $s.s4;
}

/**
 * The expected heterozygosity of a locus, `1 - Σ p²` (Nei 1973, gene
 * diversity).
 *
 * The probability that two alleles drawn at random from the population differ.
 * This is the plug-in form; `unbiasedHeterozygosity` corrects it for the
 * sample the frequencies were estimated from.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the expected heterozygosity, in [0, 1)
 * @throws {Error} kind "popstats" if the table has no such locus, or its
 *   frequencies sum below `MIN_FREQUENCY_SUM`
 */
export func expectedHeterozygosity(fdb as forensics.FrequencyDb, locus as string) {
    return 1.0 - sumsAt($fdb, $locus).s2;
}

/**
 * The expected heterozygosity with Nei's small-sample correction,
 * `2N / (2N - 1) · (1 - Σ p²)` (Nei 1978, eq. 2; Nei & Roychoudhury 1974).
 *
 * N is the locus's sample size in individuals, so 2N alleles were drawn; the
 * plug-in `1 - Σ p²` underestimates the population value by the factor
 * `(2N - 1) / 2N`.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the unbiased expected heterozygosity
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`, or if the locus
 *   has no sample size (`n == 0`)
 */
export func unbiasedHeterozygosity(fdb as forensics.FrequencyDb, locus as string) {
    def lf as forensics.LocusFrequencies init locusFrequencies($fdb, $locus);
    if ($lf.n <= 0) {
        fail("locus '" + $locus + "' has no sample size, so the small-sample " +
            "correction cannot be applied; use expectedHeterozygosity");
    }
    def twoN as float init 2.0 * convert.toFloat($lf.n);
    return $twoN / ($twoN - 1.0) * (1.0 - powerSums($lf).s2);
}

/**
 * The match probability of a locus: the chance that two unrelated people
 * share a genotype, `Σ P(g)²` over every genotype with Hardy-Weinberg
 * frequencies, which is `2 (Σ p²)² - Σ p⁴` (Fisher 1951; Jones 1972).
 *
 * Not the random match probability of a given profile (`forensics.
 * randomMatchProbability`); this averages over all profiles.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the match probability, in (0, 1]
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func matchProbability(fdb as forensics.FrequencyDb, locus as string) {
    return matchFromSums(sumsAt($fdb, $locus));
}

/**
 * The power of discrimination of a locus, `1 - matchProbability` (Fisher
 * 1951; Jones 1972).
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the power of discrimination, in [0, 1)
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func powerOfDiscrimination(fdb as forensics.FrequencyDb, locus as string) {
    return 1.0 - matchFromSums(sumsAt($fdb, $locus));
}

/**
 * The polymorphism information content of a locus (Botstein et al. 1980):
 * `1 - Σ p_i² - Σ_{i<j} 2 p_i² p_j²`.
 *
 * The chance that an offspring's allele reveals which of a heterozygous
 * parent's alleles it received; always below `expectedHeterozygosity`.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the PIC, in [0, 1)
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func polymorphismInformationContent(fdb as forensics.FrequencyDb, locus as string) {
    return picFromSums(sumsAt($fdb, $locus));
}

/**
 * The typical paternity index of a locus from its expected homozygosity,
 * `1 / (2 Σ p²)` (Brenner & Morris 1990, eq. V.C.3, `TPI = (H + h) / 2h` with `h` the
 * homozygosity and `H + h = 1`).
 *
 * `observedTypicalPaternityIndex` is the same formula on a sample's observed
 * homozygosity, the value PowerStats reports.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the typical paternity index, at least 0.5
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func typicalPaternityIndex(fdb as forensics.FrequencyDb, locus as string) {
    return 1.0 / (2.0 * sumsAt($fdb, $locus).s2);
}

/**
 * The power of exclusion of a locus in a trio: the share of random men a
 * locus excludes when the mother's genotype is known (Jamieson & Taylor 1997,
 * eq. 1a):
 *
 * ```
 * 1 - 2 S2 + S3 + 2 S4 - 3 S5 - 2 S2² + 3 S2 S3        Sk = Σ p^k
 * ```
 *
 * Exact under Hardy-Weinberg for any number of alleles; 3/16 for two equally
 * frequent alleles.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the trio power of exclusion
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func powerOfExclusionTrio(fdb as forensics.FrequencyDb, locus as string) {
    return trioFromSums(sumsAt($fdb, $locus));
}

/**
 * The power of exclusion of a locus in a duo: the share of random men a locus
 * excludes when only the child and the man are typed (Jamieson & Taylor 1997,
 * eq. 2a, after Garber & Morris 1983):
 *
 * ```
 * 1 - 4 S2 + 2 S2² + 4 S3 - 3 S4                      Sk = Σ p^k
 * ```
 *
 * A man is excluded only when he shares no allele with the child, so this is
 * always below the trio value.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {float} the duo power of exclusion
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func powerOfExclusionDuo(fdb as forensics.FrequencyDb, locus as string) {
    return duoFromSums(sumsAt($fdb, $locus));
}

/**
 * Every frequency statistic of one locus, from a single pass over its
 * alleles.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @param locus {string} the locus name
 * @return {LocusStats} the locus's statistics
 * @throws {Error} kind "popstats" as `expectedHeterozygosity`
 */
export func locusStats(fdb as forensics.FrequencyDb, locus as string) {
    def lf as forensics.LocusFrequencies init locusFrequencies($fdb, $locus);
    return statsOf($lf);
}

# statsOf is locusStats once the locus is in hand.
func statsOf(lf as forensics.LocusFrequencies) {
    def s as Sums init powerSums($lf);
    def mp as float init matchFromSums($s);
    return LocusStats{
        locus: $lf.locus,
        n: $lf.n,
        alleleCount: $s.k,
        expectedHeterozygosity: 1.0 - $s.s2,
        matchProbability: $mp,
        powerOfDiscrimination: 1.0 - $mp,
        pic: picFromSums($s),
        typicalPaternityIndex: 1.0 / (2.0 * $s.s2),
        powerOfExclusionTrio: trioFromSums($s),
        powerOfExclusionDuo: duoFromSums($s)
    };
}

/**
 * `locusStats` for every locus of a table, in the table's locus order.
 * @param fdb {forensics.FrequencyDb} the frequency table
 * @return {list of LocusStats} one row per locus
 * @throws {Error} kind "popstats" if the table has no loci, or a locus sums
 *   below `MIN_FREQUENCY_SUM`
 */
export func statsTable(fdb as forensics.FrequencyDb) {
    if (len($fdb.loci) == 0) {
        fail("frequency database '" + $fdb.name + "' has no loci");
    }
    def out as list of LocusStats init [];
    for (def locus in $fdb.loci) {
        $out[] = statsOf($fdb.loci[$locus]);
    }
    return $out;
}

# --- combining loci ---------------------------------------------------------

# checkUnit refuses a per-locus probability outside [0, 1].
func checkUnit(x as float) {
    if ($x < 0.0 or $x > 1.0) {
        fail("a power must be in [0, 1], got {$x}");
    }
}

/**
 * Per-locus statistics combined across loci, as products over independent
 * loci: `Π MP`, `1 - Π MP`, `1 - Π (1 - PE)` and `Π TPI`, each product also as
 * a sum of log10 terms.
 * @param table {list of LocusStats} the loci to combine, e.g. from `statsTable`
 * @return {CombinedStats} the combined values
 * @throws {Error} kind "popstats" on an empty table
 */
export func combinedStats(table as list of LocusStats) {
    if (len($table) == 0) {
        fail("combinedStats needs at least one locus");
    }
    def out as CombinedStats init CombinedStats{
        loci: len($table),
        matchProbability: 1.0,
        log10MatchProbability: 0.0,
        powerOfDiscrimination: 0.0,
        powerOfExclusionTrio: 0.0,
        log10NonExclusionTrio: 0.0,
        powerOfExclusionDuo: 0.0,
        log10NonExclusionDuo: 0.0,
        typicalPaternityIndex: 1.0,
        log10TypicalPaternityIndex: 0.0
    };
    def keepTrio as float init 1.0;
    def keepDuo as float init 1.0;
    for (def i in 0..len($table)) {
        def row as LocusStats init $table[$i];
        $out.matchProbability = $out.matchProbability * $row.matchProbability;
        $out.log10MatchProbability = $out.log10MatchProbability + math.log10($row.matchProbability);
        $keepTrio = $keepTrio * (1.0 - $row.powerOfExclusionTrio);
        $out.log10NonExclusionTrio = $out.log10NonExclusionTrio +
            math.log10(1.0 - $row.powerOfExclusionTrio);
        $keepDuo = $keepDuo * (1.0 - $row.powerOfExclusionDuo);
        $out.log10NonExclusionDuo = $out.log10NonExclusionDuo +
            math.log10(1.0 - $row.powerOfExclusionDuo);
        $out.typicalPaternityIndex = $out.typicalPaternityIndex * $row.typicalPaternityIndex;
        $out.log10TypicalPaternityIndex = $out.log10TypicalPaternityIndex +
            math.log10($row.typicalPaternityIndex);
    }
    $out.powerOfDiscrimination = 1.0 - $out.matchProbability;
    $out.powerOfExclusionTrio = 1.0 - $keepTrio;
    $out.powerOfExclusionDuo = 1.0 - $keepDuo;
    return $out;
}

/**
 * Per-locus powers (of discrimination, of exclusion) combined over
 * independent loci: `1 - Π (1 - x)`.
 *
 * For values that do not come from a `LocusStats`, such as
 * `observedPowerOfExclusion` per locus.
 * @param powers {list of float} one power per locus, each in [0, 1]
 * @return {float} the combined power
 * @throws {Error} kind "popstats" on an empty list or a value outside [0, 1]
 */
export func combinedPower(powers as list of float) {
    if (len($powers) == 0) {
        fail("combinedPower needs at least one locus");
    }
    def keep as float init 1.0;
    for (def x in $powers) {
        checkUnit($x);
        $keep = $keep * (1.0 - $x);
    }
    return 1.0 - $keep;
}

/**
 * The base-10 logarithm of `1 - combinedPower(powers)`, `Σ log10 (1 - x)`:
 * the combined non-discrimination or non-exclusion probability, which
 * underflows as a product on a long panel.
 * @param powers {list of float} one power per locus, each in [0, 1)
 * @return {float} log10 of the combined complement (zero or negative)
 * @throws {Error} kind "popstats" on an empty list or a value outside [0, 1)
 */
export func log10CombinedComplement(powers as list of float) {
    if (len($powers) == 0) {
        fail("log10CombinedComplement needs at least one locus");
    }
    def out as float init 0.0;
    for (def x in $powers) {
        checkUnit($x);
        if ($x >= 1.0) {
            fail("a power of 1 has no finite log10 complement");
        }
        $out = $out + math.log10(1.0 - $x);
    }
    return $out;
}

# --- observed genotypes -----------------------------------------------------

# encode checks a sample (non-empty, one locus, no empty allele) and indexes
# its alleles in first-seen order.
func encode(genotypes as list of forensics.Genotype) {
    if (len($genotypes) == 0) {
        fail("the sample has no genotypes");
    }
    def out as Sample init Sample{locus: $genotypes[0].locus, labels: [], first: [], second: []};
    def index as map of string to int init {};
    for (def g in $genotypes) {
        if ($g.locus != $out.locus) {
            fail("the sample mixes loci '" + $out.locus + "' and '" + $g.locus + "'");
        }
        if (len($g.a) == 0 or len($g.b) == 0) {
            fail("locus '" + $out.locus + "': allele labels must not be empty");
        }
        for (def a in [$g.a, $g.b]) {
            if (not maps.has($index, $a)) {
                $index[$a] = len($out.labels);
                $out.labels[] = $a;
            }
        }
        $out.first[] = $index[$g.a];
        $out.second[] = $index[$g.b];
    }
    return $out;
}

/**
 * The share of a sample that is heterozygous, `H_obs = (heterozygotes) / n`.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} the observed heterozygosity, in [0, 1]
 * @throws {Error} kind "popstats" on an empty sample, mixed loci or an empty
 *   allele label
 */
export func observedHeterozygosity(genotypes as list of forensics.Genotype) {
    def s as Sample init encode($genotypes);
    return convert.toFloat(heterozygotes($s)) / convert.toFloat(len($s.first));
}

func heterozygotes(s as Sample) {
    def out as int init 0;
    for (def i in 0..len($s.first)) {
        if ($s.first[$i] != $s.second[$i]) {
            $out = $out + 1;
        }
    }
    return $out;
}

/**
 * How many copies of each allele a sample carries: a homozygote counts its
 * allele twice. The counts sum to 2n.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {map of string to int} allele label -> copies, in first-seen order
 * @throws {Error} kind "popstats" as `observedHeterozygosity`
 */
export func alleleCounts(genotypes as list of forensics.Genotype) {
    def s as Sample init encode($genotypes);
    def out as map of string to int init {};
    for (def a in $s.labels) {
        $out[$a] = 0;
    }
    for (def i in 0..len($s.first)) {
        def x as string init $s.labels[$s.first[$i]];
        def y as string init $s.labels[$s.second[$i]];
        $out[$x] = $out[$x] + 1;
        $out[$y] = $out[$y] + 1;
    }
    return $out;
}

/**
 * Allele frequencies estimated from a sample by gene counting: each allele's
 * copies over the 2n alleles sampled (the maximum-likelihood estimate for
 * codominant alleles).
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {map of string to float} allele label -> frequency; sums to 1
 * @throws {Error} kind "popstats" as `observedHeterozygosity`
 */
export func estimateFrequencies(genotypes as list of forensics.Genotype) {
    def counts as map of string to int init alleleCounts($genotypes);
    def twoN as float init 2.0 * convert.toFloat(len($genotypes));
    def out as map of string to float init {};
    for (def a in $counts) {
        $out[$a] = convert.toFloat($counts[$a]) / $twoN;
    }
    return $out;
}

/**
 * A copy of `fdb` with the sample's locus added or replaced: the counted
 * frequencies with `n` set to the sample size, so the 5/(2N) floor and
 * `unbiasedHeterozygosity` see the right N.
 * @param fdb {forensics.FrequencyDb} the table to extend
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {forensics.FrequencyDb} the extended copy
 * @throws {Error} kind "popstats" as `observedHeterozygosity`
 */
export func withEstimatedLocus(
    fdb as forensics.FrequencyDb,
    genotypes as list of forensics.Genotype) {
    def freqs as map of string to float init estimateFrequencies($genotypes);
    return forensics.withLocus($fdb, $genotypes[0].locus, len($genotypes), $freqs);
}

# genotypeCounts counts each distinct genotype in a sample, keyed by its two
# allele indices in ascending order.
func genotypeCounts(s as Sample) {
    def out as map of string to int init {};
    for (def i in 0..len($s.first)) {
        def key as string init pairKey($s.first[$i], $s.second[$i]);
        if (maps.has($out, $key)) {
            $out[$key] = $out[$key] + 1;
        } else {
            $out[$key] = 1;
        }
    }
    return $out;
}

func pairKey(x as int, y as int) {
    if ($x > $y) {
        return "{$y}/{$x}";
    }
    return "{$x}/{$y}";
}

/**
 * The match probability observed in a sample: `Σ x_g²` over the genotypes
 * seen, `x_g` each genotype's share of the sample (Fisher 1951; the form
 * PowerStats reports). It does not assume Hardy-Weinberg proportions.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} the observed match probability, in (0, 1]
 * @throws {Error} kind "popstats" as `observedHeterozygosity`
 */
export func observedMatchProbability(genotypes as list of forensics.Genotype) {
    def counts as map of string to int init genotypeCounts(encode($genotypes));
    def n as float init convert.toFloat(len($genotypes));
    def out as float init 0.0;
    for (def key in $counts) {
        def x as float init convert.toFloat($counts[$key]) / $n;
        $out = $out + $x * $x;
    }
    return $out;
}

/**
 * The typical paternity index from a sample's observed homozygosity `h`,
 * `1 / (2h)` (Brenner & Morris 1990, eq. V.C.3).
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} the typical paternity index
 * @throws {Error} kind "popstats" as `observedHeterozygosity`, or if the
 *   sample has no homozygote (the index is unbounded)
 */
export func observedTypicalPaternityIndex(genotypes as list of forensics.Genotype) {
    def h as float init 1.0 - observedHeterozygosity($genotypes);
    if ($h <= 0.0) {
        fail("the sample has no homozygote, so its typical paternity index is unbounded");
    }
    return 1.0 / (2.0 * $h);
}

/**
 * The power of exclusion from a sample's observed heterozygosity `H` and
 * homozygosity `h = 1 - H`, `H² (1 - 2 H h²)` (Brenner & Morris 1990, eq. V.B.2; the
 * value PowerStats reports).
 *
 * An approximation that uses only the heterozygosity. It equals the exact
 * trio value, `powerOfExclusionTrio`, for a biallelic locus in Hardy-Weinberg
 * proportions with equal frequencies (3/16), and is close for typical STRs.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} the power of exclusion
 * @throws {Error} kind "popstats" as `observedHeterozygosity`
 */
export func observedPowerOfExclusion(genotypes as list of forensics.Genotype) {
    def het as float init observedHeterozygosity($genotypes);
    def hom as float init 1.0 - $het;
    return $het * $het * (1.0 - 2.0 * $het * $hom * $hom);
}

# --- Hardy-Weinberg tests ---------------------------------------------------

/**
 * The Hardy-Weinberg exact test p-value of a biallelic sample (Wigginton,
 * Cutler & Abecasis 2005).
 *
 * Conditional on the allele counts, the probability of each heterozygote
 * count is the Levene (1949) / Haldane (1954) distribution; Wigginton et al.
 * evaluate it by a recurrence from its mode. The p-value sums the
 * probabilities of every heterozygote count no more likely than the one
 * observed (ties within a relative 1e-7 count as equal). A monomorphic sample
 * returns 1.0. For three or more alleles use `hwePermutationTest`.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} the two-sided exact p-value, in (0, 1]
 * @throws {Error} kind "popstats" as `observedHeterozygosity`, or if the
 *   sample shows more than two alleles
 */
export func hweExactTest(genotypes as list of forensics.Genotype) {
    def s as Sample init encode($genotypes);
    def k as int init len($s.labels);
    if ($k == 1) {
        return 1.0;
    }
    if ($k > 2) {
        fail("locus '" + $s.locus + "' shows {$k} alleles; the exact test here is " +
            "biallelic - use hwePermutationTest");
    }
    def het as int init heterozygotes($s);
    def hom0 as int init 0;
    for (def i in 0..len($s.first)) {
        if ($s.first[$i] == 0 and $s.second[$i] == 0) {
            $hom0 = $hom0 + 1;
        }
    }
    def hom1 as int init len($s.first) - $het - $hom0;
    return wiggintonP($het, math.min($hom0, $hom1), math.max($hom0, $hom1));
}

# wiggintonP is the SNP-HWE algorithm of Wigginton et al. (2005): the
# heterozygote-count distribution built outward from its mode by the
# recurrences
#   P(h - 2) = P(h) h (h - 1) / (4 (homr + 1) (homc + 1))
#   P(h + 2) = P(h) 4 homr homc / ((h + 2) (h + 1))
# then normalised. Every term is at most the mode's, so nothing overflows.
func wiggintonP(het as int, homr as int, homc as int) {
    def rare as int init 2 * $homr + $het;
    def n as int init $het + $homr + $homc;
    def probs as list of float init [];
    for (def i in 0..($rare + 1)) {
        $probs[] = 0.0;
    }
    def mid as int init ($rare * (2 * $n - $rare)) // (2 * $n);
    if ($mid % 2 != $rare % 2) {
        $mid = $mid + 1;
    }
    $probs[$mid] = 1.0;
    def total as float init 1.0;
    def h as int init $mid;
    def r as int init ($rare - $mid) // 2;
    def c as int init $n - $mid - $r;
    while ($h > 1) {
        def num as float init convert.toFloat($h * ($h - 1));
        $probs[$h - 2] = $probs[$h] * $num / convert.toFloat(4 * ($r + 1) * ($c + 1));
        $total = $total + $probs[$h - 2];
        $h = $h - 2;
        $r = $r + 1;
        $c = $c + 1;
    }
    $h = $mid;
    $r = ($rare - $mid) // 2;
    $c = $n - $mid - $r;
    while ($h <= $rare - 2) {
        def up as float init convert.toFloat(4 * $r * $c);
        $probs[$h + 2] = $probs[$h] * $up / convert.toFloat(($h + 2) * ($h + 1));
        $total = $total + $probs[$h + 2];
        $h = $h + 2;
        $r = $r - 1;
        $c = $c - 1;
    }
    def limit as float init $probs[$het] * (1.0 + TIE_TOLERANCE);
    def p as float init 0.0;
    for (def x in $probs) {
        if ($x > 0.0 and $x <= $limit) {
            $p = $p + $x;
        }
    }
    return math.min(1.0, $p / $total);
}

# tableScore is the part of log P(table | allele counts) that varies between
# tables with the same allele counts (Levene 1949; Guo & Thompson 1992):
#   P = n! Π n_i! 2^H / ((2n)! Π_{i<=j} n_ij!)
# so log P = const + H ln 2 - Σ ln n_ij!.
func tableScore(s as Sample) {
    def counts as map of string to int init genotypeCounts($s);
    def out as float init convert.toFloat(heterozygotes($s)) * math.ln(2.0);
    for (def key in $counts) {
        $out = $out - math.lgamma(convert.toFloat($counts[$key] + 1));
    }
    return $out;
}

/**
 * The Hardy-Weinberg test p-value of a sample with any number of alleles, by
 * permutation (Guo & Thompson 1992, their Monte Carlo test).
 *
 * The 2n sampled alleles are shuffled and re-paired `permutations` times,
 * which draws genotype tables from their exact distribution given the allele
 * counts. Each table is scored by its Levene-Haldane probability; the p-value
 * is `(r + 1) / (permutations + 1)`, `r` the number of shuffled tables no more
 * likely than the observed one (Phipson & Smyth 2010: never zero, and a valid
 * p-value at any permutation count). Its standard error is about
 * `sqrt(p (1 - p) / permutations)`.
 *
 * The shuffle (Fisher-Yates) draws from a Park-Miller generator seeded with
 * `seed`, so a given seed always gives the same answer and the caller's
 * `math.rand` state is untouched. On a biallelic sample it estimates
 * `hweExactTest`.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @param permutations {int} shuffles to draw, 1 to `MAX_PERMUTATIONS`
 * @param seed {int} the generator seed, 1 to 2147483646
 * @return {float} the permutation p-value, in (0, 1]
 * @throws {Error} kind "popstats" as `observedHeterozygosity`, or on a
 *   permutation count or seed out of range
 */
export func hwePermutationTest(
    genotypes as list of forensics.Genotype,
    permutations as int,
    seed as int) {
    if ($permutations < 1 or $permutations > MAX_PERMUTATIONS) {
        fail("permutations must be between 1 and 1000000, got {$permutations}");
    }
    if ($seed < 1 or $seed >= RAND_MODULUS) {
        fail("seed must be between 1 and 2147483646, got {$seed}");
    }
    def s as Sample init encode($genotypes);
    if (len($s.labels) == 1) {
        return 1.0;
    }
    def observed as float init tableScore($s) + math.log1p(TIE_TOLERANCE);
    def n as int init len($s.first);
    def k as int init len($s.labels);
    def pool as list of int init [];
    for (def i in 0..$n) {
        $pool[] = $s.first[$i];
        $pool[] = $s.second[$i];
    }
    def lnCount as list of float init lnTable($n);
    def accept as list of int init acceptLimits(2 * $n);
    # Genotype counts live in a flat k*k list; a cell is current only when its
    # stamp is this draw's number, so no pass is spent resetting it.
    def cells as list of int init [];
    def stamp as list of int init [];
    for (def i in 0..($k * $k)) {
        $cells[] = 0;
        $stamp[] = 0;
    }
    def ln2 as float init math.ln(2.0);
    def state as int init $seed;
    def hits as int init 0;
    def draw as int init 0;
    # Fisher-Yates from the top fixes positions j and j - 1 per pass, so the
    # pair (j - 1, j) is final and is scored in the same pass: going from c to
    # c + 1 copies of a genotype subtracts ln(c + 1) from the score.
    while ($draw < $permutations) {
        $draw = $draw + 1;
        def score as float init 0.0;
        def j as int init 2 * $n - 1;
        while ($j > 0) {
            $state = ($state * RAND_MULTIPLIER) % RAND_MODULUS;
            if ($state - 1 >= $accept[$j + 1]) {
                $state = redraw($state, $accept[$j + 1]);
            }
            def pick as int init ($state - 1) % ($j + 1);
            def x as int init $pool[$pick];
            $pool[$pick] = $pool[$j];
            $pool[$j] = $x;
            $state = ($state * RAND_MULTIPLIER) % RAND_MODULUS;
            if ($state - 1 >= $accept[$j]) {
                $state = redraw($state, $accept[$j]);
            }
            $pick = ($state - 1) % $j;
            def y as int init $pool[$pick];
            $pool[$pick] = $pool[$j - 1];
            $pool[$j - 1] = $y;
            def cell as int init $x * $k + $y;
            if ($x > $y) {
                $cell = $y * $k + $x;
            }
            if ($x != $y) {
                $score = $score + $ln2;
            }
            if ($stamp[$cell] != $draw) {
                $stamp[$cell] = $draw;
                $cells[$cell] = 0;
            }
            $cells[$cell] = $cells[$cell] + 1;
            $score = $score - $lnCount[$cells[$cell]];
            $j = $j - 2;
        }
        if ($score <= $observed) {
            $hits = $hits + 1;
        }
    }
    return convert.toFloat($hits + 1) / convert.toFloat($permutations + 1);
}

# lnTable lists ln(c) for c in 0..n, with ln(0) unused and set to 0.
func lnTable(n as int) {
    def out as list of float init [0.0];
    for (def c in 1..($n + 1)) {
        $out[] = math.ln(convert.toFloat($c));
    }
    return $out;
}

# acceptLimits lists, for each bound b in 0..m, the largest multiple of b not
# above the generator's 2^31 - 2 outputs; entry 0 is unused.
func acceptLimits(m as int) {
    def span as int init RAND_MODULUS - 1;
    def out as list of int init [0];
    for (def b in 1..($m + 1)) {
        $out[] = $span - $span % $b;
    }
    return $out;
}

# nextState is one Park-Miller step, x <- 48271 x mod (2^31 - 1).
func nextState(state as int) {
    return ($state * RAND_MULTIPLIER) % RAND_MODULUS;
}

# redraw steps the generator until (state - 1) < accept, the largest multiple
# of the draw's bound below the period, so `(state - 1) % bound` is exactly
# uniform: rejection rather than a modulo-biased draw. Rejection has
# probability below bound / 2^31, so this is rarely called.
func redraw(state as int, accept as int) {
    def x as int init nextState($state);
    while ($x - 1 >= $accept) {
        $x = nextState($x);
    }
    return $x;
}

# --- F-statistics -----------------------------------------------------------

/**
 * Weir & Cockerham (1984) variance components of one locus, or of several
 * loci summed, from which `fst`, `fit` and `fis` are ratios.
 *
 * `a` is the component between populations, `b` between individuals within
 * populations, `c` between gametes within individuals, each summed over
 * alleles (and over loci for `combineFStatistics`). Components are unbiased
 * estimates and may be negative.
 * @field locus {string} the locus name; "" for a combination of loci
 * @field loci {int} the number of loci summed
 * @field populations {int} the number of population samples, r
 * @field a {float} the between-population component
 * @field b {float} the between-individual, within-population component
 * @field c {float} the within-individual component
 */
export def struct FStatistics {
    locus as string,
    loci as int,
    populations as int,
    a as float,
    b as float,
    c as float
};

# A ratio whose denominator is within this of zero is undefined: every
# population is monomorphic for the same allele, or (for fis) every
# population is fixed.
def const DEGENERATE as float init 0.000000000001;

# Tally is one locus's population samples reduced to counts over one shared
# allele index: copies[i * k + u] and hets[i * k + u] are allele u's copies,
# and the heterozygotes carrying it, in sample i.
def struct Tally {
    locus as string,
    k as int,
    sizes as list of int,
    copies as list of int,
    hets as list of int
};

# tally encodes every sample (one locus each, as `encode` checks), checks they
# share the locus, and counts allele copies and heterozygotes per sample.
func tally(samples as list of list of forensics.Genotype) {
    if (len($samples) == 0) {
        fail("F-statistics need at least one population sample");
    }
    def encoded as list of Sample init [];
    def index as map of string to int init {};
    for (def i in 0..len($samples)) {
        def s as Sample init encode($samples[$i]);
        if ($i > 0 and $s.locus != $encoded[0].locus) {
            fail("population samples mix loci '" + $encoded[0].locus + "' and '" + $s.locus + "'");
        }
        for (def a in $s.labels) {
            if (not maps.has($index, $a)) {
                $index[$a] = len($index);
            }
        }
        $encoded[] = $s;
    }
    def k as int init len($index);
    def out as Tally init Tally{locus: $encoded[0].locus, k: $k, sizes: [], copies: [], hets: []};
    for (def i in 0..(len($samples) * $k)) {
        $out.copies[] = 0;
        $out.hets[] = 0;
    }
    for (def i in 0..len($encoded)) {
        # Map the sample's own first-seen indices onto the shared index.
        def local as list of int init [];
        for (def a in $encoded[$i].labels) {
            $local[] = $index[$a];
        }
        def base as int init $i * $k;
        def n as int init len($encoded[$i].first);
        for (def j in 0..$n) {
            def x as int init $base + $local[$encoded[$i].first[$j]];
            def y as int init $base + $local[$encoded[$i].second[$j]];
            $out.copies[$x] = $out.copies[$x] + 1;
            $out.copies[$y] = $out.copies[$y] + 1;
            if ($x != $y) {
                $out.hets[$x] = $out.hets[$x] + 1;
                $out.hets[$y] = $out.hets[$y] + 1;
            }
        }
        $out.sizes[] = $n;
    }
    return $out;
}

/**
 * The Weir & Cockerham (1984) variance components `a`, `b`, `c` of one
 * locus, from a sample of genotypes in each of r populations.
 *
 * Per allele, with `n_i` individuals in sample i, `p_i` the allele's
 * frequency there and `h_i` the share of individuals heterozygous for it:
 *
 * ```
 * n̄   = Σ n_i / r                 n_c = (r n̄ - Σ n_i² / (r n̄)) / (r - 1)
 * p̄   = Σ n_i p_i / (r n̄)         h̄   = Σ n_i h_i / (r n̄)
 * s²  = Σ n_i (p_i - p̄)² / ((r - 1) n̄)
 * a   = n̄ / n_c · [s² - (p̄ (1 - p̄) - (r - 1) s² / r - h̄ / 4) / (n̄ - 1)]
 * b   = n̄ / (n̄ - 1) · [p̄ (1 - p̄) - (r - 1) s² / r - (2 n̄ - 1) h̄ / (4 n̄)]
 * c   = h̄ / 2
 * ```
 *
 * summed over the alleles seen in any sample. With one sample (r = 1) the
 * `s²` terms vanish and `a` is 0: only `fis` is defined, the within-sample
 * inbreeding coefficient (`inbreedingCoefficient`). An individual untyped at
 * the locus is left out of its sample.
 * @param samples {list of list of forensics.Genotype} one sample per
 *   population, every genotype at one locus
 * @return {FStatistics} the locus's components
 * @throws {Error} kind "popstats" on no samples, an empty sample, mixed loci,
 *   an empty allele label, or a mean sample size of 1 or less
 */
export func fStatistics(samples as list of list of forensics.Genotype) {
    def t as Tally init tally($samples);
    def r as int init len($t.sizes);
    def total as int init 0;
    def squares as int init 0;
    for (def n in $t.sizes) {
        $total = $total + $n;
        $squares = $squares + $n * $n;
    }
    def rf as float init convert.toFloat($r);
    def nbar as float init convert.toFloat($total) / $rf;
    if ($nbar <= 1.0) {
        fail("locus '" + $t.locus + "': the mean sample size must exceed one individual");
    }
    def nc as float init 0.0;
    if ($r > 1) {
        $nc = (convert.toFloat($total) - convert.toFloat($squares) / convert.toFloat($total)) /
            ($rf - 1.0);
    }
    def out as FStatistics init FStatistics{
        locus: $t.locus,
        loci: 1,
        populations: $r,
        a: 0.0,
        b: 0.0,
        c: 0.0
    };
    for (def u in 0..$t.k) {
        def pbar as float init 0.0;
        def hbar as float init 0.0;
        for (def i in 0..$r) {
            # n_i p_i = copies / 2, n_i h_i = heterozygotes
            $pbar = $pbar + convert.toFloat($t.copies[$i * $t.k + $u]) / 2.0;
            $hbar = $hbar + convert.toFloat($t.hets[$i * $t.k + $u]);
        }
        $pbar = $pbar / convert.toFloat($total);
        $hbar = $hbar / convert.toFloat($total);
        def s2 as float init 0.0;
        if ($r > 1) {
            for (def i in 0..$r) {
                def ni as float init convert.toFloat($t.sizes[$i]);
                def copies as float init convert.toFloat($t.copies[$i * $t.k + $u]);
                def d as float init $copies / (2.0 * $ni) - $pbar;
                $s2 = $s2 + $ni * $d * $d;
            }
            $s2 = $s2 / (($rf - 1.0) * $nbar);
        }
        def pq as float init $pbar * (1.0 - $pbar) - ($rf - 1.0) / $rf * $s2;
        if ($r > 1) {
            $out.a = $out.a + $nbar / $nc * ($s2 - ($pq - $hbar / 4.0) / ($nbar - 1.0));
        }
        $out.b = $out.b +
            $nbar / ($nbar - 1.0) * ($pq - (2.0 * $nbar - 1.0) / (4.0 * $nbar) * $hbar);
        $out.c = $out.c + $hbar / 2.0;
    }
    return $out;
}

/**
 * `fStatistics` for each named locus of a set of profiles grouped by
 * population. At each locus a population's sample is its profiles typed
 * there; untyped profiles are left out.
 * @param populations {list of list of forensics.Profile} one profile list per
 *   population
 * @param loci {list of string} the loci to analyse, in output order
 * @return {list of FStatistics} one entry per locus
 * @throws {Error} kind "popstats" on no loci, or a population with nobody
 *   typed at a locus, or as `fStatistics`
 */
export func fStatisticsTable(
    populations as list of list of forensics.Profile,
    loci as list of string) {
    if (len($loci) == 0) {
        fail("fStatisticsTable needs at least one locus");
    }
    def out as list of FStatistics init [];
    for (def locus in $loci) {
        def samples as list of list of forensics.Genotype init [];
        for (def i in 0..len($populations)) {
            def sample as list of forensics.Genotype init typedAt($populations[$i], $locus);
            if (len($sample) == 0) {
                fail("locus '" + $locus + "': population {$i} has nobody typed there");
            }
            $samples[] = $sample;
        }
        $out[] = fStatistics($samples);
    }
    return $out;
}

# typedAt lists the genotypes at `locus` of the profiles typed there.
func typedAt(profiles as list of forensics.Profile, locus as string) {
    def out as list of forensics.Genotype init [];
    for (def j in 0..len($profiles)) {
        if (maps.has($profiles[$j].loci, $locus)) {
            $out[] = $profiles[$j].loci[$locus];
        }
    }
    return $out;
}

/**
 * The multilocus components: `a`, `b`, `c` summed over loci, so the ratios
 * of the result are Weir & Cockerham's (1984) multilocus estimates
 * `Σ a / Σ (a + b + c)` etc. (not an average of per-locus ratios).
 * @param table {list of FStatistics} per-locus components from the same
 *   populations, e.g. from `fStatisticsTable`
 * @return {FStatistics} the sums, with `locus` "" and `loci` the loci summed
 * @throws {Error} kind "popstats" on an empty table, or entries with
 *   different population counts
 */
export func combineFStatistics(table as list of FStatistics) {
    if (len($table) == 0) {
        fail("combineFStatistics needs at least one locus");
    }
    def out as FStatistics init FStatistics{
        locus: "",
        loci: 0,
        populations: $table[0].populations,
        a: 0.0,
        b: 0.0,
        c: 0.0
    };
    for (def i in 0..len($table)) {
        if ($table[$i].populations != $out.populations) {
            fail("cannot combine loci analysed over {$out.populations} and " +
                "{$table[$i].populations} populations");
        }
        $out.loci = $out.loci + $table[$i].loci;
        $out.a = $out.a + $table[$i].a;
        $out.b = $out.b + $table[$i].b;
        $out.c = $out.c + $table[$i].c;
    }
    return $out;
}

# ratio divides, refusing a denominator within DEGENERATE of zero.
func ratio(num as float, den as float, what as string, locus as string) {
    if (math.abs($den) < DEGENERATE) {
        def where as string init "the combined loci";
        if (len($locus) > 0) {
            $where = "locus '" + $locus + "'";
        }
        fail($what + " is undefined at " + $where + ": its denominator is zero " +
            "(no variation to partition)");
    }
    return $num / $den;
}

func needPopulations(s as FStatistics, what as string) {
    if ($s.populations < 2) {
        fail($what + " needs at least two populations; one sample has only fis");
    }
}

/**
 * Weir & Cockerham's (1984) θ, the estimate of F_ST: `a / (a + b + c)`.
 *
 * The coancestry of two alleles within a population relative to the whole;
 * the quantity NRC II's `forensics.THETA_*` values stand in for. Unbiased in
 * the ratio's parts, so small samples from one population give values near 0
 * of either sign.
 * @param stats {FStatistics} components from `fStatistics` or
 *   `combineFStatistics`
 * @return {float} θ; negative values are possible
 * @throws {Error} kind "popstats" for one population, or a zero denominator
 */
export func fst(stats as FStatistics) {
    needPopulations($stats, "fst");
    return ratio($stats.a, $stats.a + $stats.b + $stats.c, "fst", $stats.locus);
}

/**
 * Weir & Cockerham's (1984) F, the estimate of F_IT: `(a + b) / (a + b + c)`.
 * @param stats {FStatistics} components from `fStatistics` or
 *   `combineFStatistics`
 * @return {float} F; negative values are possible
 * @throws {Error} kind "popstats" for one population, or a zero denominator
 */
export func fit(stats as FStatistics) {
    needPopulations($stats, "fit");
    return ratio($stats.a + $stats.b, $stats.a + $stats.b + $stats.c, "fit", $stats.locus);
}

/**
 * Weir & Cockerham's (1984) f, the estimate of F_IS: `b / (b + c)`.
 *
 * Defined for one population too, where it is the within-sample inbreeding
 * coefficient.
 * @param stats {FStatistics} components from `fStatistics` or
 *   `combineFStatistics`
 * @return {float} f; negative values (heterozygote excess) are possible
 * @throws {Error} kind "popstats" on a zero denominator (every population
 *   fixed for one allele)
 */
export func fis(stats as FStatistics) {
    return ratio($stats.b, $stats.b + $stats.c, "fis", $stats.locus);
}

/**
 * The inbreeding coefficient f of one sample, Weir & Cockerham's (1984)
 * estimator with r = 1: `fis(fStatistics([genotypes]))`.
 *
 * In terms of the observed heterozygosity `Ho` and `He = 1 - Σ p²` of the
 * counted frequencies, `1 - 2 (n - 1) Ho / (2 n He - Ho)`.
 * @param genotypes {list of forensics.Genotype} the sample, all at one locus
 * @return {float} f; 1 with no heterozygote, negative for a heterozygote excess
 * @throws {Error} kind "popstats" as `observedHeterozygosity`, on a sample of
 *   one, or a monomorphic sample
 */
export func inbreedingCoefficient(genotypes as list of forensics.Genotype) {
    return fis(fStatistics([$genotypes]));
}

# --- linkage disequilibrium -------------------------------------------------

/**
 * One pair of loci from `ldPermutationTests`.
 * @field first {string} the first locus, in name order
 * @field second {string} the second locus
 * @field n {int} individuals typed at both
 * @field pValue {float} the permutation p-value; 1.0 when no table can be
 *   formed (fewer than two individuals, or one genotype at either locus)
 */
export def struct LdTest {
    first as string,
    second as string,
    n as int,
    pValue as float
};

# LdRun is a permutation run's exceedance count and the generator state it
# ended in, so a series of tests can continue one stream.
def struct LdRun {
    hits as int,
    state as int
};

# Paired is two loci's genotypes of the individuals typed at both, aligned.
def struct Paired {
    first as list of forensics.Genotype,
    second as list of forensics.Genotype
};

# typedAtBoth collects the genotypes at `a` and `b` of every profile typed at
# both, in profile order.
func typedAtBoth(profiles as list of forensics.Profile, a as string, b as string) {
    def out as Paired init Paired{first: [], second: []};
    for (def i in 0..len($profiles)) {
        if (maps.has($profiles[$i].loci, $a) and maps.has($profiles[$i].loci, $b)) {
            $out.first[] = $profiles[$i].loci[$a];
            $out.second[] = $profiles[$i].loci[$b];
        }
    }
    return $out;
}

# checkRun applies hwePermutationTest's limits on permutations and seed.
func checkRun(permutations as int, seed as int) {
    if ($permutations < 1 or $permutations > MAX_PERMUTATIONS) {
        fail("permutations must be between 1 and 1000000, got {$permutations}");
    }
    if ($seed < 1 or $seed >= RAND_MODULUS) {
        fail("seed must be between 1 and 2147483646, got {$seed}");
    }
}

# genotypeClasses numbers each individual's genotype by first appearance.
func genotypeClasses(s as Sample) {
    def index as map of string to int init {};
    def out as list of int init [];
    for (def i in 0..len($s.first)) {
        def key as string init pairKey($s.first[$i], $s.second[$i]);
        if (not maps.has($index, $key)) {
            $index[$key] = len($index);
        }
        $out[] = $index[$key];
    }
    return $out;
}

# classCount is the number of classes genotypeClasses numbered.
func classCount(classes as list of int) {
    def out as int init 0;
    for (def x in $classes) {
        $out = math.max($out, $x + 1);
    }
    return $out;
}

# ldRun draws `permutations` re-pairings of column genotypes with the fixed
# row genotypes and counts those whose table is no more probable than the
# observed one. Given both margins a table's probability is
#   Π r_i! Π c_j! / (n! Π n_ij!)
# so it varies only through -Σ ln n_ij!, scored as hwePermutationTest does.
func ldRun(
    rows as list of int,
    cols as list of int,
    ncols as int,
    nrows as int,
    permutations as int,
    seed as int) {
    def n as int init len($rows);
    def cells as list of int init [];
    def stamp as list of int init [];
    for (def i in 0..($nrows * $ncols)) {
        $cells[] = 0;
        $stamp[] = 0;
    }
    def lnCount as list of float init lnTable($n);
    # observed table, stamped as draw 0's predecessor (-1)
    def observed as float init 0.0;
    for (def i in 0..$n) {
        def cell as int init $rows[$i] * $ncols + $cols[$i];
        if ($stamp[$cell] != -1) {
            $stamp[$cell] = -1;
            $cells[$cell] = 0;
        }
        $cells[$cell] = $cells[$cell] + 1;
        $observed = $observed - $lnCount[$cells[$cell]];
    }
    $observed = $observed + math.log1p(TIE_TOLERANCE);
    def accept as list of int init acceptLimits($n);
    def pool as list of int init $cols;
    def state as int init $seed;
    def hits as int init 0;
    def draw as int init 0;
    # Fisher-Yates from the top fixes position j per pass; its pair with row j
    # is scored in the same pass, and position 0, fixed last, after the loop.
    while ($draw < $permutations) {
        $draw = $draw + 1;
        def score as float init 0.0;
        def j as int init $n - 1;
        while ($j > 0) {
            $state = ($state * RAND_MULTIPLIER) % RAND_MODULUS;
            if ($state - 1 >= $accept[$j + 1]) {
                $state = redraw($state, $accept[$j + 1]);
            }
            def pick as int init ($state - 1) % ($j + 1);
            def x as int init $pool[$pick];
            $pool[$pick] = $pool[$j];
            $pool[$j] = $x;
            def cell as int init $rows[$j] * $ncols + $x;
            if ($stamp[$cell] != $draw) {
                $stamp[$cell] = $draw;
                $cells[$cell] = 0;
            }
            $cells[$cell] = $cells[$cell] + 1;
            $score = $score - $lnCount[$cells[$cell]];
            $j = $j - 1;
        }
        def last as int init $rows[0] * $ncols + $pool[0];
        if ($stamp[$last] != $draw) {
            $stamp[$last] = $draw;
            $cells[$last] = 0;
        }
        $score = $score - $lnCount[$cells[$last] + 1];
        if ($score <= $observed) {
            $hits = $hits + 1;
        }
    }
    return LdRun{hits: $hits, state: $state};
}

# ldPair runs the test on two aligned samples (one entry per individual) from
# `seed`, returning the hits and the generator state for a following test.
func ldPair(
    first as list of forensics.Genotype,
    second as list of forensics.Genotype,
    permutations as int,
    seed as int) {
    def na as int init len($first);
    def nb as int init len($second);
    if ($na != $nb) {
        fail("the two loci's samples differ in length ({$na} and {$nb}); they must " +
            "list the same individuals in the same order");
    }
    def sa as Sample init encode($first);
    def sb as Sample init encode($second);
    if ($sa.locus == $sb.locus) {
        fail("both samples are at locus '" + $sa.locus + "'; the test needs two loci");
    }
    def rows as list of int init genotypeClasses($sa);
    def cols as list of int init genotypeClasses($sb);
    def nrows as int init classCount($rows);
    def ncols as int init classCount($cols);
    if ($nrows == 1 or $ncols == 1) {
        return LdRun{hits: $permutations, state: $seed};
    }
    return ldRun($rows, $cols, $ncols, $nrows, $permutations, $seed);
}

/**
 * The genotypic linkage-disequilibrium test of two loci in one sample, by
 * permutation: a test of association between the two loci's diploid
 * genotypes that assumes no gametic phase and no Hardy-Weinberg proportions
 * (the genotypic-disequilibrium test of Genepop, Raymond & Rousset 1995,
 * there by Markov chain).
 *
 * The individuals' genotypes form a contingency table, rows the genotypes at
 * `first`, columns those at `second`. Under independence, given both
 * margins, the table's probability is `Π r_i! Π c_j! / (n! Π n_ij!)`;
 * shuffling the `second` genotypes among individuals draws tables from
 * exactly that distribution. Zaykin, Zhivotovsky & Weir (1995) permute
 * alleles, testing within- and between-locus association together; whole
 * single-locus genotypes are permuted here, so only the association between
 * loci is tested. The p-value is `(r + 1) / (permutations + 1)`, `r`
 * the shuffled tables no more probable than the observed one (Phipson & Smyth
 * 2010), from the same seeded Park-Miller generator as `hwePermutationTest`.
 *
 * Returns 1.0 when either locus shows one genotype.
 * @param first {list of forensics.Genotype} genotypes at one locus, one per
 *   individual
 * @param second {list of forensics.Genotype} the same individuals, in the same
 *   order, at another locus
 * @param permutations {int} shuffles to draw, 1 to `MAX_PERMUTATIONS`
 * @param seed {int} the generator seed, 1 to 2147483646
 * @return {float} the permutation p-value, in (0, 1]
 * @throws {Error} kind "popstats" as `observedHeterozygosity` for either
 *   sample, on samples of different lengths or at the same locus, or on a
 *   permutation count or seed out of range
 */
export func ldPermutationTest(
    first as list of forensics.Genotype,
    second as list of forensics.Genotype,
    permutations as int,
    seed as int) {
    checkRun($permutations, $seed);
    def run as LdRun init ldPair($first, $second, $permutations, $seed);
    return convert.toFloat($run.hits + 1) / convert.toFloat($permutations + 1);
}

/**
 * `ldPermutationTest` for every pair of loci in a set of profiles (one
 * population), each pair over the individuals typed at both.
 *
 * Loci are taken in name order and pairs as (i, j) with i before j. The
 * tests share one generator stream: the first starts from `seed`, each next
 * one where the previous stopped. With k loci there are k (k - 1) / 2 tests,
 * so about 5% fall below 0.05 under independence; judge them against a
 * multiple-testing correction (Bonferroni: 0.05 / count), not one by one.
 * Cost is about pairs × permutations × individuals.
 * @param profiles {list of forensics.Profile} the sample, one profile per
 *   individual
 * @param permutations {int} shuffles per pair, 1 to `MAX_PERMUTATIONS`
 * @param seed {int} the generator seed, 1 to 2147483646
 * @return {list of LdTest} one entry per locus pair
 * @throws {Error} kind "popstats" on no profiles, fewer than two loci, an
 *   empty allele label, or a permutation count or seed out of range
 */
export func ldPermutationTests(
    profiles as list of forensics.Profile,
    permutations as int,
    seed as int) {
    checkRun($permutations, $seed);
    if (len($profiles) == 0) {
        fail("ldPermutationTests needs at least one profile");
    }
    def seen as map of string to bool init {};
    for (def i in 0..len($profiles)) {
        for (def locus in $profiles[$i].loci) {
            $seen[$locus] = true;
        }
    }
    def names as list of string init lists.sort(maps.keys($seen));
    if (len($names) < 2) {
        fail("ldPermutationTests needs at least two loci");
    }
    def out as list of LdTest init [];
    def state as int init $seed;
    for (def x in 0..len($names)) {
        for (def y in ($x + 1)..len($names)) {
            def pair as Paired init typedAtBoth($profiles, $names[$x], $names[$y]);
            def n as int init len($pair.first);
            def p as float init 1.0;
            if ($n >= 2) {
                def run as LdRun init ldPair($pair.first, $pair.second, $permutations, $state);
                $state = $run.state;
                $p = convert.toFloat($run.hits + 1) / convert.toFloat($permutations + 1);
            }
            $out[] = LdTest{first: $names[$x], second: $names[$y], n: $n, pValue: $p};
        }
    }
    return $out;
}
