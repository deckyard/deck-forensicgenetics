# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * X-chromosomal STRs: hemizygous male and diploid female profiles, match
 * probabilities, pairwise kinship LRs from X-chromosomal IBD coefficients,
 * and paternity indices for a daughter.
 *
 * A male carries one X allele per locus, received from his mother; a father
 * passes his whole X to every daughter and none to a son. A pair's
 * X-chromosomal relationship is therefore its X-IBD coefficients together
 * with both persons' sexes: `(k0, k1, k2)` for two females, `(k0, k1)` for a
 * pair including a male (Pinto, Gusmão & Amorim 2011). Relationships that
 * autosomal STRs cannot separate often differ here: paternal half sisters
 * share an X allele IBD at every locus, maternal half sisters at half of them.
 *
 * Loci are scored one of two ways. The independent model multiplies
 * single-locus probabilities from a core `FrequencyDb`. Tightly linked loci
 * (the four three-locus clusters of the Investigator Argus X-12 kit) are
 * neither independent in the population nor in inheritance; a `LinkageGroup`
 * supplies haplotype frequencies for a cluster, and the cluster is scored as
 * one marker whose alleles are haplotypes, inherited without recombination.
 * Groups are multiplied as if unlinked. Female phase is unknown and summed
 * over. Everything uses theta 0.
 * @module xstr
 * @see https://doi.org/10.1016/j.fsigen.2010.01.011 - Pinto, Gusmão & Amorim
 *   (2011), the X-chromosomal IBD approach
 * @see https://doi.org/10.3389/fgene.2020.00926 - Gomes et al. (2020), Table 5,
 *   X-chromosomal coancestry by sex for common relationships
 * @example
 * import "@deckyard/forensicgenetics/" as forensicgenetics;
 * import "@deckyard/forensicgenetics/xstr.j" as xstr;
 *
 * def model as xstr.XFrequencies init xstr.independentModel($db);
 * def af as xstr.XProfile init xstr.maleProfile("AF", {"DXS7132": "13"});
 * def girl as xstr.XProfile init xstr.femaleProfile("C", [
 *     forensicgenetics.genotype("DXS7132", "13", "14")
 * ]);
 * def lr as float init xstr.kinshipLR($af, $girl, xstr.FATHER_DAUGHTER, $model);
 */
use math;
use maps;
use strings;

import "./forensicgenetics.j" as forensics;
import "./kinship.j" as kinship;

/**
 * The sex label of a male: one X allele per locus.
 */
export def const MALE as string init "male";

/**
 * The sex label of a female: two X alleles per locus.
 */
export def const FEMALE as string init "female";

/**
 * The separator `haplotypeKey` joins a haplotype's alleles with.
 */
export def const HAPLOTYPE_SEPARATOR as string init "-";

/**
 * An X-STR profile and the sex that decides how to read it.
 *
 * A male's genotype at a locus has `a == b`: the one allele he carries, not a
 * homozygote. Build with `maleProfile` or `femaleProfile`.
 * @field id {string} a caller-chosen label
 * @field sex {string} `MALE` or `FEMALE`
 * @field loci {map of string to Genotype} locus name -> the genotype seen there
 */
export def struct XProfile {
    id as string,
    sex as string,
    loci as map of string to forensics.Genotype
};

/**
 * A relationship as X-chromosomal IBD coefficients for an ordered pair of
 * sexes: person A of `sexA`, person B of `sexB`.
 *
 * `k0 + k1 + k2` is 1; `k2` is 0 unless both are female.
 * @field name {string} a label, read A to B ("father-daughter")
 * @field sexA {string} the first person's sex, `MALE` or `FEMALE`
 * @field sexB {string} the second person's sex
 * @field k0 {float} the probability that no X allele is shared IBD
 * @field k1 {float} the probability that exactly one is
 * @field k2 {float} the probability that two are (two females only)
 */
export def struct XRelationship {
    name as string,
    sexA as string,
    sexB as string,
    k0 as float,
    k1 as float,
    k2 as float
};

/**
 * Haplotype frequencies for a cluster of tightly linked X-STR loci.
 *
 * A haplotype is the alleles at `loci`, in that order, joined by
 * `HAPLOTYPE_SEPARATOR` (`haplotypeKey`). A haplotype not listed, or listed
 * below `minFrequency`, is scored at `minFrequency`; with `minFrequency` 0 an
 * unlisted or zero-frequency haplotype throws.
 * @field name {string} the group's label, e.g. "LG1"
 * @field loci {list of string} the linked loci, in key order
 * @field minFrequency {float} the floor for rare and unseen haplotypes, in [0, 1)
 * @field haplotypes {map of string to float} haplotype key -> frequency
 */
export def struct LinkageGroup {
    name as string,
    loci as list of string,
    minFrequency as float,
    haplotypes as map of string to float
};

/**
 * The frequency model: single-locus frequencies for loci in no group, and
 * haplotype frequencies for each linkage group.
 * @field db {FrequencyDb} single-locus X-STR allele frequencies
 * @field groups {list of LinkageGroup} the linkage groups, disjoint
 */
export def struct XFrequencies {
    db as forensics.FrequencyDb,
    groups as list of LinkageGroup
};

/**
 * The most heterozygous loci a female may have within one linkage group. Her
 * phase is summed over `2^(h-1)` haplotype pairs for `h` heterozygous loci,
 * so cost doubles per locus; at the cap one group takes about 0.15 s in a
 * pairwise LR and 0.2 s in a trio index. The Argus X-12 groups have three
 * loci. Exceeding it throws kind "xstr".
 */
export def const MAX_GROUP_HETEROZYGOUS_LOCI as int init 12;

# Phase is one resolution of a person at a unit: two haplotype keys for a
# female, the same key twice for a male.
def struct Phase {
    a as string,
    b as string
};

# Unit is what is scored as one marker: a single locus (group -1) or a whole
# linkage group (its index in XFrequencies.groups).
def struct Unit {
    name as string,
    loci as list of string,
    group as int
};

# --- relationships ----------------------------------------------------------

/**
 * Father and daughter: the father's X is the daughter's paternal X, (0, 1, 0).
 */
export def const FATHER_DAUGHTER as XRelationship init XRelationship{
    name: "father-daughter",
    sexA: "male",
    sexB: "female",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Father and son: a son's X is his mother's, so nothing is shared, (1, 0, 0).
 */
export def const FATHER_SON as XRelationship init XRelationship{
    name: "father-son",
    sexA: "male",
    sexB: "male",
    k0: 1.0,
    k1: 0.0,
    k2: 0.0
};

/**
 * Mother and son, (0, 1, 0).
 */
export def const MOTHER_SON as XRelationship init XRelationship{
    name: "mother-son",
    sexA: "female",
    sexB: "male",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Mother and daughter, (0, 1, 0).
 */
export def const MOTHER_DAUGHTER as XRelationship init XRelationship{
    name: "mother-daughter",
    sexA: "female",
    sexB: "female",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Full sisters: the paternal X always shared, the maternal half the time,
 * (0, 1/2, 1/2).
 */
export def const FULL_SISTERS as XRelationship init XRelationship{
    name: "full sisters",
    sexA: "female",
    sexB: "female",
    k0: 0.0,
    k1: 0.5,
    k2: 0.5
};

/**
 * Full brothers, (1/2, 1/2, 0); the same as maternal half brothers.
 */
export def const FULL_BROTHERS as XRelationship init XRelationship{
    name: "full brothers",
    sexA: "male",
    sexB: "male",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Brother and full sister, (1/2, 1/2, 0).
 */
export def const BROTHER_SISTER as XRelationship init XRelationship{
    name: "brother-sister",
    sexA: "male",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Half sisters through the father, (0, 1, 0): as mother and daughter.
 */
export def const PATERNAL_HALF_SISTERS as XRelationship init XRelationship{
    name: "paternal half sisters",
    sexA: "female",
    sexB: "female",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Half sisters through the mother, (1/2, 1/2, 0).
 */
export def const MATERNAL_HALF_SISTERS as XRelationship init XRelationship{
    name: "maternal half sisters",
    sexA: "female",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Half brothers through the mother, (1/2, 1/2, 0). Paternal half brothers
 * share no X (`FATHER_SON` coefficients).
 */
export def const MATERNAL_HALF_BROTHERS as XRelationship init XRelationship{
    name: "maternal half brothers",
    sexA: "male",
    sexB: "male",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Paternal grandmother and granddaughter, (0, 1, 0): the father's X is one of
 * his mother's.
 */
export def const PATERNAL_GRANDMOTHER_GRANDDAUGHTER as XRelationship init XRelationship{
    name: "paternal grandmother-granddaughter",
    sexA: "female",
    sexB: "female",
    k0: 0.0,
    k1: 1.0,
    k2: 0.0
};

/**
 * Maternal grandmother and granddaughter, (1/2, 1/2, 0).
 */
export def const MATERNAL_GRANDMOTHER_GRANDDAUGHTER as XRelationship init XRelationship{
    name: "maternal grandmother-granddaughter",
    sexA: "female",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Maternal grandfather and granddaughter, (1/2, 1/2, 0).
 */
export def const MATERNAL_GRANDFATHER_GRANDDAUGHTER as XRelationship init XRelationship{
    name: "maternal grandfather-granddaughter",
    sexA: "male",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Paternal aunt (the father's sister) and niece, (1/2, 1/2, 0).
 */
export def const PATERNAL_AUNT_NIECE as XRelationship init XRelationship{
    name: "paternal aunt-niece",
    sexA: "female",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Paternal uncle (the father's brother) and niece, (1/2, 1/2, 0).
 */
export def const PATERNAL_UNCLE_NIECE as XRelationship init XRelationship{
    name: "paternal uncle-niece",
    sexA: "male",
    sexB: "female",
    k0: 0.5,
    k1: 0.5,
    k2: 0.0
};

/**
 * Maternal aunt (the mother's full sister) and niece, (1/4, 3/4, 0).
 */
export def const MATERNAL_AUNT_NIECE as XRelationship init XRelationship{
    name: "maternal aunt-niece",
    sexA: "female",
    sexB: "female",
    k0: 0.25,
    k1: 0.75,
    k2: 0.0
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "xstr", message: $msg, file: "", line: 0, col: 0};
}

# FLOAT_MAX is the largest finite float; the across-unit products check
# against it before multiplying (the interpreter errors on overflow).
def const FLOAT_MAX as float init 1.7976931348623157e308;

# timesLR multiplies one unit's factor into a running product, throwing
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

func checkSex(sex as string, who as string) {
    if ($sex != MALE and $sex != FEMALE) {
        fail($who + ": sex must be \"male\" or \"female\", not '" + $sex + "'");
    }
}

/**
 * A relationship from its X-chromosomal IBD coefficients.
 * @param name {string} a label, read A to B
 * @param sexA {string} the first person's sex, `MALE` or `FEMALE`
 * @param sexB {string} the second person's sex
 * @param k0 {float} the probability that no X allele is shared IBD
 * @param k1 {float} the probability that exactly one is
 * @param k2 {float} the probability that two are; 0 unless both are female
 * @return {XRelationship} the relationship
 * @throws {Error} kind "xstr" on an unknown sex, a coefficient outside [0, 1],
 *   coefficients not summing to 1, or `k2 > 0` for a pair with a male
 */
export func xRelationship(
    name as string,
    sexA as string,
    sexB as string,
    k0 as float,
    k1 as float,
    k2 as float) {
    def r as XRelationship init XRelationship{
        name: $name,
        sexA: $sexA,
        sexB: $sexB,
        k0: $k0,
        k1: $k1,
        k2: $k2
    };
    checkRelationship($r);
    return $r;
}

# checkRelationship is xRelationship's validation, also applied to a struct
# literal passed to the kinship functions.
func checkRelationship(r as XRelationship) {
    def what as string init "relationship '" + $r.name + "'";
    checkSex($r.sexA, $what);
    checkSex($r.sexB, $what);
    for (def k in [$r.k0, $r.k1, $r.k2]) {
        if ($k < 0.0 or $k > 1.0) {
            fail($what + ": IBD coefficients must be in [0, 1]");
        }
    }
    if (math.abs($r.k0 + $r.k1 + $r.k2 - 1.0) > 0.000000001) {
        fail($what + ": IBD coefficients must sum to 1");
    }
    if ($r.k2 > 0.0 and ($r.sexA == MALE or $r.sexB == MALE)) {
        fail($what + ": a male has one X allele, so k2 must be 0");
    }
}

/**
 * The unrelated alternative for a pair of sexes, (1, 0, 0).
 * @param sexA {string} the first person's sex
 * @param sexB {string} the second person's sex
 * @return {XRelationship} the relationship "unrelated"
 * @throws {Error} kind "xstr" on an unknown sex
 */
export func unrelated(sexA as string, sexB as string) {
    return xRelationship("unrelated", $sexA, $sexB, 1.0, 0.0, 0.0);
}

/**
 * The X-chromosomal coancestry (kinship) coefficient: the probability that an
 * X allele drawn at random from each person is IBD. A male's draw is his one
 * allele, so it is `k2/2 + k1/4` for two females, `k1` for two males and
 * `k1/2` for a mixed pair.
 * @param rel {XRelationship} the relationship
 * @return {float} the coancestry coefficient
 * @throws {Error} kind "xstr" on an unknown sex
 */
export func coancestry(rel as XRelationship) {
    checkSex($rel.sexA, "relationship '" + $rel.name + "'");
    checkSex($rel.sexB, "relationship '" + $rel.name + "'");
    if ($rel.sexA == MALE and $rel.sexB == MALE) {
        return $rel.k1;
    }
    if ($rel.sexA == MALE or $rel.sexB == MALE) {
        return $rel.k1 / 2.0;
    }
    return $rel.k2 / 2.0 + $rel.k1 / 4.0;
}

# --- profiles ---------------------------------------------------------------

/**
 * A male X-STR profile: one allele per locus.
 * @param id {string} a label for the profile
 * @param alleles {map of string to string} locus name -> the allele
 * @return {XProfile} the profile; each genotype has `a == b`
 * @throws {Error} kind "xstr" on an empty locus name or allele label
 */
export func maleProfile(id as string, alleles as map of string to string) {
    def out as XProfile init XProfile{id: $id, sex: MALE, loci: {}};
    for (def locus in $alleles) {
        def allele as string init $alleles[$locus];
        if (len($locus) == 0 or len($allele) == 0) {
            fail("profile '" + $id + "': locus names and alleles must not be empty");
        }
        $out.loci[$locus] = forensics.Genotype{locus: $locus, a: $allele, b: $allele};
    }
    return $out;
}

/**
 * A female X-STR profile: one genotype per locus, as for an autosome.
 * @param id {string} a label for the profile
 * @param genotypes {list of Genotype} one genotype per locus
 * @return {XProfile} the profile
 * @throws {Error} kind "xstr" if a locus appears twice
 */
export func femaleProfile(id as string, genotypes as list of forensics.Genotype) {
    def out as XProfile init XProfile{id: $id, sex: FEMALE, loci: {}};
    for (def g in $genotypes) {
        if (maps.has($out.loci, $g.locus)) {
            fail("profile '" + $id + "' lists locus '" + $g.locus + "' twice");
        }
        $out.loci[$g.locus] = $g;
    }
    return $out;
}

# checkProfile rejects an unknown sex and a male genotype with two alleles,
# which a struct literal can produce and which would otherwise be read as a
# female's.
func checkProfile(p as XProfile) {
    checkSex($p.sex, "profile '" + $p.id + "'");
    if ($p.sex == MALE) {
        for (def locus in $p.loci) {
            if ($p.loci[$locus].a != $p.loci[$locus].b) {
                fail("profile '" + $p.id + "' is male but has two alleles at '" + $locus + "'");
            }
        }
    }
}

# --- linkage groups ---------------------------------------------------------

/**
 * The key of a haplotype: its alleles joined by `HAPLOTYPE_SEPARATOR`.
 * @param alleles {list of string} the alleles, in the group's locus order
 * @return {string} the key, e.g. "17-24-12"
 * @throws {Error} kind "xstr" on no alleles, or an allele that is empty or
 *   contains the separator
 */
export func haplotypeKey(alleles as list of string) {
    if (len($alleles) == 0) {
        fail("a haplotype needs at least one allele");
    }
    for (def a in $alleles) {
        if (len($a) == 0 or strings.contains($a, HAPLOTYPE_SEPARATOR)) {
            fail("haplotype allele '" + $a + "' is empty or contains '" + HAPLOTYPE_SEPARATOR +
                "'");
        }
    }
    return strings.join($alleles, HAPLOTYPE_SEPARATOR);
}

/**
 * A linkage group with its haplotype frequency table.
 *
 * Frequencies are as published (male haplotype counts over the number of
 * males). The floor for unseen haplotypes is the caller's choice; none is
 * imposed.
 * @param name {string} the group's label
 * @param loci {list of string} the linked loci, in key order; distinct
 * @param minFrequency {float} the floor for rare and unseen haplotypes, in [0, 1)
 * @param haplotypes {map of string to float} `haplotypeKey` -> frequency
 * @return {LinkageGroup} the group
 * @throws {Error} kind "xstr" on no loci, a repeated or empty locus, a floor
 *   outside [0, 1), a key with the wrong number of alleles or an empty allele,
 *   a frequency outside [0, 1], or frequencies summing past
 *   `MAX_FREQUENCY_SUM`
 */
export func linkageGroup(
    name as string,
    loci as list of string,
    minFrequency as float,
    haplotypes as map of string to float) {
    if (len($loci) == 0) {
        fail("linkage group '" + $name + "' has no loci");
    }
    def seen as map of string to bool init {};
    for (def locus in $loci) {
        if (len($locus) == 0 or maps.has($seen, $locus)) {
            fail("linkage group '" + $name + "': loci must be non-empty and distinct");
        }
        $seen[$locus] = true;
    }
    if ($minFrequency < 0.0 or $minFrequency >= 1.0) {
        fail("linkage group '" + $name + "': minFrequency must be in [0, 1)");
    }
    def total as float init 0.0;
    for (def key in $haplotypes) {
        def parts as list of string init strings.split($key, HAPLOTYPE_SEPARATOR);
        if (len($parts) != len($loci)) {
            fail("linkage group '" + $name + "': haplotype '" + $key + "' does not have " +
                "one allele per locus");
        }
        for (def part in $parts) {
            if (len($part) == 0) {
                fail("linkage group '" + $name + "': haplotype '" + $key + "' has an empty allele");
            }
        }
        def f as float init $haplotypes[$key];
        if ($f < 0.0 or $f > 1.0) {
            fail("linkage group '" + $name + "': haplotype '" + $key +
                "' has a frequency outside [0, 1]");
        }
        $total = $total + $f;
    }
    if ($total > forensics.MAX_FREQUENCY_SUM) {
        fail("linkage group '" + $name + "': haplotype frequencies sum past " +
            "MAX_FREQUENCY_SUM - a misread table?");
    }
    return LinkageGroup{
        name: $name,
        loci: $loci,
        minFrequency: $minFrequency,
        haplotypes: $haplotypes
    };
}

/**
 * The frequency of one haplotype: the listed value, floored at
 * `minFrequency`.
 * @param group {LinkageGroup} the group
 * @param alleles {list of string} the alleles, in the group's locus order
 * @return {float} the frequency, in (0, 1]
 * @throws {Error} kind "xstr" on the wrong number of alleles, an allele that
 *   is empty or contains the separator, or an unlisted or zero-frequency
 *   haplotype in a group with `minFrequency` 0
 */
export func haplotypeFrequency(group as LinkageGroup, alleles as list of string) {
    if (len($alleles) != len($group.loci)) {
        fail("linkage group '" + $group.name + "' needs one allele per locus");
    }
    return keyFrequency($group, haplotypeKey($alleles));
}

# keyFrequency is haplotypeFrequency for a key already built.
func keyFrequency(group as LinkageGroup, key as string) {
    def p as float init 0.0;
    if (maps.has($group.haplotypes, $key)) {
        $p = $group.haplotypes[$key];
    }
    if ($p < $group.minFrequency) {
        $p = $group.minFrequency;
    }
    if ($p <= 0.0) {
        fail("linkage group '" + $group.name + "': haplotype '" + $key +
            "' has no frequency and the group has no floor");
    }
    return $p;
}

/**
 * The Investigator Argus X-12 linkage groups, LG1 to LG4, three loci each,
 * as listed in the QIAGEN Investigator Argus X-12 QS Handbook (07/2022),
 * Table 1 (read 2026-11-05). Locus names only; no frequencies are bundled.
 * @return {map of string to list of string} group name -> its loci, in the
 *   handbook's order
 */
export func argusX12Groups() {
    return {
        "LG1": ["DXS8378", "DXS10135", "DXS10148"],
        "LG2": ["DXS7132", "DXS10074", "DXS10079"],
        "LG3": ["HPRTB", "DXS10101", "DXS10103"],
        "LG4": ["DXS7423", "DXS10134", "DXS10146"]
    };
}

/**
 * The independent-marker model: every locus scored on its own from `db`.
 * @param db {FrequencyDb} single-locus X-STR allele frequencies
 * @return {XFrequencies} the model, with no linkage groups
 */
export func independentModel(db as forensics.FrequencyDb) {
    return XFrequencies{db: $db, groups: []};
}

/**
 * A copy of the model with a linkage group added: its loci are then scored
 * together, from its haplotype table, and no longer from `db`.
 * @param xf {XFrequencies} the model
 * @param group {LinkageGroup} the group to add
 * @return {XFrequencies} the updated copy
 * @throws {Error} kind "xstr" if the group's name or a locus is already in a
 *   group of the model
 */
export func withLinkageGroup(xf as XFrequencies, group as LinkageGroup) {
    for (def i in 0..len($xf.groups)) {
        if ($xf.groups[$i].name == $group.name) {
            fail("linkage group '" + $group.name + "' is already in the model");
        }
        for (def locus in $group.loci) {
            if (contains($xf.groups[$i].loci, $locus)) {
                fail("locus '" + $locus + "' is already in linkage group '" +
                    $xf.groups[$i].name + "'");
            }
        }
    }
    $xf.groups[] = $group;
    return $xf;
}

func contains(xs as list of string, x as string) {
    for (def y in $xs) {
        if ($y == $x) {
            return true;
        }
    }
    return false;
}

# --- units: what is scored as one marker ------------------------------------

# groupIndex maps each grouped locus to its group's index.
func groupIndex(xf as XFrequencies) {
    def out as map of string to int init {};
    for (def i in 0..len($xf.groups)) {
        for (def locus in $xf.groups[$i].loci) {
            $out[$locus] = $i;
        }
    }
    return $out;
}

# typedCount is how many of `loci` the profile is typed at.
func typedCount(p as XProfile, loci as list of string) {
    def n as int init 0;
    for (def locus in $loci) {
        if (maps.has($p.loci, $locus)) {
            $n = $n + 1;
        }
    }
    return $n;
}

# checkGroupTyping throws when a person is typed at part of a linkage group:
# the haplotype is incomplete, and dropping or splitting the group silently
# would change the answer.
func checkGroupTyping(xf as XFrequencies, people as list of XProfile) {
    for (def g in 0..len($xf.groups)) {
        def size as int init len($xf.groups[$g].loci);
        for (def i in 0..len($people)) {
            def n as int init typedCount($people[$i], $xf.groups[$g].loci);
            if ($n > 0 and $n < $size) {
                fail("profile '" + $people[$i].id + "' is typed at only part of linkage group '" +
                    $xf.groups[$g].name + "'; type the whole group or score with a model" +
                    " that does not group these loci");
            }
        }
    }
}

# typedByAll is whether every person carries every one of `loci`.
func typedByAll(people as list of XProfile, loci as list of string) {
    for (def i in 0..len($people)) {
        if (typedCount($people[$i], $loci) < len($loci)) {
            return false;
        }
    }
    return true;
}

# unitsFor returns the units every person is typed at, in the first person's
# locus order; a group sits where the first of its loci does. A single locus
# every person carries must be in the database.
func unitsFor(xf as XFrequencies, people as list of XProfile, who as string) {
    checkGroupTyping($xf, $people);
    def gi as map of string to int init groupIndex($xf);
    def done as map of int to bool init {};
    def out as list of Unit init [];
    for (def locus in $people[0].loci) {
        def g as int init -1;
        if (maps.has($gi, $locus)) {
            $g = $gi[$locus];
        }
        if ($g >= 0 and maps.has($done, $g)) {
            continue;
        }
        def loci as list of string init [$locus];
        def name as string init $locus;
        if ($g >= 0) {
            $done[$g] = true;
            $loci = $xf.groups[$g].loci;
            $name = $xf.groups[$g].name;
        }
        if (not typedByAll($people, $loci)) {
            continue;
        }
        if ($g < 0 and not forensics.hasLocus($xf.db, $locus)) {
            fail($who + " are all typed at locus '" + $locus + "', which frequency database '" +
                $xf.db.name + "' does not cover and no linkage group holds; it will not be" +
                " dropped silently");
        }
        $out[] = Unit{name: $name, loci: $loci, group: $g};
    }
    if (len($out) == 0) {
        fail($who + " share no X-STR locus");
    }
    return $out;
}

# frequency is the population frequency of an allele (single locus) or a
# haplotype key (group) at a unit.
func frequency(xf as XFrequencies, u as Unit, key as string) {
    if ($u.group < 0) {
        return forensics.alleleFrequency($xf.db, $u.name, $key);
    }
    return keyFrequency($xf.groups[$u.group], $key);
}

# phases returns a person's resolutions at a unit: one for a male or at a
# single locus; for a female at a group, every unordered pair of haplotypes
# her genotypes allow, 2^(h-1) for h heterozygous loci. The resolutions are
# disjoint events, so probabilities over them add.
func phases(p as XProfile, u as Unit) {
    if ($u.group < 0) {
        def g as forensics.Genotype init $p.loci[$u.name];
        return [Phase{a: $g.a, b: $g.b}];
    }
    def hets as int init 0;
    for (def locus in $u.loci) {
        if ($p.loci[$locus].a != $p.loci[$locus].b) {
            $hets = $hets + 1;
        }
    }
    if ($hets > MAX_GROUP_HETEROZYGOUS_LOCI) {
        fail("profile '" + $p.id + "' is heterozygous at {$hets} loci of linkage group '" +
            $u.name + "'; MAX_GROUP_HETEROZYGOUS_LOCI is {MAX_GROUP_HETEROZYGOUS_LOCI} (the" +
            " phase sum has 2^(h-1) terms)");
    }
    def firsts as list of list of string init [[]];
    def seconds as list of list of string init [[]];
    $hets = 0;
    for (def locus in $u.loci) {
        def a as string init $p.loci[$locus].a;
        def b as string init $p.loci[$locus].b;
        # Validated once here, so the keys below are joined without haplotypeKey.
        haplotypeKey([$a, $b]);
        def n as int init len($firsts);
        if ($a != $b and $hets > 0) {
            for (def i in 0..$n) {
                $firsts[] = $firsts[$i];
                $seconds[] = $seconds[$i];
                $firsts[len($firsts) - 1][] = $b;
                $seconds[len($seconds) - 1][] = $a;
            }
        }
        for (def i in 0..$n) {
            $firsts[$i][] = $a;
            $seconds[$i][] = $b;
        }
        if ($a != $b) {
            $hets = $hets + 1;
        }
    }
    def out as list of Phase init [];
    for (def i in 0..len($firsts)) {
        $out[] = Phase{
            a: strings.join($firsts[$i], HAPLOTYPE_SEPARATOR),
            b: strings.join($seconds[$i], HAPLOTYPE_SEPARATOR)
        };
    }
    return $out;
}

# --- probabilities at one unit ----------------------------------------------

# phaseProbability is P(resolution): p for a male, p^2 or 2 p q for a female.
func phaseProbability(xf as XFrequencies, u as Unit, sex as string, r as Phase) {
    def pa as float init frequency($xf, $u, $r.a);
    if ($sex == MALE) {
        return $pa;
    }
    if ($r.a == $r.b) {
        return $pa * $pa;
    }
    return 2.0 * $pa * frequency($xf, $u, $r.b);
}

# PhaseIndex summarises a female's resolutions at a unit for the pair sums:
# `total` is Σ P(resolution); `other[x]` is the haplotype paired with `x` in
# the one resolution that carries `x` (a genotype fixes the complement of each
# of its haplotypes, so there is at most one), `otherFreq[x]` its frequency.
def struct PhaseIndex {
    total as float,
    other as map of string to string,
    otherFreq as map of string to float
};

# Joint is the joint probability of two observations at a unit split by IBD
# state: k0 c0 + k1 c1 + k2 c2 for any relationship of the pair's sexes.
def struct Joint {
    c0 as float,
    c1 as float,
    c2 as float
};

# phaseIndex builds the PhaseIndex of a female's resolutions in one pass.
func phaseIndex(xf as XFrequencies, u as Unit, rs as list of Phase) {
    def out as PhaseIndex init PhaseIndex{total: 0.0, other: {}, otherFreq: {}};
    for (def i in 0..len($rs)) {
        def a as string init $rs[$i].a;
        def b as string init $rs[$i].b;
        def fa as float init frequency($xf, $u, $a);
        def fb as float init $fa;
        if ($a != $b) {
            $fb = frequency($xf, $u, $b);
            $out.total = $out.total + 2.0 * $fa * $fb;
        } else {
            $out.total = $out.total + $fa * $fa;
        }
        $out.other[$a] = $b;
        $out.otherFreq[$a] = $fb;
        $out.other[$b] = $a;
        $out.otherFreq[$b] = $fa;
    }
    return $out;
}

# heldWith is Σ over a female's resolutions of Q(G, x): P(her observation |
# one of her haplotypes is fixed at `x`, the other drawn from the population).
func heldWith(ix as PhaseIndex, x as string) {
    if (maps.has($ix.otherFreq, $x)) {
        return $ix.otherFreq[$x];
    }
    return 0.0;
}

# jointParts sums the pair terms over both persons' resolutions,
#   male-male      k0 pa pb + k1 [a == b] pa
#   male-female    k0 px P(Gb) + k1 px Q(Gb, x)
#   female-female  k0 P(Ga) P(Gb) + k1 P(Ga) (Q(Gb, a1) + Q(Gb, a2)) / 2
#                  + k2 [Ga == Gb] P(Ga)
# factorised: the k0 part is P(A) P(B), and the k1 / k2 parts look each of A's
# resolutions up in B's PhaseIndex, so the cost is O(|A| + |B|) resolutions,
# not their product.
func jointParts(xf as XFrequencies, u as Unit, a as XProfile, b as XProfile) {
    def ra as list of Phase init phases($a, $u);
    def rb as list of Phase init phases($b, $u);
    if ($a.sex == MALE and $b.sex == MALE) {
        def pa as float init frequency($xf, $u, $ra[0].a);
        def same as float init 0.0;
        if ($ra[0].a == $rb[0].a) {
            $same = $pa;
        }
        return Joint{c0: $pa * frequency($xf, $u, $rb[0].a), c1: $same, c2: 0.0};
    }
    if ($a.sex == MALE) {
        def ib as PhaseIndex init phaseIndex($xf, $u, $rb);
        def px as float init frequency($xf, $u, $ra[0].a);
        return Joint{c0: $px * $ib.total, c1: $px * heldWith($ib, $ra[0].a), c2: 0.0};
    }
    if ($b.sex == MALE) {
        def ia as PhaseIndex init phaseIndex($xf, $u, $ra);
        def py as float init frequency($xf, $u, $rb[0].a);
        return Joint{c0: $py * $ia.total, c1: $py * heldWith($ia, $rb[0].a), c2: 0.0};
    }
    def ib as PhaseIndex init phaseIndex($xf, $u, $rb);
    def out as Joint init Joint{c0: 0.0, c1: 0.0, c2: 0.0};
    def total as float init 0.0;
    for (def i in 0..len($ra)) {
        def a1 as string init $ra[$i].a;
        def a2 as string init $ra[$i].b;
        def pga as float init phaseProbability($xf, $u, FEMALE, $ra[$i]);
        $total = $total + $pga;
        $out.c1 = $out.c1 + $pga * (0.5 * heldWith($ib, $a1) + 0.5 * heldWith($ib, $a2));
        if (maps.has($ib.other, $a1) and $ib.other[$a1] == $a2) {
            $out.c2 = $out.c2 + $pga;
        }
    }
    $out.c0 = $total * $ib.total;
    return $out;
}

# jointUnder is the joint probability of the two observations under `rel`.
func jointUnder(j as Joint, rel as XRelationship) {
    return $rel.k0 * $j.c0 + $rel.k1 * $j.c1 + $rel.k2 * $j.c2;
}

# unitJoint is the joint probability of A's and B's observations at a unit
# under `rel`, summed over both persons' resolutions.
func unitJoint(xf as XFrequencies, u as Unit, a as XProfile, b as XProfile, rel as XRelationship) {
    return jointUnder(jointParts($xf, $u, $a, $b), $rel);
}

# unitProbability is P(person's observation at a unit), summed over phase.
func unitProbability(xf as XFrequencies, u as Unit, p as XProfile) {
    def out as float init 0.0;
    for (def r in phases($p, $u)) {
        $out = $out + phaseProbability($xf, $u, $p.sex, $r);
    }
    return $out;
}

# --- match probability ------------------------------------------------------

# profileUnits checks a profile and returns its units, every one covered.
func profileUnits(xf as XFrequencies, p as XProfile) {
    checkProfile($p);
    if (len($p.loci) == 0) {
        fail("profile '" + $p.id + "' has no locus");
    }
    return unitsFor($xf, [$p], "profile '" + $p.id + "'");
}

/**
 * The probability that an unrelated person of the same sex has the profile,
 * theta 0: a male's allele or haplotype frequency, a female's Hardy-Weinberg
 * genotype frequency (summed over phase within a linkage group), multiplied
 * across loci and groups. On a long panel the product can underflow to 0;
 * use `log10MatchProbability`.
 * @param xf {XFrequencies} the frequency model
 * @param p {XProfile} the profile
 * @return {float} the match probability
 * @throws {Error} kind "xstr" on an empty or malformed profile, a partly typed
 *   linkage group, a female heterozygous at more than
 *   `MAX_GROUP_HETEROZYGOUS_LOCI` loci of a group, or a locus neither the
 *   database nor a group covers
 */
export func matchProbability(xf as XFrequencies, p as XProfile) {
    def out as float init 1.0;
    for (def u in profileUnits($xf, $p)) {
        $out = $out * unitProbability($xf, $u, $p);
    }
    return $out;
}

/**
 * `log10` of `matchProbability`, summed per locus or group so a long panel
 * does not underflow.
 * @param xf {XFrequencies} the frequency model
 * @param p {XProfile} the profile
 * @return {float} log10 of the match probability
 * @throws {Error} kind "xstr" as `matchProbability`
 */
export func log10MatchProbability(xf as XFrequencies, p as XProfile) {
    def out as float init 0.0;
    for (def u in profileUnits($xf, $p)) {
        $out = $out + math.log10(unitProbability($xf, $u, $p));
    }
    return $out;
}

# --- pairwise kinship -------------------------------------------------------

# pairSetup checks two profiles against the relationships and returns the
# units both carry.
func pairSetup(
    xf as XFrequencies,
    a as XProfile,
    b as XProfile,
    rel as XRelationship,
    alt as XRelationship) {
    checkProfile($a);
    checkProfile($b);
    for (def r in [$rel, $alt]) {
        checkRelationship($r);
        if ($r.sexA != $a.sex or $r.sexB != $b.sex) {
            fail("relationship '" + $r.name + "' is for a " + $r.sexA + " and a " + $r.sexB +
                ", but profiles '" + $a.id + "' and '" + $b.id + "' are " + $a.sex + " and " +
                $b.sex);
        }
    }
    return unitsFor($xf, [$a, $b], "profiles '" + $a.id + "' and '" + $b.id + "'");
}

# unitLR is the likelihood ratio at one unit.
func unitLR(
    xf as XFrequencies,
    u as Unit,
    a as XProfile,
    b as XProfile,
    rel as XRelationship,
    alt as XRelationship) {
    def j as Joint init jointParts($xf, $u, $a, $b);
    def under as float init jointUnder($j, $alt);
    if ($under <= 0.0) {
        fail("'" + $u.name + "': the alternative '" + $alt.name +
            "' gives these observations probability zero, so no ratio is defined");
    }
    return jointUnder($j, $rel) / $under;
}

/**
 * The loci and linkage groups a pair is scored over, in the first profile's
 * order: a group's name stands for its loci.
 * @param xf {XFrequencies} the frequency model
 * @param a {XProfile} one profile
 * @param b {XProfile} the other
 * @return {list of string} the unit names
 * @throws {Error} kind "xstr" as `kinshipLR`
 */
export func pairUnits(xf as XFrequencies, a as XProfile, b as XProfile) {
    checkProfile($a);
    checkProfile($b);
    def out as list of string init [];
    for (def u in unitsFor($xf, [$a, $b], "profiles '" + $a.id + "' and '" + $b.id + "'")) {
        $out[] = $u.name;
    }
    return $out;
}

/**
 * The X-STR kinship likelihood ratio of a relationship against unrelated.
 *
 * At each locus or group the joint probability of the two observations
 * under the X-IBD coefficients, over the same with `k0 = 1`; multiplied
 * across loci and groups. A locus only one profile carries is skipped.
 * @param a {XProfile} person A, of `rel.sexA`
 * @param b {XProfile} person B, of `rel.sexB`
 * @param rel {XRelationship} the relationship to test
 * @param xf {XFrequencies} the frequency model
 * @return {float} the likelihood ratio
 * @throws {Error} kind "xstr" if the relationship fails `xRelationship`'s
 *   checks or its sexes do not match the profiles, the profiles share no
 *   locus, share a locus nothing covers, one is typed at part of a linkage
 *   group, a female is heterozygous at more than
 *   `MAX_GROUP_HETEROZYGOUS_LOCI` loci of a group, or the LR exceeds the
 *   float range (use `log10KinshipLRAgainst`)
 */
export func kinshipLR(a as XProfile, b as XProfile, rel as XRelationship, xf as XFrequencies) {
    return kinshipLRAgainst($a, $b, $rel, unrelated($rel.sexA, $rel.sexB), $xf);
}

/**
 * The X-STR kinship likelihood ratio of one relationship against another,
 * for the same pair of sexes.
 * @param a {XProfile} person A
 * @param b {XProfile} person B
 * @param rel {XRelationship} the numerator hypothesis
 * @param alt {XRelationship} the denominator hypothesis
 * @param xf {XFrequencies} the frequency model
 * @return {float} the likelihood ratio
 * @throws {Error} kind "xstr" as `kinshipLR`, or if `alt` gives an
 *   observation probability zero
 */
export func kinshipLRAgainst(
    a as XProfile,
    b as XProfile,
    rel as XRelationship,
    alt as XRelationship,
    xf as XFrequencies) {
    def out as float init 1.0;
    for (def u in pairSetup($xf, $a, $b, $rel, $alt)) {
        $out = timesLR(
            $out,
            unitLR($xf, $u, $a, $b, $rel, $alt),
            "X-STR kinship LR",
            $a.id,
            $b.id,
            "; use log10KinshipLRAgainst");
    }
    return $out;
}

/**
 * `log10` of `kinshipLRAgainst`, summed per locus or group.
 * @param a {XProfile} person A
 * @param b {XProfile} person B
 * @param rel {XRelationship} the numerator hypothesis
 * @param alt {XRelationship} the denominator hypothesis
 * @param xf {XFrequencies} the frequency model
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "xstr" as `kinshipLRAgainst`, or if the ratio is 0 at
 *   some locus (an exclusion has no logarithm)
 */
export func log10KinshipLRAgainst(
    a as XProfile,
    b as XProfile,
    rel as XRelationship,
    alt as XRelationship,
    xf as XFrequencies) {
    def out as float init 0.0;
    for (def u in pairSetup($xf, $a, $b, $rel, $alt)) {
        def lr as float init unitLR($xf, $u, $a, $b, $rel, $alt);
        if ($lr <= 0.0) {
            fail("'" + $u.name + "': the ratio is 0 (an exclusion), which has no logarithm");
        }
        $out = $out + math.log10($lr);
    }
    return $out;
}

# --- paternity of a daughter ------------------------------------------------

# fatherTransmits maps each haplotype in the daughter's resolutions to the
# probability that a man with haplotype `f` passes it on: per-locus stepwise
# mutation, independent across a group's loci, each locus's probabilities
# computed once. A father passes his single X whole.
func fatherTransmits(
    model as kinship.MutationModel,
    u as Unit,
    f as string,
    daughter as XProfile,
    dphases as list of Phase) {
    def fs as list of string init [$f];
    if ($u.group >= 0) {
        $fs = strings.split($f, HAPLOTYPE_SEPARATOR);
    }
    def table as list of map of string to float init [];
    for (def i in 0..len($u.loci)) {
        def locus as string init $u.loci[$i];
        def g as forensics.Genotype init $daughter.loci[$locus];
        def row as map of string to float init {};
        for (def x in [$g.a, $g.b]) {
            $row[$x] = mutationStep($model, $locus, $fs[$i], $x);
        }
        $table[] = $row;
    }
    def out as map of string to float init {};
    for (def r in 0..len($dphases)) {
        for (def x in [$dphases[$r].a, $dphases[$r].b]) {
            if (maps.has($out, $x)) {
                continue;
            }
            def xs as list of string init [$x];
            if ($u.group >= 0) {
                $xs = strings.split($x, HAPLOTYPE_SEPARATOR);
            }
            def p as float init 1.0;
            for (def i in 0..len($xs)) {
                $p = $p * $table[$i][$xs[$i]];
            }
            $out[$x] = $p;
        }
    }
    return $out;
}

# mutationStep is kinship.mutationProbability with its "kinship" errors (a
# non-numeric allele under a non-zero rate) rethrown as this module's.
func mutationStep(
    model as kinship.MutationModel,
    locus as string,
    source as string,
    target as string) {
    try {
        return kinship.mutationProbability($model, $locus, $source, $target);
    } catch (e) {
        if ($e.kind != "kinship") {
            throw $e;
        }
        fail($e.message);
    }
    return 0.0;
}

# motherTransmits maps each haplotype the mother can pass on to the
# probability that she does: Σ over her resolutions of P(resolution | her
# observation) times the haplotype's share of the resolution (1/2 per copy).
func motherTransmits(xf as XFrequencies, u as Unit, mphases as list of Phase) {
    def weights as list of float init [];
    def total as float init 0.0;
    for (def i in 0..len($mphases)) {
        def w as float init phaseProbability($xf, $u, FEMALE, $mphases[$i]);
        $weights[] = $w;
        $total = $total + $w;
    }
    def out as map of string to float init {};
    for (def i in 0..len($mphases)) {
        def half as float init $weights[$i] / $total / 2.0;
        for (def x in [$mphases[$i].a, $mphases[$i].b]) {
            if (maps.has($out, $x)) {
                $out[$x] = $out[$x] + $half;
            } else {
                $out[$x] = $half;
            }
        }
    }
    return $out;
}

# transmitted reads a haplotype's probability from motherTransmits' map.
func transmitted(tm as map of string to float, x as string) {
    if (maps.has($tm, $x)) {
        return $tm[$x];
    }
    return 0.0;
}

# trioUnit is the trio index at one unit: X / Y over the daughter's
# resolutions (c1, c2),
#   X = Σ T_M(c1) F(c2) + T_M(c2) F(c1)     (one term if c1 == c2)
#   Y = Σ T_M(c1) p(c2) + T_M(c2) p(c1)
func trioUnit(
    xf as XFrequencies,
    u as Unit,
    mother as XProfile,
    daughter as XProfile,
    father as XProfile,
    model as kinship.MutationModel) {
    def tm as map of string to float init motherTransmits($xf, $u, phases($mother, $u));
    def dph as list of Phase init phases($daughter, $u);
    def ft as map of string to float init fatherTransmits(
        $model,
        $u,
        phases($father, $u)[0].a,
        $daughter,
        $dph);
    def x as float init 0.0;
    def y as float init 0.0;
    for (def r in $dph) {
        def t1 as float init transmitted($tm, $r.a);
        def t2 as float init transmitted($tm, $r.b);
        $x = $x + $t1 * $ft[$r.b];
        $y = $y + $t1 * frequency($xf, $u, $r.b);
        if ($r.a != $r.b) {
            $x = $x + $t2 * $ft[$r.a];
            $y = $y + $t2 * frequency($xf, $u, $r.a);
        }
    }
    if ($y <= 0.0) {
        fail("'" + $u.name + "': the mother cannot have transmitted either of the daughter's" +
            " alleles - this is a maternal exclusion, not a paternity index");
    }
    return $x / $y;
}

# duoUnit is the duo index at one unit, over the daughter's resolutions:
#   X = Σ F(c1) p(c2) + F(c2) p(c1)     (one term if c1 == c2)
#   Y = Σ P(c1, c2)
func duoUnit(
    xf as XFrequencies,
    u as Unit,
    daughter as XProfile,
    father as XProfile,
    model as kinship.MutationModel) {
    def dph as list of Phase init phases($daughter, $u);
    def ft as map of string to float init fatherTransmits(
        $model,
        $u,
        phases($father, $u)[0].a,
        $daughter,
        $dph);
    def x as float init 0.0;
    def y as float init 0.0;
    for (def r in $dph) {
        $x = $x + $ft[$r.a] * frequency($xf, $u, $r.b);
        if ($r.a != $r.b) {
            $x = $x + $ft[$r.b] * frequency($xf, $u, $r.a);
        }
        $y = $y + phaseProbability($xf, $u, FEMALE, $r);
    }
    return $x / $y;
}

func requireSex(p as XProfile, sex as string, role as string) {
    checkProfile($p);
    if ($p.sex != $sex) {
        fail("the " + $role + " profile '" + $p.id + "' must be " + $sex);
    }
}

/**
 * The combined paternity index for a daughter with her mother typed.
 *
 * The father passes his X whole, so his transmission is his allele
 * (haplotype) with stepwise mutation per locus; the mother passes one of hers
 * (phase averaged within a group). `X / Y` per locus or group, multiplied,
 * over those all three carry.
 * @param mother {XProfile} the mother, female
 * @param daughter {XProfile} the daughter, female
 * @param allegedFather {XProfile} the alleged father, male
 * @param xf {XFrequencies} the frequency model
 * @param model {MutationModel} the alleged father's mutation model
 * @return {float} the combined paternity index
 * @throws {Error} kind "xstr" on a wrong sex, no shared locus, a locus nothing
 *   covers, a partly typed group, a female heterozygous at more than
 *   `MAX_GROUP_HETEROZYGOUS_LOCI` loci of a group, a maternal exclusion, a
 *   non-numeric allele where `model` needs a step count, or an index past the
 *   float range
 */
export func daughterTrioIndex(
    mother as XProfile,
    daughter as XProfile,
    allegedFather as XProfile,
    xf as XFrequencies,
    model as kinship.MutationModel) {
    requireSex($mother, FEMALE, "mother");
    requireSex($daughter, FEMALE, "daughter");
    requireSex($allegedFather, MALE, "alleged father");
    def units as list of Unit init unitsFor(
        $xf,
        [$daughter, $mother, $allegedFather],
        "profiles '" + $mother.id + "', '" + $daughter.id + "' and '" + $allegedFather.id + "'");
    def out as float init 1.0;
    for (def u in $units) {
        $out = timesLR(
            $out,
            trioUnit($xf, $u, $mother, $daughter, $allegedFather, $model),
            "daughter paternity index",
            $daughter.id,
            $allegedFather.id,
            "");
    }
    return $out;
}

/**
 * The combined paternity index for a daughter without the mother: the
 * daughter's paternal allele from the alleged father, her other from the
 * population. With `NO_MUTATION` it equals `kinshipLR` for `FATHER_DAUGHTER`.
 * @param daughter {XProfile} the daughter, female
 * @param allegedFather {XProfile} the alleged father, male
 * @param xf {XFrequencies} the frequency model
 * @param model {MutationModel} the alleged father's mutation model
 * @return {float} the combined paternity index
 * @throws {Error} kind "xstr" on a wrong sex, no shared locus, a locus nothing
 *   covers, a partly typed group, a daughter heterozygous at more than
 *   `MAX_GROUP_HETEROZYGOUS_LOCI` loci of a group, a non-numeric allele where
 *   `model` needs a step count, or an index past the float range
 */
export func daughterDuoIndex(
    daughter as XProfile,
    allegedFather as XProfile,
    xf as XFrequencies,
    model as kinship.MutationModel) {
    requireSex($daughter, FEMALE, "daughter");
    requireSex($allegedFather, MALE, "alleged father");
    def units as list of Unit init unitsFor(
        $xf,
        [$daughter, $allegedFather],
        "profiles '" + $daughter.id + "' and '" + $allegedFather.id + "'");
    def out as float init 1.0;
    for (def u in $units) {
        $out = timesLR(
            $out,
            duoUnit($xf, $u, $daughter, $allegedFather, $model),
            "daughter paternity index",
            $daughter.id,
            $allegedFather.id,
            "");
    }
    return $out;
}
