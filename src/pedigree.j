# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Whole-pedigree likelihoods by Elston-Stewart peeling.
 *
 * For three or more relatives at once, where `kinship.j`'s pairwise IBD form
 * loses information: a missing person identified through several relatives, or
 * a paternity case made through the alleged father's parents.
 *
 * The likelihood at one locus sums over every genotype of every untyped member.
 * Peeling passes vectors over genotypes along the edges of the
 * individual/family tree (the message-passing form), which is linear in the
 * number of people instead of exponential.
 *
 * Simple pedigrees only: `validate` rejects a loop (a consanguineous mating or
 * a repeated marriage), since loops need cutset conditioning, which is not
 * implemented.
 *
 * Alleles no typed member shows are lumped into one class carrying their
 * combined frequency. The lumping is exact, because untyped members' alleles
 * are exchangeable. Theta is 0 throughout, as in `kinship.j`.
 * @module pedigree
 * @example
 * import "@deckyard/forensicgenetics/pedigree.j" as pedigree;
 *
 * # Is the child the son of AF, or of an unrelated, untyped man UN?
 * def base as pedigree.Pedigree init pedigree.withFounder(
 *     pedigree.withFounder(pedigree.withFounder(pedigree.pedigree("trio"), "AF", "m"),
 *     "MO", "f"), "UN", "m");
 * def h1 as pedigree.Pedigree init pedigree.withChild($base, "CH", "m", "AF", "MO");
 * def h2 as pedigree.Pedigree init pedigree.withChild($base, "CH", "m", "UN", "MO");
 * def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
 */
use math;
use maps;
use lists;

import "./forensicgenetics.j" as forensics;

/**
 * One member of a pedigree. A founder has `""` for both parents; everyone else
 * names both.
 * @field id {string} the member's identifier, unique within the pedigree
 * @field sex {string} "m", "f" or ""; not used by the autosomal calculation
 * @field father {string} the father's id, or "" for a founder
 * @field mother {string} the mother's id, or "" for a founder
 */
export def struct Person {
    id as string,
    sex as string,
    father as string,
    mother as string
};

/**
 * A pedigree: a set of members and the parentages joining them.
 * @field name {string} a label for the pedigree, used in error messages
 * @field members {map of string to Person} id -> member
 */
export def struct Pedigree {
    name as string,
    members as map of string to Person
};

# Frame is everything the peeling needs at one locus, built once by
# `locusLikelihood` and threaded through the message recursion.
#
# Genotypes are indexed, not named: `gi[k]` and `gj[k]` are the two allele
# indices of genotype k, with gi <= gj so each unordered pair appears once.
# Allele indices run over the alleles actually observed at this locus, plus -
# when their frequencies leave room - one final lumped class standing for
# every allele nobody showed.
def struct Frame {
    alleleFreq as list of float,
    nAlleles as int,
    pairIdx as list of list of int,
    gi as list of int,
    gj as list of int,
    famFather as list of string,
    famMother as list of string,
    famChildren as list of list of string,
    incidence as map of string to list of int,
    founder as map of string to bool,
    observed as map of string to int
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "pedigree", message: $msg, file: "", line: 0, col: 0};
}

# --- building a pedigree ----------------------------------------------------

/**
 * An empty pedigree.
 * @param name {string} a label for the pedigree
 * @return {Pedigree} a pedigree with no members
 */
export func pedigree(name as string) {
    return Pedigree{name: $name, members: {}};
}

/**
 * A copy of `ped` with a founder added: a member whose parents are outside the
 * pedigree.
 * @param ped {Pedigree} the pedigree to extend
 * @param id {string} the new member's id
 * @param sex {string} "m", "f" or ""
 * @return {Pedigree} the extended copy
 * @throws {Error} kind "pedigree" on an empty or duplicate id, or an unknown sex
 */
export func withFounder(ped as Pedigree, id as string, sex as string) {
    checkNewMember($ped, $id, $sex);
    def out as Pedigree init $ped;
    $out.members[$id] = Person{id: $id, sex: $sex, father: "", mother: ""};
    return $out;
}

/**
 * A copy of `ped` with a child added.
 *
 * Both parents must be named. They need not be members yet, but `validate`
 * and every likelihood function require them to be.
 * @param ped {Pedigree} the pedigree to extend
 * @param id {string} the new member's id
 * @param sex {string} "m", "f" or ""
 * @param father {string} the father's id
 * @param mother {string} the mother's id
 * @return {Pedigree} the extended copy
 * @throws {Error} kind "pedigree" on an empty or duplicate id, an unknown sex,
 *   an unnamed parent, a member named as their own parent, or the same id for
 *   both parents
 */
export func withChild(
    ped as Pedigree,
    id as string,
    sex as string,
    father as string,
    mother as string) {
    checkNewMember($ped, $id, $sex);
    if (len($father) == 0 or len($mother) == 0) {
        fail("member '" + $id + "' must name both parents or neither");
    }
    if ($father == $id or $mother == $id) {
        fail("member '" + $id + "' cannot be their own parent");
    }
    if ($father == $mother) {
        fail("member '" + $id + "' names '" + $father + "' as both parents");
    }
    def out as Pedigree init $ped;
    $out.members[$id] = Person{id: $id, sex: $sex, father: $father, mother: $mother};
    return $out;
}

func checkNewMember(ped as Pedigree, id as string, sex as string) {
    if (len($id) == 0) {
        fail("a pedigree member needs a non-empty id");
    }
    if (maps.has($ped.members, $id)) {
        fail("pedigree '" + $ped.name + "' already has a member '" + $id + "'");
    }
    if ($sex != "m" and $sex != "f" and $sex != "") {
        fail("member '" + $id + "': sex must be \"m\", \"f\" or \"\"");
    }
}

/**
 * The members of a pedigree, in insertion order.
 * @param ped {Pedigree} the pedigree
 * @return {list of string} the member ids
 */
export func members(ped as Pedigree) {
    return maps.keys($ped.members);
}

/**
 * Whether a pedigree contains a member.
 * @param ped {Pedigree} the pedigree
 * @param id {string} the member id
 * @return {bool} true if present
 */
export func hasMember(ped as Pedigree, id as string) {
    return maps.has($ped.members, $id);
}

/**
 * One member of a pedigree.
 * @param ped {Pedigree} the pedigree
 * @param id {string} the member id
 * @return {Person} the member
 * @throws {Error} kind "pedigree" if there is no such member
 */
export func personOf(ped as Pedigree, id as string) {
    if (not maps.has($ped.members, $id)) {
        fail("pedigree '" + $ped.name + "' has no member '" + $id + "'");
    }
    return $ped.members[$id];
}

/**
 * Whether a member is a founder, whose genotype prior is the population
 * frequency.
 * @param ped {Pedigree} the pedigree
 * @param id {string} the member id
 * @return {bool} true for a founder
 * @throws {Error} kind "pedigree" if there is no such member
 */
export func isFounder(ped as Pedigree, id as string) {
    return len(personOf($ped, $id).father) == 0;
}

/**
 * The children a member has within the pedigree, in member order. An unknown
 * id has none.
 * @param ped {Pedigree} the pedigree
 * @param id {string} the member id
 * @return {list of string} the children's ids
 */
export func childrenOf(ped as Pedigree, id as string) {
    def out as list of string init [];
    for (def other in $ped.members) {
        def p as Person init $ped.members[$other];
        if ($p.father == $id or $p.mother == $id) {
            $out[] = $other;
        }
    }
    return $out;
}

# --- structural validation --------------------------------------------------

/**
 * Check that a pedigree is well formed and peelable: every named parent is a
 * member, nobody is their own ancestor, and the individual/family graph is a
 * forest. A consanguineous mating is reported as a loop, not as an ancestor
 * cycle.
 * @param ped {Pedigree} the pedigree
 * @throws {Error} kind "pedigree" naming the specific problem
 */
export func validate(ped as Pedigree) {
    for (def id in $ped.members) {
        def p as Person init $ped.members[$id];
        if (len($p.father) > 0) {
            requireMember($ped, $p.father, $id, "father");
            requireMember($ped, $p.mother, $id, "mother");
        }
    }
    checkNoAncestorCycle($ped);
    checkForest($ped);
}

func requireMember(ped as Pedigree, parent as string, child as string, role as string) {
    if (not maps.has($ped.members, $parent)) {
        fail("member '" + $child + "' names a " + $role + " '" + $parent +
            "' that is not in pedigree '" + $ped.name + "'");
    }
}

# checkNoAncestorCycle proves that nobody is their own ancestor. One
# depth-first pass over the whole parent graph settles it: a member on the
# current path is open, and meeting an open member again closes a cycle.
#
# Walking up separately from every member would revisit the same ancestors over
# and over - quadratic on a deep pedigree - where this settles each member once.
# Meeting an already-settled ancestor by a second route is a consanguinity
# diamond, not a cycle, and is left to `checkForest` to report in its own words.
func checkNoAncestorCycle(ped as Pedigree) {
    def state as map of string to int init {};
    for (def id in $ped.members) {
        if (not maps.has($state, $id)) {
            $state = settleAncestors($ped, $state, $id);
        }
    }
}

# settleAncestors walks up from `start`, marking each member open on the way up
# and settled on the way back down. `state` is 1 while a member is on the
# current path and 2 once its ancestors are all accounted for.
func settleAncestors(ped as Pedigree, state as map of string to int, start as string) {
    def out as map of string to int init $state;
    # The stack is a list plus a live-entry count. Pushing reuses a slot an
    # earlier pop left behind before it grows the list, and popping only moves
    # the count - where lists.pop would return a fresh copy of the whole stack
    # on every pass. A member stays on the stack while its parents are walked,
    # which is what closes it on the way back down.
    def stack as list of string init [$start];
    def top as int init 1;
    while ($top > 0) {
        def cur as string init $stack[$top - 1];
        if (maps.has($out, $cur)) {
            if ($out[$cur] == 1) {
                $out[$cur] = 2;
            }
            $top = $top - 1;
            continue;
        }
        $out[$cur] = 1;
        for (def p in unsettledParents($ped, $out, $cur)) {
            if ($top < len($stack)) {
                $stack[$top] = $p;
            } else {
                $stack[] = $p;
            }
            $top = $top + 1;
        }
    }
    return $out;
}

# unsettledParents lists the parents that still need walking, and raises if one
# is already open - that is the member appearing among its own ancestors.
func unsettledParents(ped as Pedigree, state as map of string to int, id as string) {
    def out as list of string init [];
    for (def p in parentIds($ped, $id)) {
        if (maps.has($state, $p)) {
            if ($state[$p] == 1) {
                fail("member '" + $p + "' is their own ancestor");
            }
        } else {
            $out[] = $p;
        }
    }
    return $out;
}

# parentIds lists a member's parents, or nothing for a founder. A parent named
# but absent is left to `validate`'s own check, which reports it precisely.
func parentIds(ped as Pedigree, id as string) {
    def p as Person init $ped.members[$id];
    if (len($p.father) == 0) {
        return [];
    }
    def out as list of string init [];
    if (maps.has($ped.members, $p.father)) {
        $out[] = $p.father;
    }
    if (maps.has($ped.members, $p.mother)) {
        $out[] = $p.mother;
    }
    return $out;
}

# Couples is the pedigree's nuclear families: the parent pairs in
# first-appearance order, plus a father -> mother -> family-number index so a
# lookup is a hash hit rather than a scan.
def struct Couples {
    pairs as list of list of string,
    index as map of string to map of string to int
};

# couple orders a parent pair so the same two people always give the same
# two-element list, whichever child names them.
func couple(father as string, mother as string) {
    if ($father < $mother) {
        return [$father, $mother];
    }
    return [$mother, $father];
}

# familyCouples groups the pedigree's children by parent couple, returning the
# couples in first-appearance order.
func familyCouples(ped as Pedigree) {
    def out as Couples init Couples{pairs: [], index: {}};
    for (def id in $ped.members) {
        def p as Person init $ped.members[$id];
        if (len($p.father) > 0) {
            $out = addCouple($out, couple($p.father, $p.mother));
        }
    }
    return $out;
}

# addCouple records a parent pair the first time it is met, numbering it in
# first-appearance order.
func addCouple(couples as Couples, pair as list of string) {
    if (hasCouple($couples, $pair[0], $pair[1])) {
        return $couples;
    }
    def out as Couples init $couples;
    def row as map of string to int init {};
    if (maps.has($out.index, $pair[0])) {
        $row = $out.index[$pair[0]];
    }
    def pairs as list of list of string init $out.pairs;
    $row[$pair[1]] = len($pairs);
    $pairs[] = $pair;
    $out.index[$pair[0]] = $row;
    $out.pairs = $pairs;
    return $out;
}

# indexOfCouple finds a parent pair's family number. The pairs are indexed by
# father then mother, so this is two hash lookups rather than a scan through
# every couple in the pedigree comparing lists.
func indexOfCouple(couples as Couples, a as string, b as string) {
    return $couples.index[$a][$b];
}

func hasCouple(couples as Couples, a as string, b as string) {
    return maps.has($couples.index, $a) and maps.has($couples.index[$a], $b);
}

# checkForest rejects any loop in the individual/family graph, using union-find
# over nodes named "i:<id>" for individuals and "f:<n>" for families. An edge
# that joins two already-connected nodes closes a loop.
#
# The unions are done in place. Passing the parent map to a helper and taking
# the result back would copy the whole map once per edge, which is quadratic in
# the pedigree and was the bulk of validation's cost.
func checkForest(ped as Pedigree) {
    def parent as map of string to string init {};
    for (def id in $ped.members) {
        $parent["i:" + $id] = "i:" + $id;
    }
    def couples as Couples init familyCouples($ped);
    for (def n in 0..len($couples.pairs)) {
        $parent["f:{$n}"] = "f:{$n}";
    }
    for (def id in $ped.members) {
        def p as Person init $ped.members[$id];
        if (len($p.father) > 0) {
            def pair as list of string init couple($p.father, $p.mother);
            def n as int init indexOfCouple($couples, $pair[0], $pair[1]);
            def ra as string init findRoot($parent, "f:{$n}");
            def rb as string init findRoot($parent, "i:" + $id);
            if ($ra == $rb) {
                failLoop($ped.name);
            }
            $parent[$ra] = $rb;
        }
    }
    for (def n in 0..len($couples.pairs)) {
        for (def side in 0..2) {
            def ra as string init findRoot($parent, "f:{$n}");
            def rb as string init findRoot($parent, "i:" + $couples.pairs[$n][$side]);
            if ($ra == $rb) {
                failLoop($ped.name);
            }
            $parent[$ra] = $rb;
        }
    }
}

func failLoop(name as string) {
    fail("pedigree '" + $name + "' contains a loop (a consanguineous mating or a repeated" +
        " marriage); Elston-Stewart peeling needs a simple pedigree");
}

# indexOfString finds a label's position, or -1.
func indexOfString(xs as list of string, want as string) {
    for (def i in 0..len($xs)) {
        if ($xs[$i] == $want) {
            return $i;
        }
    }
    return -1;
}

func findRoot(parent as map of string to string, x as string) {
    def cur as string init $x;
    while ($parent[$cur] != $cur) {
        $cur = $parent[$cur];
    }
    return $cur;
}

/**
 * The non-throwing form of `validate`.
 * @param ped {Pedigree} the pedigree
 * @return {bool} true if `validate` would succeed
 */
export func isSimple(ped as Pedigree) {
    try {
        validate($ped);
    } catch (e) {
        return false;
    }
    return true;
}

# --- per-locus frame --------------------------------------------------------

# observedAlleles lists, in first-appearance order, every allele label a typed
# member of the pedigree shows at this locus. Everything else is lumped.
func observedAlleles(ped as Pedigree, genotypes as map of string to forensics.Genotype) {
    def out as list of string init [];
    for (def id in $genotypes) {
        if (maps.has($ped.members, $id)) {
            $out = addLabels($out, $genotypes[$id]);
        }
    }
    return $out;
}

func addLabels(labels as list of string, g as forensics.Genotype) {
    def out as list of string init $labels;
    for (def a in [$g.a, $g.b]) {
        if (not lists.contains($out, $a)) {
            $out[] = $a;
        }
    }
    return $out;
}

# buildFrame assembles everything the peeling needs at one locus: the lumped
# allele frequencies, the genotype index, the nuclear families and their
# incidence, and each typed member's observed genotype index.
func buildFrame(
    ped as Pedigree,
    db as forensics.FrequencyDb,
    locus as string,
    genotypes as map of string to forensics.Genotype) {
    def labels as list of string init observedAlleles($ped, $genotypes);
    def freq as list of float init [];
    def total as float init 0.0;
    for (def a in $labels) {
        def p as float init forensics.alleleFrequency($db, $locus, $a);
        $freq[] = $p;
        $total = $total + $p;
    }
    # The observed alleles' floored frequencies can overshoot 1. Rescale rather
    # than carry a negative lumped class; otherwise the leftover mass becomes
    # the lumped allele standing for everything nobody showed.
    if ($total >= 1.0) {
        for (def i in 0..len($freq)) {
            $freq[$i] = $freq[$i] / $total;
        }
    } else {
        $freq[] = 1.0 - $total;
    }
    def gi as list of int init [];
    def gj as list of int init [];
    for (def i in 0..len($freq)) {
        for (def j in $i..len($freq)) {
            $gi[] = $i;
            $gj[] = $j;
        }
    }
    def observed as map of string to int init {};
    for (def id in $genotypes) {
        if (maps.has($ped.members, $id)) {
            def g as forensics.Genotype init $genotypes[$id];
            def ia as int init indexOfString($labels, $g.a);
            def ib as int init indexOfString($labels, $g.b);
            $observed[$id] = genotypeIndex(len($freq), $ia, $ib);
        }
    }
    return withFamilies(
        $ped,
        Frame{
            alleleFreq: $freq,
            nAlleles: len($freq),
            pairIdx: pairTable(len($freq)),
            gi: $gi,
            gj: $gj,
            famFather: [],
            famMother: [],
            famChildren: [],
            incidence: {},
            founder: {},
            observed: $observed
        });
}

# pairTable precomputes every allele pair's genotype index, so the innermost
# loops read it instead of calling `genotypeIndex` - a table lookup where a
# function call would otherwise sit inside a doubly-nested sum.
func pairTable(nAlleles as int) {
    def out as list of list of int init [];
    for (def a in 0..$nAlleles) {
        def row as list of int init [];
        for (def b in 0..$nAlleles) {
            $row[] = genotypeIndex($nAlleles, $a, $b);
        }
        $out[] = $row;
    }
    return $out;
}

# genotypeIndex maps an unordered allele pair to its position in the genotype
# list. Genotypes are enumerated as (i, j) with i <= j, row by row, so the
# position is the start of row `lo` plus the offset within it - a closed form,
# not a search.
func genotypeIndex(nAlleles as int, a as int, b as int) {
    def lo as int init math.min($a, $b);
    def hi as int init math.max($a, $b);
    return $lo * $nAlleles - ($lo * ($lo - 1)) // 2 + ($hi - $lo);
}

# withFamilies fills in the nuclear families and the individual-to-family
# incidence lists the message passing walks.
func withFamilies(ped as Pedigree, frame as Frame) {
    def out as Frame init $frame;
    for (def id in $ped.members) {
        $out.incidence[$id] = [];
        $out.founder[$id] = len($ped.members[$id].father) == 0;
    }
    def couples as Couples init familyCouples($ped);
    def fathers as list of string init [];
    def mothers as list of string init [];
    def kids as list of list of string init [];
    for (def n in 0..len($couples.pairs)) {
        $fathers[] = $couples.pairs[$n][0];
        $mothers[] = $couples.pairs[$n][1];
        $kids[] = [];
        $out = addIncidence($out, $couples.pairs[$n][0], $n);
        $out = addIncidence($out, $couples.pairs[$n][1], $n);
    }
    for (def id in $ped.members) {
        def p as Person init $ped.members[$id];
        if (len($p.father) > 0) {
            def pair as list of string init couple($p.father, $p.mother);
            def n as int init indexOfCouple($couples, $pair[0], $pair[1]);
            def row as list of string init $kids[$n];
            $row[] = $id;
            $kids[$n] = $row;
            $out = addIncidence($out, $id, $n);
        }
    }
    $out.famFather = $fathers;
    $out.famMother = $mothers;
    $out.famChildren = $kids;
    return $out;
}

func addIncidence(frame as Frame, id as string, famIdx as int) {
    def out as Frame init $frame;
    def cur as list of int init $out.incidence[$id];
    $cur[] = $famIdx;
    $out.incidence[$id] = $cur;
    return $out;
}

# --- the peeling ------------------------------------------------------------

# founderProbability is the Hardy-Weinberg prior for a founder's genotype.
func founderProbability(frame as Frame, g as int) {
    def pa as float init $frame.alleleFreq[$frame.gi[$g]];
    def pb as float init $frame.alleleFreq[$frame.gj[$g]];
    if ($frame.gi[$g] == $frame.gj[$g]) {
        return $pa * $pa;
    }
    return 2.0 * $pa * $pb;
}

# A child takes one allele from each parent, so only four genotypes are
# reachable from a parental pair - one per combination of parental slots, each
# with probability 1/4. Every sum over "all genotypes weighted by transmission"
# is therefore a sum over those four, which is what these two helpers do. It is
# exactly the transmission probability
#   P(child (x,y) | gf, gm) = T(gf,x) T(gm,y) + T(gf,y) T(gm,x)
# written from the parents' side rather than the child's, and it turns the inner
# loop from O(genotypes) into O(1).

# penetrance is 1 for an untyped member and, for a typed one, 1 only on the
# genotype actually observed.
func penetrance(frame as Frame, id as string, g as int) {
    if (not maps.has($frame.observed, $id)) {
        return 1.0;
    }
    if ($frame.observed[$id] == $g) {
        return 1.0;
    }
    return 0.0;
}

# messageToFamily summarises everything hanging off `id` other than family
# `skip` into a vector over `id`'s genotypes: its own observation, its founder
# prior if it has one, and the messages from its other families.
func messageToFamily(frame as Frame, id as string, skip as int) {
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        def v as float init penetrance($frame, $id, $g);
        if ($frame.founder[$id]) {
            $v = $v * founderProbability($frame, $g);
        }
        $out[] = $v;
    }
    for (def fam in $frame.incidence[$id]) {
        if ($fam != $skip) {
            def inbound as list of float init messageToIndividual($frame, $fam, $id);
            for (def g in 0..len($out)) {
                $out[$g] = $out[$g] * $inbound[$g];
            }
        }
    }
    return $out;
}

# messageToIndividual summarises family `fam` - its other members and
# everything hanging off them - into a vector over `id`'s genotypes. `id` is
# either one of the two parents or one of the children, and the two cases sum
# over different variables.
func messageToIndividual(frame as Frame, fam as int, id as string) {
    if ($id == $frame.famFather[$fam] or $id == $frame.famMother[$fam]) {
        return messageToParent($frame, $fam, $id);
    }
    return messageToChild($frame, $fam, $id);
}

# childWeights builds, for one nuclear family, the matrix
#
#   prod[gf][gm] = product over the given children of
#                  P(that child's whole subtree | parents gf, gm)
#
# Computing it once per family is the point. The per-child factor is a sum over
# the four genotypes a child of (gf, gm) can carry - one allele from each
# parent, a quarter each - which is the transmission probability
#   P(child (x,y) | gf, gm) = T(gf,x) T(gm,y) + T(gf,y) T(gm,x)
# written from the parents' side, and is exact.
#
# The children loop is outermost deliberately. Binding a child's message vector
# inside the doubly-nested sum copies the whole vector on every iteration, and
# that copying - not the arithmetic - was the module's dominant cost.
func childWeights(frame as Frame, kids as list of list of float) {
    def size as int init len($frame.gi);
    def prod as list of list of float init [];
    for (def gf in 0..$size) {
        def row as list of float init [];
        for (def gm in 0..$size) {
            $row[] = 1.0;
        }
        $prod[] = $row;
    }
    def out as list of list of float init $prod;
    for (def msg in $kids) {
        $out = foldChild($frame, $out, $msg);
    }
    return $out;
}

# foldChild multiplies one child's contribution into every cell of the matrix.
func foldChild(frame as Frame, prod as list of list of float, msg as list of float) {
    def out as list of list of float init $prod;
    def size as int init len($frame.gi);
    def idx as list of list of int init $frame.pairIdx;
    for (def gf in 0..$size) {
        def fa as int init $frame.gi[$gf];
        def fb as int init $frame.gj[$gf];
        for (def gm in 0..$size) {
            def ma as int init $frame.gi[$gm];
            def mb as int init $frame.gj[$gm];
            def share as float init $msg[$idx[$fa][$ma]] + $msg[$idx[$fa][$mb]] +
                $msg[$idx[$fb][$ma]] + $msg[$idx[$fb][$mb]];
            $out[$gf][$gm] = $out[$gf][$gm] * 0.25 * $share;
        }
    }
    return $out;
}

# messageToParent: sum over the spouse's genotypes, weighting each by the
# spouse's own message and by every child's summarised transmission.
func messageToParent(frame as Frame, fam as int, id as string) {
    def spouse as string init $frame.famMother[$fam];
    if ($id == $spouse) {
        $spouse = $frame.famFather[$fam];
    }
    def up as list of float init messageToFamily($frame, $spouse, $fam);
    def prod as list of list of float init childWeights($frame, childMessages($frame, $fam, ""));
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        def total as float init 0.0;
        for (def gs in 0..len($frame.gi)) {
            if ($up[$gs] > 0.0) {
                $total = $total + $up[$gs] * $prod[$g][$gs];
            }
        }
        $out[] = $total;
    }
    return $out;
}

# messageToChild: sum over both parents' genotypes, weighting each pair by the
# parents' own messages, by the child's transmission, and by every sibling's
# summarised transmission.
func messageToChild(frame as Frame, fam as int, id as string) {
    def upFather as list of float init messageToFamily($frame, $frame.famFather[$fam], $fam);
    def upMother as list of float init messageToFamily($frame, $frame.famMother[$fam], $fam);
    def prod as list of list of float init childWeights($frame, childMessages($frame, $fam, $id));
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        $out[] = 0.0;
    }
    for (def gf in 0..len($frame.gi)) {
        if ($upFather[$gf] > 0.0) {
            $out = overMotherGenotypes($frame, $out, $prod, $upFather[$gf], $upMother, $gf);
        }
    }
    return $out;
}

# overMotherGenotypes folds in every maternal genotype for one fixed paternal
# genotype, which keeps `messageToChild`'s double sum flat.
#
# The scatter is written out here rather than delegated. `acc` is mutated, so
# handing it to a helper and taking the result back would copy the whole vector
# once per maternal genotype; held here it is copied once per paternal genotype.
func overMotherGenotypes(
    frame as Frame,
    acc as list of float,
    prod as list of list of float,
    wFather as float,
    upMother as list of float,
    gf as int) {
    def out as list of float init $acc;
    def idx as list of list of int init $frame.pairIdx;
    def fa as int init $frame.gi[$gf];
    def fb as int init $frame.gj[$gf];
    for (def gm in 0..len($frame.gi)) {
        def w as float init $wFather * $upMother[$gm] * $prod[$gf][$gm];
        if ($w > 0.0) {
            # A quarter of the weight to each of the four genotypes a child of
            # this pair can carry.
            def share as float init 0.25 * $w;
            def ma as int init $frame.gi[$gm];
            def mb as int init $frame.gj[$gm];
            def one as int init $idx[$fa][$ma];
            def two as int init $idx[$fa][$mb];
            def three as int init $idx[$fb][$ma];
            def four as int init $idx[$fb][$mb];
            $out[$one] = $out[$one] + $share;
            $out[$two] = $out[$two] + $share;
            $out[$three] = $out[$three] + $share;
            $out[$four] = $out[$four] + $share;
        }
    }
    return $out;
}

# childMessages collects the inbound message of every child of `fam` except
# `exclude`.
func childMessages(frame as Frame, fam as int, exclude as string) {
    def out as list of list of float init [];
    for (def c in $frame.famChildren[$fam]) {
        if ($c != $exclude) {
            $out[] = messageToFamily($frame, $c, $fam);
        }
    }
    return $out;
}

# componentRoots picks one member from each connected component of the
# pedigree, so a pedigree that is several unconnected families is peeled once
# per family and the results multiplied.
func componentRoots(ped as Pedigree, frame as Frame) {
    def seen as map of string to bool init {};
    def roots as list of string init [];
    for (def id in $ped.members) {
        if (not maps.has($seen, $id)) {
            $roots[] = $id;
            $seen = markComponent($frame, $seen, $id);
        }
    }
    return $roots;
}

# markComponent adds every member reachable from `start` through shared
# families to `seen`.
func markComponent(frame as Frame, seen as map of string to bool, start as string) {
    def out as map of string to bool init $seen;
    # The frontier is consumed by advancing an index rather than by lists.pop,
    # which copies the whole list on every pass. Reachability does not care
    # about the order, so a queue marks the same members a stack would.
    def queue as list of string init [$start];
    def head as int init 0;
    while ($head < len($queue)) {
        def cur as string init $queue[$head];
        $head = $head + 1;
        if (not maps.has($out, $cur)) {
            $out[$cur] = true;
            for (def id in familyMembers($frame, $cur)) {
                $queue[] = $id;
            }
        }
    }
    return $out;
}

# familyMembers lists everyone sharing a nuclear family with `id` - the
# neighbours a component walk visits next.
func familyMembers(frame as Frame, id as string) {
    def out as list of string init [];
    for (def fam in $frame.incidence[$id]) {
        $out[] = $frame.famFather[$fam];
        $out[] = $frame.famMother[$fam];
        for (def c in $frame.famChildren[$fam]) {
            $out[] = $c;
        }
    }
    return $out;
}

/**
 * The likelihood of a pedigree at one locus: the probability of the observed
 * genotypes given the parentages, with every untyped member summed out.
 *
 * Members absent from `genotypes` are untyped; ids that are not members are
 * ignored. A pedigree in several unconnected parts is the product of their
 * likelihoods.
 * @param ped {Pedigree} the pedigree
 * @param genotypes {map of string to Genotype} member id -> observed genotype
 *   at this locus
 * @param db {FrequencyDb} the reference frequencies
 * @param locus {string} the locus name
 * @return {float} the likelihood at this locus
 * @throws {Error} kind "pedigree" if the pedigree is malformed, looped or
 *   empty, or a member's genotype is at another locus; kind "forensicgenetics"
 *   if the database does not cover the locus
 */
export func locusLikelihood(
    ped as Pedigree,
    genotypes as map of string to forensics.Genotype,
    db as forensics.FrequencyDb,
    locus as string) {
    validate($ped);
    if (len($ped.members) == 0) {
        fail("pedigree '" + $ped.name + "' has no members");
    }
    for (def id in $genotypes) {
        if (maps.has($ped.members, $id) and $genotypes[$id].locus != $locus) {
            fail("member '" + $id + "' has a genotype at locus '" + $genotypes[$id].locus +
                "' but the likelihood was asked for at '" + $locus + "'");
        }
    }
    def frame as Frame init buildFrame($ped, $db, $locus, $genotypes);
    def out as float init 1.0;
    for (def root in componentRoots($ped, $frame)) {
        def msg as list of float init messageToFamily($frame, $root, -1);
        def total as float init 0.0;
        for (def g in 0..len($frame.gi)) {
            $total = $total + $msg[$g];
        }
        $out = $out * $total;
    }
    return $out;
}

# --- across a panel ---------------------------------------------------------

/**
 * The loci the database covers and at least one pedigree member is typed at.
 * @param ped {Pedigree} the pedigree
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {list of string} the usable locus names, in database order
 */
export func pedigreeLoci(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def out as list of string init [];
    for (def locus in forensics.loci($db)) {
        for (def id in $profiles) {
            if (maps.has($ped.members, $id) and forensics.hasGenotype($profiles[$id], $locus) and
                not lists.contains($out, $locus)) {
                $out[] = $locus;
            }
        }
    }
    return $out;
}

# genotypesAt gathers the members typed at one locus.
func genotypesAt(ped as Pedigree, profiles as map of string to forensics.Profile, locus as string) {
    def out as map of string to forensics.Genotype init {};
    for (def id in $profiles) {
        if (maps.has($ped.members, $id) and forensics.hasGenotype($profiles[$id], $locus)) {
            $out[$id] = forensics.genotypeAt($profiles[$id], $locus);
        }
    }
    return $out;
}

/**
 * The likelihood of a pedigree across a panel: the product of
 * `locusLikelihood` over `pedigreeLoci`. A long panel can underflow; use
 * `log10PedigreeLikelihood` or `pedigreeLR`.
 * @param ped {Pedigree} the pedigree
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood
 * @throws {Error} kind "pedigree" if no locus is usable, or if the pedigree is
 *   malformed or looped
 */
export func pedigreeLikelihood(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def loci as list of string init pedigreeLoci($ped, $profiles, $db);
    if (len($loci) == 0) {
        fail("pedigree '" + $ped.name +
            "': no locus is both covered by the database and typed in a member");
    }
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = $out * locusLikelihood($ped, genotypesAt($ped, $profiles, $locus), $db, $locus);
    }
    return $out;
}

/**
 * The base-10 logarithm of `pedigreeLikelihood`, summed per locus so it does
 * not underflow.
 * @param ped {Pedigree} the pedigree
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} log10 of the likelihood
 * @throws {Error} kind "pedigree" as for `pedigreeLikelihood`, or if a locus
 *   has likelihood zero (the pedigree is impossible, and has no logarithm)
 */
export func log10PedigreeLikelihood(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def loci as list of string init pedigreeLoci($ped, $profiles, $db);
    if (len($loci) == 0) {
        fail("pedigree '" + $ped.name +
            "': no locus is both covered by the database and typed in a member");
    }
    def out as float init 0.0;
    for (def locus in $loci) {
        def l as float init locusLikelihood(
            $ped,
            genotypesAt($ped, $profiles, $locus),
            $db,
            $locus);
        if ($l <= 0.0) {
            fail("pedigree '" + $ped.name + "' has likelihood zero at locus '" + $locus +
                "', so it has no logarithm");
        }
        $out = $out + math.log10($l);
    }
    return $out;
}

/**
 * The likelihood ratio of two pedigrees over the same typed people.
 *
 * Taken locus by locus over the numerator's `pedigreeLoci`, so it stays in
 * range where either likelihood would underflow. Both pedigrees must type the
 * same members; untyped members may differ, which is how a hypothesis adds an
 * unknown alternative father.
 * @param numerator {Pedigree} the first hypothesis
 * @param denominator {Pedigree} the second hypothesis
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood ratio
 * @throws {Error} kind "pedigree" if the two pedigrees type different members,
 *   if no locus is usable, or if the denominator has likelihood zero
 */
export func pedigreeLR(
    numerator as Pedigree,
    denominator as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def typedA as list of string init typedMembers($numerator, $profiles);
    def typedB as list of string init typedMembers($denominator, $profiles);
    if (lists.sort($typedA) != lists.sort($typedB)) {
        fail("pedigrees '" + $numerator.name + "' and '" + $denominator.name +
            "' type different members, so their likelihoods are not comparable");
    }
    def loci as list of string init pedigreeLoci($numerator, $profiles, $db);
    if (len($loci) == 0) {
        fail("no locus is both covered by the database and typed in a member");
    }
    def out as float init 1.0;
    for (def locus in $loci) {
        def under as float init locusLikelihood(
            $denominator,
            genotypesAt($denominator, $profiles, $locus),
            $db,
            $locus);
        if ($under <= 0.0) {
            fail("pedigree '" + $denominator.name + "' has likelihood zero at locus '" + $locus +
                "', so no ratio is defined there");
        }
        def over as float init locusLikelihood(
            $numerator,
            genotypesAt($numerator, $profiles, $locus),
            $db,
            $locus);
        $out = $out * ($over / $under);
    }
    return $out;
}

/**
 * The members of a pedigree that have a profile with at least one genotype.
 * @param ped {Pedigree} the pedigree
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @return {list of string} the typed member ids, in pedigree order
 */
export func typedMembers(ped as Pedigree, profiles as map of string to forensics.Profile) {
    def out as list of string init [];
    for (def id in $ped.members) {
        if (maps.has($profiles, $id) and len(forensics.profileLoci($profiles[$id])) > 0) {
            $out[] = $id;
        }
    }
    return $out;
}
