# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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

/**
 * The most distinct alleles the typed members may show at one locus, excluding
 * the lumped class.
 *
 * A nuclear family whose parents are both untyped costs children * G^2, with
 * G = A (A + 1) / 2 genotypes over A = alleles + 1. At this limit one untyped
 * couple with ten typed children takes about 3 s per locus.
 */
export def const MAX_PEDIGREE_ALLELES as int init 20;

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
    checkParents($id, $father, $mother);
    def out as Pedigree init $ped;
    $out.members[$id] = Person{id: $id, sex: $sex, father: $father, mother: $mother};
    return $out;
}

/**
 * A pedigree built from its members in one pass.
 *
 * A person with `""` for both parents is a founder and is checked as by
 * `withFounder`; anyone else as by `withChild`. Folding `withFounder` /
 * `withChild` copies the whole pedigree per member, which is quadratic (23 s
 * for 4000 members); this is linear.
 * @param name {string} a label for the pedigree
 * @param people {list of Person} the members, in member order
 * @return {Pedigree} the pedigree
 * @throws {Error} kind "pedigree" as `withFounder` and `withChild` do
 */
export func pedigreeOf(name as string, people as list of Person) {
    def out as Pedigree init Pedigree{name: $name, members: {}};
    for (def p in $people) {
        checkNewMember($out, $p.id, $p.sex);
        if (len($p.father) > 0 or len($p.mother) > 0) {
            checkParents($p.id, $p.father, $p.mother);
        }
        $out.members[$p.id] = $p;
    }
    return $out;
}

# checkParents rejects a child's parent pair: one of them unnamed, the child
# as its own parent, or one id for both.
func checkParents(id as string, father as string, mother as string) {
    if (len($father) == 0 or len($mother) == 0) {
        fail("member '" + $id + "' must name both parents or neither");
    }
    if ($father == $id or $mother == $id) {
        fail("member '" + $id + "' cannot be their own parent");
    }
    if ($father == $mother) {
        fail("member '" + $id + "' names '" + $father + "' as both parents");
    }
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
    def p as Person init personOf($ped, $id);
    return len($p.father) == 0 and len($p.mother) == 0;
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

# Shape is a validated pedigree in index form, built once per call by
# `shapeOf` and shared by every locus. Members are numbered in member order.
# Nodes of the individual/family graph are numbered members first (0..n-1),
# then nuclear families (n..n+F-1).
#
# A family's two parents are `parentA` and `parentB` in lexicographic id order,
# not father and mother: autosomal transmission is symmetric, and a fixed
# order fixes the summation order of the peeling.
#
# `order` lists every node, component by component, each component in
# depth-first pre-order from its root; `up` is a node's neighbour towards the
# root (-1 for a root). Every node comes after its `up`, so walking `order`
# backwards has every message a node needs ready before the node itself.
def struct Shape {
    ids as list of string,
    index as map of string to int,
    founder as list of bool,
    parentA as list of int,
    parentB as list of int,
    children as list of list of int,
    incidence as list of list of int,
    order as list of int,
    up as list of int,
    roots as list of int
};

/**
 * Check that a pedigree is well formed and peelable: every named parent is a
 * member, nobody is their own ancestor, and the individual/family graph is a
 * forest. A consanguineous mating is reported as a loop, not as an ancestor
 * cycle.
 * @param ped {Pedigree} the pedigree
 * @throws {Error} kind "pedigree" naming the specific problem
 */
export func validate(ped as Pedigree) {
    shapeOf($ped);
}

func requireMember(ped as Pedigree, parent as string, child as string, role as string) {
    if (not maps.has($ped.members, $parent)) {
        fail("member '" + $child + "' names a " + $role + " '" + $parent +
            "' that is not in pedigree '" + $ped.name + "'");
    }
}

# Parentage is each member's parents by member number (-1 for a founder).
def struct Parentage {
    founder as list of bool,
    father as list of int,
    mother as list of int
};

# Families is the nuclear families, numbered in first-appearance order, and
# the family each member is a child of (-1 for a founder), and each member's
# incidence list (see `incidenceOf`). Helpers take it whole: a struct argument
# that is only read is borrowed, where a bare list argument is copied per call.
def struct Families {
    parentA as list of int,
    parentB as list of int,
    children as list of list of int,
    childFam as list of int,
    incidence as list of list of int
};

# Walk is the traversal order of the individual/family forest (see `Shape`).
def struct Walk {
    order as list of int,
    up as list of int,
    roots as list of int
};

# shapeOf validates a pedigree and returns it in index form. Every structure is
# a local list or map written in place and every helper only reads its struct
# arguments (borrowed, not copied), so the pass is linear in members plus
# families, with a log factor from union-find.
func shapeOf(ped as Pedigree) {
    def ids as list of string init maps.keys($ped.members);
    def index as map of string to int init {};
    for (def i in 0..len($ids)) {
        $index[$ids[$i]] = $i;
    }
    def par as Parentage init parentage($ped, $ids, $index);
    checkNoAncestorCycle($ids, $par);
    def fams as Families init families($ids, $par);
    $fams.incidence = incidenceOf(len($ids), $fams);
    checkForest($ped.name, len($ids), $fams);
    def w as Walk init walk(len($ids), $fams);
    return Shape{
        ids: $ids,
        index: $index,
        founder: $par.founder,
        parentA: $fams.parentA,
        parentB: $fams.parentB,
        children: $fams.children,
        incidence: $fams.incidence,
        order: $w.order,
        up: $w.up,
        roots: $w.roots
    };
}

# parentage resolves every member's parents to member numbers, failing on a
# parent who is not a member. A `Person` literal skips `withChild`'s checks, so
# they are repeated here: one named parent would otherwise make a founder and
# drop the other parent's link (and any exclusion it carries).
func parentage(ped as Pedigree, ids as list of string, index as map of string to int) {
    def out as Parentage init Parentage{founder: [], father: [], mother: []};
    for (def id in $ids) {
        def p as Person init $ped.members[$id];
        if (len($p.father) == 0 and len($p.mother) == 0) {
            $out.founder[] = true;
            $out.father[] = -1;
            $out.mother[] = -1;
            continue;
        }
        checkParents($id, $p.father, $p.mother);
        requireMember($ped, $p.father, $id, "father");
        requireMember($ped, $p.mother, $id, "mother");
        $out.founder[] = false;
        $out.father[] = $index[$p.father];
        $out.mother[] = $index[$p.mother];
    }
    return $out;
}

# families groups the children by parent couple. A couple's key is its two
# member numbers in lexicographic id order, so both parents' roles and the
# order the children name them in give the same family.
func families(ids as list of string, par as Parentage) {
    def out as Families init Families{
        parentA: [],
        parentB: [],
        children: [],
        childFam: [],
        incidence: []
    };
    def famIndex as map of string to int init {};
    for (def i in 0..len($ids)) {
        if ($par.founder[$i]) {
            $out.childFam[] = -1;
            continue;
        }
        def a as int init $par.father[$i];
        def b as int init $par.mother[$i];
        if ($ids[$b] < $ids[$a]) {
            $a = $par.mother[$i];
            $b = $par.father[$i];
        }
        def key as string init "{$a},{$b}";
        if (not maps.has($famIndex, $key)) {
            $famIndex[$key] = len($out.parentA);
            $out.parentA[] = $a;
            $out.parentB[] = $b;
            $out.children[] = [];
        }
        def f as int init $famIndex[$key];
        $out.children[$f][] = $i;
        $out.childFam[] = $f;
    }
    return $out;
}

# incidenceOf lists each member's families: those it parents, in family order,
# then the one it is a child of. The order fixes the order inbound messages
# are multiplied in.
func incidenceOf(n as int, fams as Families) {
    def out as list of list of int init [];
    for (def i in 0..$n) {
        $out[] = [];
    }
    for (def f in 0..len($fams.parentA)) {
        $out[$fams.parentA[$f]][] = $f;
        $out[$fams.parentB[$f]][] = $f;
    }
    for (def i in 0..$n) {
        if ($fams.childFam[$i] >= 0) {
            $out[$i][] = $fams.childFam[$i];
        }
    }
    return $out;
}

# walk orders the individual/family forest depth-first from the first
# unvisited member of each component. The graph is known to be a forest, so
# every node is pushed once and a stack of n + F slots suffices.
func walk(n as int, fams as Families) {
    def total as int init $n + len($fams.parentA);
    def out as Walk init Walk{order: [], up: [], roots: []};
    def seen as list of bool init [];
    def stack as list of int init [];
    for (def k in 0..$total) {
        $seen[] = false;
        $out.up[] = -1;
        $stack[] = 0;
    }
    for (def s in 0..$n) {
        if ($seen[$s]) {
            continue;
        }
        $out.roots[] = $s;
        $seen[$s] = true;
        $stack[0] = $s;
        def top as int init 1;
        while ($top > 0) {
            $top = $top - 1;
            def cur as int init $stack[$top];
            $out.order[] = $cur;
            for (def v in neighbours($n, $fams, $cur)) {
                if ($seen[$v]) {
                    continue;
                }
                $seen[$v] = true;
                $out.up[$v] = $cur;
                $stack[$top] = $v;
                $top = $top + 1;
            }
        }
    }
    return $out;
}

# neighbours lists a node's neighbours in the individual/family graph: a
# member's families, or a family's two parents and its children.
func neighbours(n as int, fams as Families, node as int) {
    def out as list of int init [];
    if ($node < $n) {
        for (def f in $fams.incidence[$node]) {
            $out[] = $n + $f;
        }
        return $out;
    }
    def f as int init $node - $n;
    $out[] = $fams.parentA[$f];
    $out[] = $fams.parentB[$f];
    for (def c in $fams.children[$f]) {
        $out[] = $c;
    }
    return $out;
}

# checkNoAncestorCycle proves that nobody is their own ancestor. One
# depth-first pass over the whole parent graph settles it: `state` is 1 while a
# member is on the current path and 2 once its ancestors are all accounted
# for, and meeting an open member again closes a cycle.
#
# Walking up separately from every member would revisit the same ancestors over
# and over - quadratic on a deep pedigree - where this settles each member once.
# Meeting an already-settled ancestor by a second route is a consanguinity
# diamond, not a cycle, and is left to `checkForest` to report in its own words.
#
# The stack is a list plus a live-entry count: popping only moves the count,
# where lists.pop would copy the whole stack. A member stays on the stack while
# its parents are walked, which is what closes it on the way back down. Each
# member opens once and pushes at most two parents, so 2n + 1 slots suffice.
func checkNoAncestorCycle(ids as list of string, par as Parentage) {
    def state as list of int init [];
    def stack as list of int init [0];
    for (def i in 0..len($ids)) {
        $state[] = 0;
        $stack[] = 0;
        $stack[] = 0;
    }
    for (def s in 0..len($ids)) {
        if ($state[$s] != 0) {
            continue;
        }
        $stack[0] = $s;
        def top as int init 1;
        while ($top > 0) {
            def cur as int init $stack[$top - 1];
            if ($state[$cur] != 0) {
                if ($state[$cur] == 1) {
                    $state[$cur] = 2;
                }
                $top = $top - 1;
                continue;
            }
            $state[$cur] = 1;
            for (def q in parentsOf($par, $cur)) {
                if ($state[$q] == 1) {
                    fail("member '" + $ids[$q] + "' is their own ancestor");
                }
                if ($state[$q] == 0) {
                    $stack[$top] = $q;
                    $top = $top + 1;
                }
            }
        }
    }
}

# parentsOf lists a member's two parents, or nothing for a founder.
func parentsOf(par as Parentage, i as int) {
    if ($par.founder[$i]) {
        return [];
    }
    return [$par.father[$i], $par.mother[$i]];
}

# UnionFind is checkForest's disjoint-set forest: each node's parent and each
# root's tree size.
def struct UnionFind {
    parent as list of int,
    size as list of int
};

# checkForest rejects any loop in the individual/family graph by union-find
# over its nodes (members 0..n-1, families n..n+F-1). An edge that joins two
# already-connected nodes closes a loop. Union by size keeps every tree
# O(log n) deep, so `findRoot` needs no path compression and can borrow the
# forest read-only. Edges run child-family first, then family-parent;
# whether a loop exists does not depend on the order.
func checkForest(name as string, n as int, fams as Families) {
    def uf as UnionFind init UnionFind{parent: [], size: []};
    for (def k in 0..($n + len($fams.parentA))) {
        $uf.parent[] = $k;
        $uf.size[] = 1;
    }
    def edgeA as list of int init [];
    def edgeB as list of int init [];
    for (def i in 0..$n) {
        if ($fams.childFam[$i] >= 0) {
            $edgeA[] = $n + $fams.childFam[$i];
            $edgeB[] = $i;
        }
    }
    for (def f in 0..len($fams.parentA)) {
        $edgeA[] = $n + $f;
        $edgeB[] = $fams.parentA[$f];
        $edgeA[] = $n + $f;
        $edgeB[] = $fams.parentB[$f];
    }
    for (def e in 0..len($edgeA)) {
        def ra as int init findRoot($uf, $edgeA[$e]);
        def rb as int init findRoot($uf, $edgeB[$e]);
        if ($ra == $rb) {
            failLoop($name);
        }
        if ($uf.size[$ra] < $uf.size[$rb]) {
            $uf.parent[$ra] = $rb;
            $uf.size[$rb] = $uf.size[$rb] + $uf.size[$ra];
        } else {
            $uf.parent[$rb] = $ra;
            $uf.size[$ra] = $uf.size[$ra] + $uf.size[$rb];
        }
    }
}

func failLoop(name as string) {
    fail("pedigree '" + $name + "' contains a loop (a consanguineous mating or a repeated" +
        " marriage); Elston-Stewart peeling needs a simple pedigree");
}

func findRoot(uf as UnionFind, x as int) {
    def cur as int init $x;
    while ($uf.parent[$cur] != $cur) {
        $cur = $uf.parent[$cur];
    }
    return $cur;
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

# Frame is the per-locus part of the peeling, built by `buildFrame`.
#
# Genotypes are indexed, not named: `gi[k]` and `gj[k]` are the two allele
# indices of genotype k, with gi <= gj so each unordered pair appears once.
# Allele indices run over the alleles actually observed at this locus, plus -
# when their frequencies leave room - one final lumped class standing for
# every allele nobody showed. `observed[i]` is member i's genotype index, or -1
# when untyped.
def struct Frame {
    alleleFreq as list of float,
    nAlleles as int,
    pairIdx as list of list of int,
    gi as list of int,
    gj as list of int,
    observed as list of int
};

# SUM_TOLERANCE absorbs the rounding of a floored-frequency sum: a table whose
# observed alleles sum to exactly 1 can add up to 1.0000000000000002.
def const SUM_TOLERANCE as float init 0.000000001;

# observedAlleles lists, in first-appearance order, every allele label a typed
# member of the pedigree shows at this locus. Everything else is lumped.
func observedAlleles(shape as Shape, genotypes as map of string to forensics.Genotype) {
    def out as list of string init [];
    for (def id in $genotypes) {
        if (maps.has($shape.index, $id)) {
            for (def a in [$genotypes[$id].a, $genotypes[$id].b]) {
                if (not lists.contains($out, $a)) {
                    $out[] = $a;
                }
            }
        }
    }
    return $out;
}

# buildFrame assembles the per-locus part of the peeling: the lumped allele
# frequencies, the genotype index and each typed member's genotype index.
func buildFrame(
    shape as Shape,
    db as forensics.FrequencyDb,
    locus as string,
    genotypes as map of string to forensics.Genotype) {
    def labels as list of string init observedAlleles($shape, $genotypes);
    def shown as int init len($labels);
    if ($shown > MAX_PEDIGREE_ALLELES) {
        fail("locus '" + $locus + "': the typed members show {$shown} distinct alleles," +
            " more than MAX_PEDIGREE_ALLELES ({MAX_PEDIGREE_ALLELES})");
    }
    def freq as list of float init [];
    def total as float init 0.0;
    for (def a in $labels) {
        def p as float init forensics.alleleFrequency($db, $locus, $a);
        $freq[] = $p;
        $total = $total + $p;
    }
    # The leftover mass becomes the lumped allele standing for everything
    # nobody showed. On a small table the observed alleles' floored
    # frequencies can overshoot 1, and then there is no distribution to sum an
    # untyped member over: kept as they are, an untyped founder's prior sums
    # to more than 1 and the ratio comes to depend on how many untyped people
    # each hypothesis names; rescaled to sum to 1, a floored allele drops back
    # below its 5/(2N) floor and this module stops agreeing with kinship.j on
    # the same pair. Neither is an answer, so the locus is refused. A sum past
    # 1 by no more than SUM_TOLERANCE is rounding: the lumped class is then
    # empty (clamped to 0) and left out.
    if ($total > 1.0 + SUM_TOLERANCE) {
        def n as int init forensics.sampleSize($db, $locus);
        fail("locus '" + $locus + "': the observed alleles' floored frequencies sum to {$total}," +
            " past 1, so there is no distribution to sum untyped members over - the table" +
            " (N = {$n}) is too small for this pedigree at this locus; leave the locus out" +
            " with forensics.restrict, or use a larger table");
    }
    if ($total < 1.0) {
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
    def observed as list of int init [];
    for (def i in 0..len($shape.ids)) {
        $observed[] = -1;
    }
    for (def id in $genotypes) {
        if (maps.has($shape.index, $id)) {
            def ia as int init indexOfString($labels, $genotypes[$id].a);
            def ib as int init indexOfString($labels, $genotypes[$id].b);
            $observed[$shape.index[$id]] = genotypeIndex(len($freq), $ia, $ib);
        }
    }
    return Frame{
        alleleFreq: $freq,
        nAlleles: len($freq),
        pairIdx: pairTable(len($freq)),
        gi: $gi,
        gj: $gj,
        observed: $observed
    };
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
# is therefore a sum over those four. It is exactly the transmission probability
#   P(child (x,y) | gf, gm) = T(gf,x) T(gm,y) + T(gf,y) T(gm,x)
# written from the parents' side rather than the child's, and it turns the inner
# loop from O(genotypes) into O(1).

# Messages holds the peeling's message vectors, indexed by node. It is a struct
# because a struct parameter that is only read is borrowed, while a bare list
# argument is deep-copied on every call - which made peeling quadratic.
def struct Messages {
    v as list of list of float
};

# peel is the likelihood at one locus: the product over components of the sum
# of the root's message. It runs the message passing bottom-up over
# `shape.order` instead of recursing, so a pedigree's depth costs no call
# stack. `msgs[node]` is the message from a node to its `up` neighbour (for a
# root, its whole component's vector), each computed exactly once.
func peel(shape as Shape, frame as Frame) {
    def n as int init len($shape.ids);
    def msgs as Messages init Messages{v: []};
    for (def k in 0..len($shape.order)) {
        $msgs.v[] = [];
    }
    def k as int init len($shape.order) - 1;
    while ($k >= 0) {
        def node as int init $shape.order[$k];
        if ($node < $n) {
            $msgs.v[$node] = messageToFamily($shape, $frame, $msgs, $node);
        } else {
            $msgs.v[$node] = messageToIndividual(
                $shape,
                $frame,
                $msgs,
                $node - $n,
                $shape.up[$node]);
        }
        $k = $k - 1;
    }
    def out as float init 1.0;
    for (def root in $shape.roots) {
        def total as float init 0.0;
        for (def g in 0..len($frame.gi)) {
            $total = $total + $msgs.v[$root][$g];
        }
        $out = $out * $total;
    }
    return $out;
}

# messageToFamily summarises everything hanging off member `i` other than its
# `up` family into a vector over `i`'s genotypes: its own observation, its
# founder prior if it has one, and the messages from its other families.
func messageToFamily(shape as Shape, frame as Frame, msgs as Messages, i as int) {
    def n as int init len($shape.ids);
    def obs as int init $frame.observed[$i];
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        def v as float init 1.0;
        if ($obs >= 0 and $obs != $g) {
            $v = 0.0;
        }
        if ($shape.founder[$i]) {
            $v = $v * founderProbability($frame, $g);
        }
        $out[] = $v;
    }
    for (def fam in $shape.incidence[$i]) {
        def node as int init $n + $fam;
        if ($node != $shape.up[$i]) {
            for (def g in 0..len($out)) {
                $out[$g] = $out[$g] * $msgs.v[$node][$g];
            }
        }
    }
    return $out;
}

# messageToIndividual summarises family `fam` - its other members and
# everything hanging off them - into a vector over member `target`'s genotypes.
# `target` is either one of the two parents or one of the children, and the two
# cases sum over different variables.
func messageToIndividual(
    shape as Shape,
    frame as Frame,
    msgs as Messages,
    fam as int,
    target as int) {
    def a as int init $shape.parentA[$fam];
    def b as int init $shape.parentB[$fam];
    if ($target == $a or $target == $b) {
        return messageToParent($shape, $frame, $msgs, $fam, $target);
    }
    return messageToChild($shape, $frame, $msgs, $fam, $target);
}

# childWeights builds, for one nuclear family, the matrix
#
#   w[r * len(cols) + c] = product over the children other than `exclude` of
#                          P(that child's whole subtree | parents rows[r], cols[c])
#
# `rows` and `cols` are the genotypes of the first and the second parent the
# caller reads; a typed parent's message is non-zero at one genotype, so only
# those cells are filled (a family with both parents typed is one cell, not
# G^2). The per-child factor is a sum over the four genotypes a child of
# (gf, gm) can carry - one allele from each parent, a quarter each - which is
# the transmission probability written from the parents' side, and is exact.
# The children's messages are read in place from `msgs`.
func childWeights(
    shape as Shape,
    frame as Frame,
    msgs as Messages,
    fam as int,
    exclude as int,
    rows as list of int,
    cols as list of int) {
    def nc as int init len($cols);
    def prod as list of float init [];
    for (def k in 0..(len($rows) * $nc)) {
        $prod[] = 1.0;
    }
    for (def c in $shape.children[$fam]) {
        if ($c != $exclude) {
            for (def r in 0..len($rows)) {
                def fa as int init $frame.gi[$rows[$r]];
                def fb as int init $frame.gj[$rows[$r]];
                for (def k in 0..$nc) {
                    def ma as int init $frame.gi[$cols[$k]];
                    def mb as int init $frame.gj[$cols[$k]];
                    def share as float init $msgs.v[$c][$frame.pairIdx[$fa][$ma]] +
                        $msgs.v[$c][$frame.pairIdx[$fa][$mb]] +
                        $msgs.v[$c][$frame.pairIdx[$fb][$ma]] +
                        $msgs.v[$c][$frame.pairIdx[$fb][$mb]];
                    def cell as int init $r * $nc + $k;
                    $prod[$cell] = $prod[$cell] * 0.25 * $share;
                }
            }
        }
    }
    return $prod;
}

# support lists the genotypes at which node `node`'s message is non-zero, in
# ascending order.
func support(msgs as Messages, node as int) {
    def out as list of int init [];
    for (def g in 0..len($msgs.v[$node])) {
        if ($msgs.v[$node][$g] > 0.0) {
            $out[] = $g;
        }
    }
    return $out;
}

# messageToParent: sum over the spouse's genotypes, weighting each by the
# spouse's own message and by every child's summarised transmission. The
# target's genotype is the first parent in `childWeights` whichever parent it
# is (transmission is symmetric), and only the spouse's non-zero genotypes are
# filled and summed.
func messageToParent(shape as Shape, frame as Frame, msgs as Messages, fam as int, target as int) {
    def spouse as int init $shape.parentB[$fam];
    if ($target == $spouse) {
        $spouse = $shape.parentA[$fam];
    }
    def all as list of int init [];
    for (def g in 0..len($frame.gi)) {
        $all[] = $g;
    }
    def cols as list of int init support($msgs, $spouse);
    def nc as int init len($cols);
    def prod as list of float init childWeights($shape, $frame, $msgs, $fam, -1, $all, $cols);
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        def total as float init 0.0;
        for (def k in 0..$nc) {
            $total = $total + $msgs.v[$spouse][$cols[$k]] * $prod[$g * $nc + $k];
        }
        $out[] = $total;
    }
    return $out;
}

# messageToChild: sum over both parents' genotypes, weighting each pair by the
# parents' own messages, by the child's transmission, and by every sibling's
# summarised transmission. A quarter of each pair's weight goes to each of the
# four genotypes a child of that pair can carry. Pairs where either parent's
# message is zero contribute nothing and are skipped.
func messageToChild(shape as Shape, frame as Frame, msgs as Messages, fam as int, target as int) {
    def a as int init $shape.parentA[$fam];
    def b as int init $shape.parentB[$fam];
    def rows as list of int init support($msgs, $a);
    def cols as list of int init support($msgs, $b);
    def nc as int init len($cols);
    def prod as list of float init childWeights($shape, $frame, $msgs, $fam, $target, $rows, $cols);
    def out as list of float init [];
    for (def g in 0..len($frame.gi)) {
        $out[] = 0.0;
    }
    for (def r in 0..len($rows)) {
        def gf as int init $rows[$r];
        def wFather as float init $msgs.v[$a][$gf];
        def fa as int init $frame.gi[$gf];
        def fb as int init $frame.gj[$gf];
        for (def k in 0..$nc) {
            def gm as int init $cols[$k];
            def w as float init $wFather * $msgs.v[$b][$gm] * $prod[$r * $nc + $k];
            if ($w > 0.0) {
                def share as float init 0.25 * $w;
                def ma as int init $frame.gi[$gm];
                def mb as int init $frame.gj[$gm];
                def one as int init $frame.pairIdx[$fa][$ma];
                def two as int init $frame.pairIdx[$fa][$mb];
                def three as int init $frame.pairIdx[$fb][$ma];
                def four as int init $frame.pairIdx[$fb][$mb];
                $out[$one] = $out[$one] + $share;
                $out[$two] = $out[$two] + $share;
                $out[$three] = $out[$three] + $share;
                $out[$four] = $out[$four] + $share;
            }
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
 *   empty, a member's genotype is at another locus, the typed members show
 *   more than `MAX_PEDIGREE_ALLELES` alleles, or the observed alleles' floored
 *   frequencies sum past 1; kind "forensicgenetics" if a member is typed and
 *   the database does not cover the locus (with nobody typed the likelihood
 *   is 1)
 */
export func locusLikelihood(
    ped as Pedigree,
    genotypes as map of string to forensics.Genotype,
    db as forensics.FrequencyDb,
    locus as string) {
    def shape as Shape init shapeOf($ped);
    if (len($ped.members) == 0) {
        fail("pedigree '" + $ped.name + "' has no members");
    }
    for (def id in $genotypes) {
        if (maps.has($ped.members, $id) and $genotypes[$id].locus != $locus) {
            fail("member '" + $id + "' has a genotype at locus '" + $genotypes[$id].locus +
                "' but the likelihood was asked for at '" + $locus + "'");
        }
    }
    return peel($shape, buildFrame($shape, $db, $locus, $genotypes));
}

# likelihoodAt is `locusLikelihood` for a pedigree already validated into
# `shape`, with genotypes from `genotypesAt` (all at `locus`).
func likelihoodAt(
    shape as Shape,
    genotypes as map of string to forensics.Genotype,
    db as forensics.FrequencyDb,
    locus as string) {
    return peel($shape, buildFrame($shape, $db, $locus, $genotypes));
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

# requireCovered throws when a member is typed at a locus the database does
# not cover. Dropping it would hide whatever that locus says about the
# pedigree - an impossible transmission, most dangerously - behind a figure
# built from the others, with no sign in the output. The caller narrows the
# profiles with forensics.restrict if the omission is deliberate.
func requireCovered(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    for (def id in $ped.members) {
        if (maps.has($profiles, $id)) {
            def locus as string init uncoveredLocus($profiles[$id], $db);
            if (len($locus) > 0) {
                fail("member '" + $id + "' of pedigree '" + $ped.name + "' is typed at locus '" +
                    $locus + "', which frequency database '" + $db.name + "' does not cover;" +
                    " it will not be dropped silently - narrow the profiles with" +
                    " forensics.restrict if that is intended");
            }
        }
    }
}

# uncoveredLocus names the first locus of a profile the database lacks, or ""
# when it covers them all.
func uncoveredLocus(p as forensics.Profile, db as forensics.FrequencyDb) {
    for (def locus in forensics.profileLoci($p)) {
        if (not forensics.hasLocus($db, $locus)) {
            return $locus;
        }
    }
    return "";
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
 * @throws {Error} kind "pedigree" if no locus is usable, if a member is typed
 *   at a locus the database does not cover, if the pedigree is malformed or
 *   looped, or as `locusLikelihood` at any locus
 */
export func pedigreeLikelihood(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def loci as list of string init panelLoci($ped, $profiles, $db);
    def shape as Shape init shapeOf($ped);
    def out as float init 1.0;
    for (def locus in $loci) {
        $out = $out * likelihoodAt($shape, genotypesAt($ped, $profiles, $locus), $db, $locus);
    }
    return $out;
}

# panelLoci is `pedigreeLoci` after the checks every panel function makes
# first: no member typed at an uncovered locus, and at least one usable locus.
func panelLoci(
    ped as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    requireCovered($ped, $profiles, $db);
    def loci as list of string init pedigreeLoci($ped, $profiles, $db);
    if (len($loci) == 0) {
        fail("pedigree '" + $ped.name +
            "': no locus is both covered by the database and typed in a member");
    }
    return $loci;
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
    def loci as list of string init panelLoci($ped, $profiles, $db);
    def shape as Shape init shapeOf($ped);
    def out as float init 0.0;
    for (def locus in $loci) {
        def l as float init likelihoodAt($shape, genotypesAt($ped, $profiles, $locus), $db, $locus);
        if ($l <= 0.0) {
            fail("pedigree '" + $ped.name + "' has likelihood zero at locus '" + $locus +
                "', so it has no logarithm");
        }
        $out = $out + math.log10($l);
    }
    return $out;
}

# FLOAT_MAX is the largest finite float. The interpreter raises a runtime error
# on a quotient or product past it, so `pedigreeLR` checks first.
def const FLOAT_MAX as float init 1.7976931348623157e308;

# LocusPair is one locus's likelihood under each hypothesis.
def struct LocusPair {
    locus as string,
    over as float,
    under as float
};

# ratioLoci validates both pedigrees once and returns, per usable locus, the
# numerator's and the denominator's likelihood. The denominator is validated
# and computed first, so its errors take precedence.
func ratioLoci(
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
    requireCovered($numerator, $profiles, $db);
    def loci as list of string init pedigreeLoci($numerator, $profiles, $db);
    if (len($loci) == 0) {
        fail("no locus is both covered by the database and typed in a member");
    }
    def shapeB as Shape init shapeOf($denominator);
    def shapeA as Shape init shapeOf($numerator);
    def out as list of LocusPair init [];
    for (def locus in $loci) {
        def under as float init likelihoodAt(
            $shapeB,
            genotypesAt($denominator, $profiles, $locus),
            $db,
            $locus);
        if ($under <= 0.0) {
            fail("pedigree '" + $denominator.name + "' has likelihood zero at locus '" + $locus +
                "', so no ratio is defined there");
        }
        def over as float init likelihoodAt(
            $shapeA,
            genotypesAt($numerator, $profiles, $locus),
            $db,
            $locus);
        $out[] = LocusPair{locus: $locus, over: $over, under: $under};
    }
    return $out;
}

/**
 * The likelihood ratio of two pedigrees over the same typed people.
 *
 * Taken locus by locus over the numerator's `pedigreeLoci`, so it stays finite
 * where either likelihood alone would underflow. The product of the per-locus
 * ratios can still exceed the float range on a long panel; that throws, and
 * `log10PedigreeLR` gives the logarithm instead. Both pedigrees must type the
 * same members; untyped members may differ, which is how a hypothesis adds an
 * unknown alternative father.
 * @param numerator {Pedigree} the first hypothesis
 * @param denominator {Pedigree} the second hypothesis
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} the likelihood ratio
 * @throws {Error} kind "pedigree" if the two pedigrees type different members,
 *   if a typed member carries a locus the database does not cover, if no
 *   locus is usable, if the denominator has likelihood zero, if either
 *   pedigree is malformed or looped, or if the ratio exceeds the float range
 */
export func pedigreeLR(
    numerator as Pedigree,
    denominator as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def out as float init 1.0;
    for (def pair in ratioLoci($numerator, $denominator, $profiles, $db)) {
        # Neither the quotient nor the product may leave the float range: an
        # under-1 divisor bounds the first check's product, an over-1
        # accumulator the second's quotient, so both checks are finite.
        if (($pair.under < 1.0 and $pair.over > FLOAT_MAX * $pair.under) or
            ($out > 1.0 and $pair.over / $pair.under > FLOAT_MAX / $out)) {
            fail("the likelihood ratio of pedigrees '" + $numerator.name + "' and '" +
                $denominator.name + "' exceeds the float range (above 1.8e308); use" +
                " log10PedigreeLR");
        }
        $out = $out * ($pair.over / $pair.under);
    }
    return $out;
}

/**
 * The base-10 logarithm of `pedigreeLR`, summed per locus so it neither
 * overflows nor underflows where the ratio would.
 * @param numerator {Pedigree} the first hypothesis
 * @param denominator {Pedigree} the second hypothesis
 * @param profiles {map of string to Profile} member id -> that member's profile
 * @param db {FrequencyDb} the reference frequencies
 * @return {float} log10 of the likelihood ratio
 * @throws {Error} kind "pedigree" as for `pedigreeLR`, or if the numerator has
 *   likelihood zero at a locus (the ratio is 0, which has no logarithm)
 */
export func log10PedigreeLR(
    numerator as Pedigree,
    denominator as Pedigree,
    profiles as map of string to forensics.Profile,
    db as forensics.FrequencyDb) {
    def out as float init 0.0;
    for (def pair in ratioLoci($numerator, $denominator, $profiles, $db)) {
        if ($pair.over <= 0.0) {
            fail("pedigree '" + $numerator.name + "' has likelihood zero at locus '" + $pair.locus +
                "', so the ratio is 0 and has no logarithm");
        }
        $out = $out + (math.log10($pair.over) - math.log10($pair.under));
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
