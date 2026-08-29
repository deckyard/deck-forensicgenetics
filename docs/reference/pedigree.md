# `pedigree`

```jennifer
import "@jennifer/forensicgenetics/pedigree.j" as pedigree;
```

Whole-pedigree likelihoods by Elston-Stewart peeling. Builds on `forensicgenetics`.
Errors are `Error{kind: "pedigree"}`.

Uses the plain product rule - `θ = 0`. Simple (loop-free) pedigrees only.

## Types

```jennifer
def struct Person {
    id as string,
    sex as string,        # "m", "f" or ""; recorded, not used by the maths
    father as string,     # "" for a founder
    mother as string      # "" for a founder
};

def struct Pedigree {
    name as string,
    members as map of string to Person
};
```

## Building

| Signature | Returns |
| --- | --- |
| `pedigree(name as string)` | `Pedigree` - empty |
| `withFounder(ped, id as string, sex as string)` | `Pedigree` - a copy with a parentless member |
| `withChild(ped, id as string, sex as string, father as string, mother as string)` | `Pedigree` - a copy; **both** parents required |

Parents need not exist yet - a pedigree can be built in any order - but must
exist by the time it is used. An empty or duplicate id, an unknown sex, one
parent named without the other, a member as their own parent, or the same id as
both parents all throw.

## Inspecting

| Signature | Returns |
| --- | --- |
| `members(ped)` | `list of string` - ids, insertion order |
| `hasMember(ped, id as string)` | `bool` |
| `personOf(ped, id as string)` | `Person` - throws if absent |
| `isFounder(ped, id as string)` | `bool` |
| `childrenOf(ped, id as string)` | `list of string` |
| `typedMembers(ped, profiles as map of string to forensics.Profile)` | `list of string` - members with a non-empty profile |

## Validation

| Signature | Returns |
| --- | --- |
| `validate(ped)` | nothing - throws, naming the problem |
| `isSimple(ped)` | `bool` - the non-throwing form |

`validate` checks that every named parent exists, that nobody is their own
ancestor, and that the individual/family graph is a **forest**. A consanguineous
mating or any other loop is rejected. Every likelihood function calls it first.

## Likelihoods

| Signature | Returns |
| --- | --- |
| `locusLikelihood(ped, genotypes as map of string to forensics.Genotype, db, locus as string)` | `float` |
| `pedigreeLoci(ped, profiles as map of string to forensics.Profile, db)` | `list of string` - loci the database covers and at least one member is typed at |
| `pedigreeLikelihood(ped, profiles, db)` | `float` - the product across loci |
| `log10PedigreeLikelihood(ped, profiles, db)` | `float` - the sum of logarithms |
| `pedigreeLR(numerator, denominator, profiles, db)` | `float` - locus-by-locus ratio |

In `locusLikelihood`, members absent from `genotypes` are untyped and summed
out; ids not in the pedigree are ignored. A genotype recorded at a *different*
locus is an error - scoring it against this locus's frequencies would return a
plausible number for a question nobody asked.

`log10PedigreeLikelihood` throws if a locus has likelihood zero - an impossible
pedigree has no logarithm.

`pedigreeLR` requires both pedigrees to have the **same typed members**, and
throws otherwise: two hypotheses scored against different data have no
meaningful ratio. Untyped members may differ freely, which is how a hypothesis
introduces an unknown alternative parent. It also throws if the denominator has
likelihood zero at some locus.

## Performance

Work is `O(edges × G²)` per locus, where `G` is the number of genotype classes
after **allele lumping** - alleles no member shows are pooled into one class,
which is exact because untyped members are unobserved. A four-typed-person,
eight-locus pedigree LR runs in well under a second.

## Errors

`Error{kind: "pedigree"}` for: an empty, duplicate or self-referential member;
one parent named without the other; an unknown sex; a parent not in the
pedigree; an ancestor cycle; a loop in the individual/family graph; a genotype
supplied for the wrong locus; an empty pedigree; no usable locus; a zero
likelihood where a logarithm or a ratio was asked for; and mismatched typed
members between two pedigrees.

A consanguineous mating is reported as a **loop**, not as an ancestor cycle -
the two are different mistakes and the message says which one it found.

See [General pedigrees](../casework/pedigrees.md) for the worked example.
