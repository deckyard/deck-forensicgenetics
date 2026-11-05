# `pedigree`

```jennifer
import "@deckyard/forensicgenetics/pedigree.j" as pedigree;
```

Whole-pedigree likelihoods by Elston-Stewart peeling. It takes the core's
`Profile`, `Genotype` and `FrequencyDb`. Signatures and struct fields are in
the [Cheatsheet](cheatsheet.md#pedigree).

## Conventions

- `θ = 0`; simple (loop-free) pedigrees only.
- A founder has `""` for both parents; everyone else names both. Parents may
  be added after their child, but must exist when the pedigree is used.
  `validate` repeats `withChild`'s parent checks, so `Person` literals get
  them too.
- `pedigreeOf(name, people)` builds in one pass; folding `withFounder` /
  `withChild` copies the pedigree per member (quadratic).
- `sex` is `"m"`, `"f"` or `""` and does not affect the result.
- Members without a genotype or profile are untyped and summed out; ids that
  are not members are ignored. A member with an empty profile is untyped.
- `pedigreeLoci` lists, in database order, the loci the table covers and at
  least one member is typed at. `pedigreeLR` and `log10PedigreeLR` use the
  numerator's. A member typed at a locus the table does not cover makes
  `pedigreeLikelihood`, `log10PedigreeLikelihood`, `pedigreeLR` and
  `log10PedigreeLR` throw rather than drop it; narrow the profiles with
  `forensicgenetics.restrict` to leave a locus out on purpose.
- Unconnected parts of a pedigree multiply.
- Each call validates each pedigree once, in `O((members + families) log
  members)` (union-find by size); peeling is iterative, so depth costs no
  call stack.
- At most `MAX_PEDIGREE_ALLELES` (20) distinct alleles per locus among the
  typed members, plus the lumped class.
- The observed alleles' floored frequencies may exceed 1 by up to `1e-9`
  (float rounding of an exact 1); the lumped class is then empty.
- Work per locus is `O(edges × G²)`, `G` being the number of genotype classes
  after lumping the alleles no typed member shows; a typed parent contributes
  one genotype, not `G`. Timings are in
  [Limitations](../project/limitations.md#peeling-slows-as-observed-alleles-accumulate).

Guide: [General pedigrees](../casework/pedigrees.md). Formulae:
[Pedigree likelihood](../project/formulae.md#pedigree-likelihood).

## Errors

`Error{kind: "pedigree"}`, thrown for:

- an empty or duplicate member id, or an unknown sex
- one parent named without the other, a member named as their own parent, or
  one id as both parents
- an unknown member in `personOf` or `isFounder`
- a named parent that is not a member
- an ancestor cycle (someone is their own ancestor)
- a loop in the individual/family graph; a consanguineous mating is reported
  as a loop, not as an ancestor cycle
- a typed member carrying a locus the table does not cover, in the panel
  functions
- more than `MAX_PEDIGREE_ALLELES` distinct alleles typed at one locus
- a locus where the observed alleles' floored frequencies sum past 1, so
  there is no distribution to sum untyped members over
- an empty pedigree
- a genotype recorded at a locus other than the one asked for
- no usable locus
- a zero likelihood where a logarithm or a ratio's denominator was needed
- a `pedigreeLR` past the float maximum (use `log10PedigreeLR`)
- two pedigrees in `pedigreeLR` that type different members

A locus the database lacks surfaces as the core's `"forensicgenetics"` error
when a member is typed there; with nobody typed, `locusLikelihood` returns 1.
