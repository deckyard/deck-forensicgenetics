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
- `sex` is `"m"`, `"f"` or `""` and does not affect the result.
- Members without a genotype or profile are untyped and summed out; ids that
  are not members are ignored. A member with an empty profile is untyped.
- `pedigreeLoci` lists, in database order, the loci the table covers and at
  least one member is typed at. `pedigreeLR` uses the numerator's.
- Unconnected parts of a pedigree multiply.
- Work per locus is `O(edges × G²)`, `G` being the number of genotype classes
  after lumping the alleles no typed member shows. Timings are in
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
- an empty pedigree
- a genotype recorded at a locus other than the one asked for
- no usable locus
- a zero likelihood where a logarithm or a ratio's denominator was needed
- two pedigrees in `pedigreeLR` that type different members

A locus the database lacks surfaces as the core's `"forensicgenetics"` error.
