# `mutation`

```jennifer
import "@deckyard/forensicgenetics/mutation.j" as mutation;
```

Familias mutation models as matrices over a locus's alleles, sex-specific
models for trio and duo paternity, silent (null) alleles in paternity and
pairwise kinship, and published per-locus mutation counts. It takes the core's
`Profile`, `Genotype` and `FrequencyDb` and `kinship.j`'s `Relationship`.
Signatures and fields are in the [Cheatsheet](cheatsheet.md#mutation); the
guide is [Mutation models and silent alleles](../casework/mutation-models.md).

## Models

| Constructor | Familias name | Off-diagonal `M[i][j]` |
| --- | --- | --- |
| `equalModel(rate)` | Equal | `rate / (n - 1)` |
| `proportionalModel(rate)` | Proportional | `alpha · p_j`, `alpha = rate / S / Σ p(1-p)` |
| `stepwiseModel(rate, range)` | Stepwise ("probability decreasing with range (equal)" in Familias 2) | `rate · range^k / Σ range^k'` within `i`'s microvariant class, 0 across classes |
| `extendedStepwiseModel(rate, range, rate2)` | Extended stepwise (Familias 3) | as Stepwise within the class; `rate2 / (alleles in other classes)` across |
| `stationaryModel(model)` | Stabilization "PM" | each row of `model` off the diagonal scaled by `d_i` |

`n` and `p` cover the detectable alleles; `S` is their summed frequency
(`1 - q` with a silent allele). Formulae, the diagonal of each model and the
PM scaling are in [Formulae and sources](../project/formulae.md#familias-mutation-models).

## Conventions

- A matrix is built per locus over the table's alleles (insertion order), then
  typed alleles the table lacks (at the `5/(2N)` floor), then
  `SILENT_ALLELE` when the locus has a silent frequency above 0. Equal and
  Proportional therefore depend on how many alleles the table lists.
- Stepwise needs every detectable allele to be a repeat count with at most
  one decimal (`15`, `15.2`), no `.0` decimal (`15.0` is written `15`), and
  no two labels with one repeat value (`10` and `010`): such a pair would be
  two alleles zero steps apart. Repeat values above `1e6` throw. The
  microvariant class is FamiliasLocus's `microgroup`:
  `round(10 (v - round(v)))` with R's `round` (halves to even), so `15.7` and
  `15.3` are different classes, and so are `15.5` (class -5) and `14.5` /
  `16.5` (class 5).
- PM needs at most `MAX_STATIONARY_ALLELES` (200) detectable alleles.
- The silent allele never mutates and is never mutated to; its row and column
  are those of the identity.
- `rate` (and a `perLocus` override) is the integer-step rate for Stepwise;
  `rate2` has no per-locus override.
- `ParentageModel.paternal` acts on the alleged father in a trio or duo,
  `maternal` on the mother in a trio. The random man of the alternative does
  not mutate. The duo does not use `maternal`: the untyped mother's allele is
  drawn from the population.
- With a silent frequency `q`, detectable frequencies are the table's floored
  values times `1 - q`; a single-peak genotype `a` is `a/a` or `a/silent`.
- `theta` applies the Balding-Nichols sampling formula as in `kinship.j`'s
  `...Theta` functions, over the alleles of the true genotypes being summed.
- With `MENDELIAN` (or every rate 0) and no silent allele, `paternityIndexAt`,
  `duoPaternityIndexAt` and `kinshipLRAt` return exactly what
  `kinship.paternityIndexAtTheta`, `kinship.duoPaternityIndexAt` (under
  `kinship.NO_MUTATION`) and `kinship.kinshipLRAtTheta` return.
- `aabb2003Counts` names loci as `panels.normalizeLocus` does (`vWA`,
  `Penta D`, `SE33`). `withLocusRates` matches names exactly; normalise the
  profiles and table with `panels.normalizeProfile` / `panels.normalizeDb`
  first if they spell loci differently.

## Errors

`Error{kind: "mutation"}`, thrown for:

- a rate or `rate2` outside `[0, 1)`, a range not above 0, `rate + rate2`
  above 1, a silent frequency outside `[0, 1)`, theta outside `[0, 1)`
- an unknown model name in a hand-built `MutationSpec`, or a rate in force at
  the locus (`rate` or its `perLocus` override) outside `[0, 1)`; for
  Stepwise also `rate2`, `range` and `rate + rate2` (checked before a zero
  rate short-cuts to the identity)
- Stepwise over an allele that is not a repeat count with at most one decimal,
  a `.0` label, a repeat value above `1e6`, two labels with one repeat value (the message names both), a
  `range` whose `range^k` for the largest step `k` within a microvariant
  class falls outside about `1e±300` (it would underflow to 0 or overflow), or
  a hand-built `MutationSpec` with `range` not above 0
- a matrix entry outside `[0, 1]` or a row not summing to 1 (a Proportional
  rate too high for the frequencies)
- PM over more than `MAX_STATIONARY_ALLELES` detectable alleles
- a PM transformation that does not exist (a scaled diagonal below 0, common
  with a rare allele) or a matrix with no unique stationary distribution
  (Stepwise with microvariants and `rate2` 0)
- a typed genotype carrying `SILENT_ALLELE`; a locus the table does not cover
- genotypes at different loci; profiles sharing no locus, or sharing one the
  table does not cover (never dropped silently)
- a child genotype impossible given the mother (a maternal exclusion when the
  maternal model has no mutation), or of probability zero in the population
  (duo); an alternative giving the observations probability zero
- a CPI or LR product past the float range (`kinshipLR` names
  `log10KinshipLR`); a locus ratio of 0 in `log10KinshipLR`
- an allele not in the matrix (`transitionProbability`), or a matrix and a
  `LocusAlleles` listing different alleles (`overallRate`, `isStationary`)
- `aabb2003Rates` with anything but `"total"`, `"paternal"` or `"maternal"`

`kinshipLRAt` with `silent` 0 delegates to `kinship.kinshipLRAtTheta` and
throws its `"kinship"` errors (mismatched loci, a zero alternative), and so do
`kinshipLR` / `log10KinshipLR` at such a locus;
`kinshipModelWithRates` throws `"kinship"` on a bad rate. A missing allele
frequency is the core's `"forensicgenetics"` error.
