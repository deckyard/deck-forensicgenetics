# `kinship`

```jennifer
import "@deckyard/forensicgenetics/kinship.j" as kinship;
```

Pairwise relationship likelihood ratios and trio paternity. It takes the
core's `Profile`, `Genotype` and `FrequencyDb`. Signatures, the relationship
and mutation constants, and struct fields are in the
[Cheatsheet](cheatsheet.md#kinship).

## Conventions

- The original functions use `θ = 0`. Each has a `...Theta` variant taking
  `theta` last (`kinshipLRTheta`, `paternityIndexWithTheta`, ...) that draws
  every non-IBD allele by the Balding-Nichols sampling formula; at `theta = 0`
  a variant returns exactly what the original returns. `duoPaternityIndex` and
  `duoPaternityIndexAt` take `model` and `theta` directly. See
  [Theta](../concepts/theta.md).
- `HALF_SIBLINGS`, `GRANDPARENT` and `AVUNCULAR` have equal coefficients and
  give equal LRs.
- `kinshipLR` and `kinshipLRAgainst` multiply over `pairLoci`;
  `paternityIndex`, `paternityIndexWith` and `inconsistentLoci` over
  `trioLoci`. A locus not every profile carries is skipped; one they all carry
  but the table does not cover throws, so an exclusion cannot be dropped
  unseen. Narrow the profiles with `forensicgenetics.restrict` to leave a
  locus out on purpose.
- `paternityIndex` uses `DEFAULT_MUTATION` (rate 0.002, decay 0.1, a
  placeholder). The model applies to the alleged father's transmission only.
  Published per-locus rates, other mutation models, a maternal model and
  silent alleles are in [`mutation`](mutation.md).
- `mutationProbability` counts `k = max(1, round(|value(s) - value(t)|))`
  steps; at rate 0 any change has probability 0 without counting steps.
- `relationship` accepts coefficients in `[0, 1]` summing to 1 within 1e-9;
  `mutationModel` and `withLocusRate` accept rates and decay in `[0, 1)`;
  `probabilityOfPaternityPrior` a prior in `(0, 1)`; the theta functions a
  theta in `[0, 1)`.

Formulae: [Formulae and sources](../project/formulae.md). Guides:
[Kinship](../casework/kinship.md), [Paternity](../casework/paternity.md).

## Errors

`Error{kind: "kinship"}`, thrown for:

- IBD coefficients outside `[0, 1]` or not summing to 1
- a mutation rate or decay outside `[0, 1)`
- a step count between alleles that are not repeat counts, when the rate is
  above 0
- a theta outside `[0, 1)` (the `...Theta` and `duo...` functions)
- a genotype pair, a duo or a trio, at different loci
- profiles sharing no locus, or sharing one the table does not cover
- an alternative hypothesis giving the observed pair probability zero
  (`kinshipLRAt` and the functions built on it)
- no shared locus for a pair, a duo or a trio
- a duo whose child genotype has probability zero in the table
- a maternal exclusion (the mother can transmit neither of the child's
  alleles)
- a negative CPI, or a prior outside `(0, 1)`

A locus the database lacks surfaces as the core's `"forensicgenetics"` error.
