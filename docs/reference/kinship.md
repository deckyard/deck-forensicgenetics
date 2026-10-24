# `kinship`

```jennifer
import "@deckyard/forensicgenetics/kinship.j" as kinship;
```

Pairwise relationship likelihood ratios and trio paternity. It takes the
core's `Profile`, `Genotype` and `FrequencyDb`. Signatures, the relationship
and mutation constants, and struct fields are in the
[Cheatsheet](cheatsheet.md#kinship).

## Conventions

- `θ = 0` throughout. See [Theta](../concepts/theta.md).
- `HALF_SIBLINGS`, `GRANDPARENT` and `AVUNCULAR` have equal coefficients and
  give equal LRs.
- `kinshipLR` and `kinshipLRAgainst` multiply over `pairLoci`;
  `paternityIndex`, `paternityIndexWith` and `inconsistentLoci` over
  `trioLoci`. Other loci are ignored.
- `paternityIndex` uses `DEFAULT_MUTATION` (rate 0.002, decay 0.1, a
  placeholder). The model applies to the alleged father's transmission only.
- `mutationProbability` counts `k = max(1, round(|value(s) - value(t)|))`
  steps; at rate 0 any change has probability 0 without counting steps.
- `relationship` accepts coefficients in `[0, 1]` summing to 1 within 1e-9;
  `mutationModel` and `withLocusRate` accept rates and decay in `[0, 1)`;
  `probabilityOfPaternityPrior` a prior in `(0, 1)`.

Formulae: [Formulae and sources](../project/formulae.md). Guides:
[Kinship](../casework/kinship.md), [Paternity](../casework/paternity.md).

## Errors

`Error{kind: "kinship"}`, thrown for:

- IBD coefficients outside `[0, 1]` or not summing to 1
- a mutation rate or decay outside `[0, 1)`
- a step count between alleles that are not repeat counts, when the rate is
  above 0
- a genotype pair, or a trio, at different loci
- an alternative hypothesis giving the observed pair probability zero
  (`kinshipLRAt` and the functions built on it)
- no shared locus for a pair or a trio
- a maternal exclusion (the mother can transmit neither of the child's
  alleles)
- a negative CPI, or a prior outside `(0, 1)`

A locus the database lacks surfaces as the core's `"forensicgenetics"` error.
