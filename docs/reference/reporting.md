# `reporting`

```jennifer
import "@deckyard/forensicgenetics/reporting.j" as reporting;
```

Reporting a single-source match: allowances for a frequency database's finite
size, and the ENFSI (2015) verbal equivalent of a likelihood ratio. Takes the
core's `FrequencyDb`, `Genotype` and `Profile`. Signatures, constants and
struct fields are in the [Cheatsheet](cheatsheet.md#reporting).

## Conventions

- Counts are read back from the table as `x = p · 2N`, unrounded; an
  unlisted allele has `x = 0`. Every profile locus needs `n > 0`.
- No `5/(2N)` floor: the size-bias pseudo-counts and the Dirichlet prior are
  the allowance for rare alleles.
- Theta enters through NRC II 4.2 at the adjusted frequencies, as in
  `forensicgenetics.genotypeFrequency`; theta in `[0, 1)`.
- The size-bias correction adds two carriers of the matching genotype to the
  database: `(x + 2) / (2N + 4)` per heterozygote allele, `(x + 4) / (2N + 4)`
  for a homozygote.
- The posterior functions give each locus `Dir(x_i + α)`, `α` = 1
  (`PRIOR_UNIFORM`) or `1/k` (`PRIOR_ONE_OVER_K`, `k` the listed alleles plus
  any profile allele the table lacks). The count a partial table leaves short
  of `2N` is one further class with no prior. Loci are independent.
- Quantiles are nearest-rank: the `ceil(q · draws)`-th smallest of `draws`
  simulated log10 RMPs. A seed in `[1, 2147483646]` fixes the result;
  `math.rand` is not touched. Draws in `[1, MAX_DRAWS]`.
- `verbalEquivalent` bands include their upper bound (`(10, 100]`), and the
  mirrored bands below LR = 1 their lower bound (`[0.01, 0.1)`). LRs in
  `(1, 2)`, which the guideline's table does not list, are weak support.

Formulae: [Sampling uncertainty](../project/formulae.md#sampling-uncertainty-of-a-match-probability)
and [Verbal equivalents](../project/formulae.md#verbal-equivalents-of-a-likelihood-ratio).
Guide: [Reporting a match](../casework/reporting.md).

## Errors

`Error{kind: "reporting"}`, thrown for:

- a profile with no genotypes, a locus the table lacks, or a locus with
  `n == 0`
- theta outside `[0, 1)`
- a prior other than the two constants, an empty quantile list or a quantile
  outside `(0, 1)`, a draw count or seed out of range
- `sizeBiasLikelihoodRatio` when the product underflows to zero;
  `likelihoodRatioLowerBound` when the bound passes `1e300`
- `verbalEquivalent` on an LR that is zero or negative
