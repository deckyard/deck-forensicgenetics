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
- `bayesianFrequency` is FSR-GUI-0012 section 11.4 literally, `k` the listed
  alleles (also for an unlisted allele): `(x + 1)/(2N + k)` (`PRIOR_UNIFORM`),
  `(x + 1/k)/(2N + 1)` (`PRIOR_ONE_OVER_K`, 11.4.1),
  `(x + 1/(k + 2))/(2N + 1)` (`PRIOR_ONE_OVER_K_PLUS_TWO`, 11.4.2). The
  `bayesian...` RMP / LR functions evaluate NRC II 4.2 at those values.
- Theta enters through NRC II 4.2 at the adjusted frequencies, as in
  `forensicgenetics.genotypeFrequency`; theta in `[0, 1)`.
- The size-bias correction adds two carriers of the matching genotype to the
  database: `(x + 2) / (2N + 4)` per heterozygote allele, `(x + 4) / (2N + 4)`
  for a homozygote.
- The posterior functions give each locus `Dir(x_i + α)`, `α` = 1
  (`PRIOR_UNIFORM`), `1/k'` (`PRIOR_ONE_OVER_K`, `k'` the listed alleles plus
  any profile allele the table lacks) or `1/(k + 2)`
  (`PRIOR_ONE_OVER_K_PLUS_TWO`, `k` the listed alleles; an unlisted profile
  allele takes one of the two unobserved slots, unused slots join the rest
  class). Both non-uniform priors have total mass 1, so an allele's posterior
  mean is `(x + α)/(2N + 1)`. Under `PRIOR_ONE_OVER_K` that is section
  11.4.1 when every profile allele is listed (`k' = k`) and 11.4.2 when both
  alleles of a heterozygote are unlisted (`k' = k + 2`); with one unlisted
  allele it is neither. The count a partial table leaves short
  of `2N` is one further class with no prior; a table summing above 1 has
  its counts scaled to total `2N` (`bayesianFrequency` does not scale). Loci
  are independent.
- Quantiles are nearest-rank: the `ceil(q · draws)`-th smallest of `draws`
  simulated log10 RMPs, `q · draws` computed in decimal from `q`'s shortest
  decimal form (0.07 · 100 is rank 7, not the float product's 8). A seed in `[1, 2147483646]` fixes the result;
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
- a prior other than the three `PRIOR_*` constants, a locus listing no
  allele (`bayesian...`), an empty quantile list or a quantile
  outside `(0, 1)`, a draw count or seed out of range
- `sizeBiasLikelihoodRatio` / `bayesianLikelihoodRatio` when the RMP is
  below `1e-300` (the ratio passes `1e300`; use the `log10...` variant);
  `likelihoodRatioLowerBound` when the bound passes `1e300`
- `verbalEquivalent` on an LR that is zero or negative
