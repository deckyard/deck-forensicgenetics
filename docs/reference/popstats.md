# `popstats`

```jennifer
import "@deckyard/forensicgenetics/popstats.j" as popstats;
```

Population statistics of a reference table and of the genotype sample behind
it. Takes the core's `FrequencyDb` and `Genotype`. Signatures, constants and
struct fields are in the [Cheatsheet](cheatsheet.md#popstats).

## Conventions

- Frequency statistics use the listed frequencies without the `5/(2N)` floor,
  rescaled to sum to 1. A locus summing below `MIN_FREQUENCY_SUM` (0.95)
  throws.
- Every frequency statistic assumes Hardy-Weinberg proportions; the
  `observed...` functions do not.
- `unbiasedHeterozygosity` needs the locus's `n`; with `n == 0` it throws.
- Combined values assume independent loci. `combinedStats` takes
  `powerOfDiscrimination` as `1 - Π MP`, so the combined PD keeps the
  precision of the match probabilities.
- A sample is a `list of forensics.Genotype`, all at one locus, one entry per
  person. `withEstimatedLocus` sets the locus's `n` to the sample size.
- `hweExactTest` is exact for one or two alleles and throws on more;
  `hwePermutationTest` takes any number. Both return 1.0 for a monomorphic
  sample. Probabilities within a relative 1e-7 count as ties.
- `hwePermutationTest` uses its own Park-Miller generator; `math.rand` is not
  touched. Seed in `[1, 2147483646]`, permutations in
  `[1, MAX_PERMUTATIONS]`.

- F-statistics take one genotype sample per population at one locus
  (`fStatistics`) or profiles grouped by population (`fStatisticsTable`,
  which leaves out a profile untyped at the locus). `fStatistics` returns the
  components `a`, `b`, `c`; `fst`, `fit` and `fis` are their ratios, and the
  multilocus estimate is the ratio of sums (`combineFStatistics`), not a mean
  of per-locus ratios. Estimates may be negative.
- With one population only `fis` is defined (`a` is 0); `fst` and `fit`
  throw. A ratio whose denominator is within `1e-12` of zero throws:
  every population monomorphic for the same allele (all three), or every
  population fixed (`fis`).
- `ldPermutationTest` takes two aligned genotype lists, one entry per
  individual; `ldPermutationTests` pairs loci in name order over the
  individuals typed at both and runs every pair on one generator stream from
  `seed`. A pair with fewer than two individuals, or one genotype at either
  locus, gets 1.0. No multiple-testing correction is applied.

Formulae: [Population statistics](../project/formulae.md#population-statistics).
Guide: [Population statistics](../concepts/population-statistics.md),
[Population structure and linkage](../concepts/population-structure.md).

## Errors

`Error{kind: "popstats"}`, thrown for:

- a locus the table lacks, a locus with no alleles, or one summing below
  `MIN_FREQUENCY_SUM`
- `unbiasedHeterozygosity` on a locus with `n == 0`
- an empty table, an empty `combinedStats` / `combinedPower` /
  `log10CombinedComplement` input, a power outside `[0, 1]`, or a power of 1
  in `log10CombinedComplement`
- an empty sample, a sample mixing loci, or an empty allele label
- `observedTypicalPaternityIndex` on a sample with no homozygote
- `hweExactTest` on more than two alleles
- a permutation count or seed out of range
- `fStatistics` with no samples, an empty sample, samples at different loci,
  or a mean sample size of one individual or less
- `fst` / `fit` with one population; `fst` / `fit` / `fis` on a zero
  denominator
- `fStatisticsTable` with no loci, or a population nobody in which is typed
  at a locus; `combineFStatistics` on an empty table or mixed population
  counts
- `ldPermutationTest` on samples of different lengths or at the same locus;
  `ldPermutationTests` with no profiles or fewer than two loci

`withEstimatedLocus` builds through `forensics.withLocus`; its checks raise the
core's `"forensicgenetics"` kind.
