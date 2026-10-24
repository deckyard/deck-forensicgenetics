# `forensicgenetics`

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
```

The entry module: profiles, reference frequencies, match probability and
mixtures. Signatures, struct fields and constant values are in the
[Cheatsheet](cheatsheet.md#forensicgenetics); this page lists the rules the
functions apply.

## Conventions

- Allele labels are strings, compared verbatim. `genotype` orders a pair by
  repeat count when both alleles are numeric, lexicographically otherwise.
- Builders (`withLocus`, `withGenotype`) return a copy and leave the argument
  unchanged.
- Loci and map keys keep insertion order; `commonLoci` follows profile order.
- `alleleFrequency` floors at `5/(2N)` and caps at 1.0; at `n = 0` there is no
  floor. `genotypeFrequency` and the PI sum are capped at 1.0.
- `theta` must be in `[0, 1)`. The mixture functions take none.
- Panel functions use every locus of the profile or mixture and throw on one
  the table lacks; narrow with `commonLoci` and `restrict`.
- The `Across` weights must be non-negative, one per table, and sum to 1
  within 0.001.

Background: [Profiles and genotypes](../concepts/profiles.md),
[Reference frequency data](../concepts/frequency-data.md),
[Theta](../concepts/theta.md), [Match probability](../casework/match-probability.md),
[Mixtures](../casework/mixtures.md).

## Errors

`Error{kind: "forensicgenetics"}`, thrown for:

- a locus the database lacks, or a profile untyped at a locus asked for
  (`genotypeAt`, `restrict`)
- in `withLocus`: an empty locus name, a negative `n`, a frequency outside
  `(0, 1]`, or frequencies summing past `MAX_FREQUENCY_SUM`
- an unlisted allele at a locus with `n = 0` (no frequency and no floor)
- an empty locus name or allele label in `genotype`; a locus listed twice in
  `profile`
- a non-numeric allele passed to `alleleValue`
- a theta outside `[0, 1)`
- an empty profile; a mixture with no loci, or a locus with no alleles
- a profile and mixture sharing no locus (`isIncluded`)
- in `combineFrequencies`: an empty list, no shared locus, or a shared locus
  with `n = 0` in every input
- in the `Across` functions: an empty list, a weight count that differs from
  the table count, a negative weight, or weights not summing to 1
- a zero RMP in a likelihood-ratio function (only through underflow)
