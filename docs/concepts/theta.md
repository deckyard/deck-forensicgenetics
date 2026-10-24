# Theta and the sub-population correction

The product rule assumes one randomly mating population. Real populations are
structured into sub-populations, and two people from the same one share
alleles more often than the overall frequencies predict. Ignoring that makes
the match probability too small and the evidence look stronger than it is.

Theta (`θ`, also written Fst) is the co-ancestry coefficient: the probability
that two alleles drawn from the same sub-population are identical by descent.
A larger theta gives a larger, more conservative match probability.

`genotypeFrequency` applies NRC II recommendation 4.2, the Balding-Nichols
sampling formula conditioned on one copy of the genotype already observed. The
formulae are in [Formulae and sources](../project/formulae.md#single-locus-genotype-frequency);
at `θ = 0` they reduce to `pa²` and `2 pa pb`.

```jennifer
forensicgenetics.genotypeFrequency($db, $g, $theta);
```

Theta must be in `[0, 1)`. At `θ = 1` every allele is identical by descent and
the formula is meaningless, although it still divides.

## Which value

| Constant | Value | When |
| --- | --- | --- |
| `THETA_NONE` | 0.0 | plain product rule; an idealised panmictic population, not casework |
| `THETA_GENERAL` | 0.01 | a broad, well-mixed population (NRC II) |
| `THETA_ISOLATED` | 0.03 | small, isolated or endogamous populations (NRC II) |

Any float in `[0, 1)` is accepted.

The eight-locus profile in `examples/matchprobability.j`, against the
synthetic table:

| theta | log10 RMP | about 1 in |
| --- | --- | --- |
| 0.00 | -16.209 | 16 000 000 000 000 000 |
| 0.01 | -14.792 | 620 000 000 000 000 |
| 0.03 | -13.095 | 12 000 000 000 000 |

Three orders of magnitude on the same profile and table. Report theta with
every number it produced.

## Where theta is applied

Applied: `genotypeFrequency`, `randomMatchProbability`,
`log10RandomMatchProbability`, `matchLikelihoodRatio` and their `Across` forms.

Not applied, by design:

- `inclusionProbability`, `cpi`, `cpe`. PI counts who cannot be excluded; it
  is not a co-ancestry-corrected match probability. See
  [Mixtures](../casework/mixtures.md).
- [`kinship.j`](../casework/kinship.md) and
  [`pedigree.j`](../casework/pedigrees.md). The correction does not extend to
  the IBD-1 term by substituting an allele frequency; it needs the full
  recursive sampling formula. Applying theta to some terms and not others
  would answer neither question, so both modules use `θ = 0` throughout. This
  is listed in [Limitations](../project/limitations.md), and matches the
  default of other relationship-testing software.

## Theta and the 5/(2N) floor

Both make an estimate more conservative, for different reasons. The floor
handles sampling error in a finite reference table (see
[Reference frequency data](frequency-data.md)); theta handles population
structure. They compose: `alleleFrequency` applies the floor, and
`genotypeFrequency` applies theta to the floored values.
