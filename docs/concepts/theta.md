# Theta and the sub-population correction

## The problem it fixes

The product rule - multiply `2 pa pb` across loci - assumes the reference
population is one big randomly-mating pool. Real populations are not. They are
made of sub-populations that mate within themselves more often than across,
and two people drawn from the same one share alleles more often than the
overall frequencies predict.

Ignore that and the match probability comes out too small: the profile looks
rarer than it is, and the evidence looks stronger than it is. Theta is the
correction.

## What theta is

Theta (`θ`, also written `Fst`) is the **co-ancestry coefficient**: the
probability that two alleles drawn from the same sub-population are identical
by descent rather than by coincidence. It is a single number standing in for
the whole population structure.

Larger theta means more shared ancestry, which means a larger - more
conservative - match probability.

## The formulae

NRC II recommendation 4.2, the Balding-Nichols sampling formula, conditioned on
one copy of the genotype already being observed. That conditioning is what
makes it a *match* probability rather than a population frequency:

```
homozygote  aa:  [2θ + (1-θ)pa] · [3θ + (1-θ)pa] / [(1+θ)(1+2θ)]

heterozygote ab: 2 · [θ + (1-θ)pa] · [θ + (1-θ)pb] / [(1+θ)(1+2θ)]
```

At `θ = 0` both collapse to the familiar Hardy-Weinberg forms, `pa²` and
`2 pa pb`. The deck implements exactly this:

```jennifer
forensics.genotypeFrequency($db, $g, $theta);
```

Theta must be in `[0, 1)`. The interval is half-open because at `θ = 1` every
allele is identical by descent and the formula stops meaning anything, even
though it still divides.

## Which value

Three named constants, because these are the choices that actually get made:

| Constant | Value | When |
| --- | --- | --- |
| `THETA_NONE` | 0.0 | The plain product rule. Correct only for an idealised panmictic population - for illustration, not casework. |
| `THETA_GENERAL` | 0.01 | A broad, well-mixed population. The usual NRC II recommendation. |
| `THETA_ISOLATED` | 0.03 | Small, isolated or endogamous populations. |

Any float in `[0, 1)` works; the constants just name the conventional ones.

## What it does to a number

From `examples/matchprobability.j`, an eight-locus profile against the bundled
synthetic table:

| theta | log10 RMP | about 1 in |
| --- | --- | --- |
| 0.00 | -16.209 | 16 000 000 000 000 000 |
| 0.01 | -14.792 | 620 000 000 000 000 |
| 0.03 | -13.095 | 12 000 000 000 000 |

Three orders of magnitude between no correction and the isolated-population
value, on the same profile and the same table. Theta is not a rounding
detail - it is a stated assumption about the population, and it belongs in
anything you write down alongside the number it produced.

## Where theta is applied - and where it is not

**Applied:** the single-source match probability.
`genotypeFrequency`, `randomMatchProbability`, `log10RandomMatchProbability`
and `matchLikelihoodRatio` all take a theta and use it.

**Not applied - by design:**

- **`inclusionProbability`, `cpi`, `cpe`.** The probability of inclusion counts
  who cannot be excluded; it is not a co-ancestry-corrected match probability,
  and applying theta to it would be mixing two different questions. See
  [Mixtures](../casework/mixtures.md).
- **[`kinship.j`](../casework/kinship.md) and
  [`pedigree.j`](../casework/pedigrees.md).** These use the plain product rule
   throughout. The sub-population correction does not extend to the IBD-1
   term - the joint probability of two related people's
  genotypes under co-ancestry needs the full recursive sampling formula, not a
  swapped-in allele frequency. Rather than apply theta to some terms and not
  others and produce a number that is neither one thing nor the other, those
  modules are uniformly theta-free and say so.

This is a real limitation, not a preference; it is listed as one in
[Limitations](../project/limitations.md). It also matches what other
relationship-testing software does by default.

## Theta and the frequency floor are different things

Both make an estimate more conservative, and they are easy to confuse:

- the **5/2N floor** handles *sampling* - the reference table is finite, so a
  rare allele's frequency is not well measured. See
  [Reference frequency data](frequency-data.md).
- **theta** handles *structure* - the population is not one pool, so alleles
  correlate within it.

They compose: `alleleFrequency` applies the floor, and `genotypeFrequency`
applies theta to whatever the floor returned.
