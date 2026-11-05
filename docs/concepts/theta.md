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
`log10RandomMatchProbability`, `matchLikelihoodRatio` and their `Across` forms;
the mixture LR in `mixture.j`, through the sampling formula conditioned on the
known contributors (see [Mixtures](../casework/mixtures.md#the-model)).

Applied on request, in [`kinship.j`](../casework/kinship.md): the `...Theta`
variants of the kinship and trio paternity functions, and the duo paternity
index. The original `kinship.j` functions keep `θ = 0`, and each variant
returns exactly the original's value at `theta = 0`:

| `θ = 0` function | with theta |
| --- | --- |
| `jointGenotypeProbability(db, ga, gb, rel)` | `jointGenotypeProbabilityTheta(db, ga, gb, rel, theta)` |
| `kinshipLRAt(db, ga, gb, rel, alt)` | `kinshipLRAtTheta(db, ga, gb, rel, alt, theta)` |
| `kinshipLR(a, b, rel, db)` | `kinshipLRTheta(a, b, rel, db, theta)` |
| `kinshipLRAgainst(a, b, rel, alt, db)` | `kinshipLRAgainstTheta(a, b, rel, alt, db, theta)` |
| `paternityIndexAt(db, mother, child, father, model)` | `paternityIndexAtTheta(db, mother, child, father, model, theta)` |
| `paternityIndex(mother, child, allegedFather, db)` | `paternityIndexTheta(mother, child, allegedFather, db, theta)` |
| `paternityIndexWith(mother, child, allegedFather, db, model)` | `paternityIndexWithTheta(mother, child, allegedFather, db, model, theta)` |
| - | `duoPaternityIndexAt(db, child, father, model, theta)`, `duoPaternityIndex(child, allegedFather, db, model, theta)` |

```jennifer
def lr as float init kinship.kinshipLRTheta($a, $b, kinship.FULL_SIBLINGS, $db,
    forensicgenetics.THETA_GENERAL);
```

Not applied, by design:

- `inclusionProbability`, `cpi`, `cpe`. PI counts who cannot be excluded; it
  is not a co-ancestry-corrected match probability. See
  [Mixtures](../casework/mixtures.md).
- [`pedigree.j`](../casework/pedigrees.md). Peeling with theta would need the
  sampling formula across every founder allele of the pedigree at once, which
  the message-passing form does not carry; it uses `θ = 0` throughout. See
  [Limitations](../project/limitations.md).

## Theta in relationship testing

Substituting the NRC II match probability for `P(G)` in the kinship formula
is wrong: the IBD-1 term would still draw B's free allele by the product
rule, so theta would apply to some terms and not others. The `...Theta`
functions instead treat every allele that is not identical by descent as a
draw from one sub-population, in sequence, with the Balding-Nichols sampling
formula (Balding & Nichols 1994):

```
P(next allele is x | n alleles drawn, m of them x) = (m θ + (1-θ) px) / (1 + (n-1) θ)
```

A pair draws four alleles under IBD 0, three under IBD 1 (B's IBD allele is a
copy of A's, not a draw) and two under IBD 2. For a trio index, `X` involves no
population draw at all; `Y` draws the random man's allele given the four
typed alleles of the mother and the alleged father, with denominator
`1 + 3θ`. The formulae are in
[Formulae and sources](../project/formulae.md#kinship-and-paternity-with-theta).

The effect is that a shared allele is less surprising in a structured
population, so an LR built on sharing a rare allele falls. Two `10/10` people,
`p10 = 0.1`:

| theta | parent/child LR | full-sibling LR |
| --- | --- | --- |
| 0.00 | 10 | 30.25 |
| 0.10 | 1.2 / 0.39 = 3.08 | 4.71 |

Theta does not always lower a number. A trio index rescued by mutation can
rise: with mother `10/10`, child `10/12` and alleged father `11/11`, the
random man's `12` is drawn given four alleles none of which is a `12`, so it
becomes rarer, `0.27 / 1.3` against `0.3`, and the index goes from `0.003` to
`0.0043` at `θ = 0.1`. A pair sharing no allele keeps its IBD-0 coefficient
as the LR at any theta. Report theta with the number either way.

## Theta and the 5/(2N) floor

Both make an estimate more conservative, for different reasons. The floor
handles sampling error in a finite reference table (see
[Reference frequency data](frequency-data.md)); theta handles population
structure. They compose: `alleleFrequency` applies the floor, and
`genotypeFrequency` applies theta to the floored values.
