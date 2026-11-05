# Kinship

A kinship likelihood ratio weighs allele sharing by descent against sharing by
chance, for a pair of people.

## IBD coefficients

A non-inbred pair's relationship is `(k0, k1, k2)`: the probabilities of
sharing 0, 1 or 2 alleles identical by descent at a locus.

| Relationship | k0 | k1 | k2 |
| --- | --- | --- | --- |
| `UNRELATED` | 1 | 0 | 0 |
| `PARENT_CHILD` | 0 | 1 | 0 |
| `FULL_SIBLINGS` | 1/4 | 1/2 | 1/4 |
| `HALF_SIBLINGS` | 1/2 | 1/2 | 0 |
| `GRANDPARENT` | 1/2 | 1/2 | 0 |
| `AVUNCULAR` | 1/2 | 1/2 | 0 |
| `FIRST_COUSINS` | 3/4 | 1/4 | 0 |
| `IDENTICAL` | 0 | 0 | 1 |

Any other relationship:

```jennifer
def dfc as kinship.Relationship init kinship.relationship(
    "double first cousins", 0.5625, 0.375, 0.0625);
```

Each coefficient must be in `[0, 1]` and the three must sum to 1 (within
1e-9).

`HALF_SIBLINGS`, `GRANDPARENT` and `AVUNCULAR` have the same coefficients and
give exactly the same LR: autosomal STRs cannot separate them. The three names
exist so a report can state the hypothesis. Separating them needs
[lineage markers](lineage.md), a denser panel, or more relatives in a
[pedigree](pedigrees.md).

## Joint probability

At one locus:

```
P(Ga, Gb) = k0 · P(Ga) P(Gb) + k1 · T1(Ga, Gb) + k2 · [Ga == Gb] · P(Ga)
```

`T1` is the probability with one of B's alleles a copy of one of A's and the
other drawn from the population; its definition is in
[Formulae and sources](../project/formulae.md#kinship-joint-genotype-probability).
`P(G)` is the Hardy-Weinberg frequency at `θ = 0`; for theta see
[With theta](#with-theta).

```jennifer
kinship.jointGenotypeProbability($db, $ga, $gb, kinship.FULL_SIBLINGS);
```

`T1` is symmetric, so the LR does not depend on which person is A; the test
suite asserts this.

## Likelihood ratios

```jennifer
def lr as float init kinship.kinshipLR($a, $b, kinship.FULL_SIBLINGS, $db);
def lr2 as float init kinship.kinshipLRAgainst($a, $b,
    kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, $db);
kinship.kinshipLRAt($db, $ga, $gb, $rel, $alt);   # one locus
kinship.pairLoci($db, $a, $b);                    # the loci used
```

`kinshipLR` tests against `UNRELATED`; `kinshipLRAgainst` against any
alternative, which is the form most real questions take ("full or half
siblings?"). Both multiply over the loci both profiles carry. A locus only one
profile carries is skipped, since nothing can be computed there; a locus both
carry but the table does not cover throws, as in `randomMatchProbability`,
because dropping it could hide an exclusion. Narrow the profiles with
`forensicgenetics.restrict` to leave a locus out on purpose. No shared locus
throws. `kinshipLRAt` throws if the alternative gives the observed pair
probability zero.

## Closed forms

At one locus, with allele frequency `p`:

| Case | LR against unrelated |
| --- | --- |
| parent/child, `aa` and `ab` | `1 / (2 pa)` |
| full sibs, both `aa` | `(1 + 1/p)² / 4` |
| full sibs, both `ab` | `1/4 + (1 + pa + pb) / (8 pa pb)` |
| identical, genotype `g` | `1 / P(g)` |
| parent/child sharing no allele | `0` |

At `p = 0.1`, two `aa` people give a sibling LR of `(1 + 10)²/4 = 30.25`. The
test suite checks these closed forms.

A parent and child must share an allele at every locus, so one locus with none
shared makes the parent/child LR zero. Full siblings need not share, so the
same locus only counts against them.

## With theta

Each kinship function has a `...Theta` variant taking the co-ancestry
coefficient last:

```jennifer
def lr as float init kinship.kinshipLRTheta($a, $b, kinship.FULL_SIBLINGS, $db, 0.01);
def lr2 as float init kinship.kinshipLRAgainstTheta($a, $b,
    kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, $db, 0.01);
kinship.kinshipLRAtTheta($db, $ga, $gb, $rel, $alt, 0.01);
kinship.jointGenotypeProbabilityTheta($db, $ga, $gb, $rel, 0.01);
```

Every allele not identical by descent is drawn from one sub-population by
the Balding-Nichols sampling formula: four draws under IBD 0, three under
IBD 1, two under IBD 2 (formulae in
[Formulae and sources](../project/formulae.md#kinship-and-paternity-with-theta)).
This is not the theta-0 formula with the NRC II match probability put in for
`P(G)`, which would leave the IBD-1 term on the product rule. At
`theta = 0` each variant returns exactly the original's value. Theta must be
in `[0, 1)`.

Closed forms with theta, `u = θ + (1-θ) pa`, `v = θ + (1-θ) pb`; the test
suite checks each at `θ = 0.1`:

| Case | LR against unrelated |
| --- | --- |
| parent/child, both `aa` | `(1 + 2θ) / (3θ + (1-θ) pa)` |
| parent/child, both `ab` | `(1 + 2θ)(u + v) / (4 u v)` |
| parent/child, `ab` and `ac` | `(1 + 2θ) / (4 u)` |
| identical, genotype `g` | `1 / genotypeFrequency(g, θ)` |
| full sibs | `1/4 + 1/2 · LR(parent/child) + 1/4 · LR(identical)` |

With `p = 0.1` and `θ = 0.1`, two `aa` people give a parent/child LR of
`1.2 / 0.39 = 3.08` (10 at `θ = 0`) and a full-sibling LR of 4.71 (30.25).

## Limits

- The original functions use `θ = 0`; the `...Theta` variants apply the
  correction. See [Theta](../concepts/theta.md).
- Non-inbred pairs only. An inbred pedigree needs the nine condensed Jacquard
  coefficients.
- A third person (the mother in a trio, a second sibling, an untyped link) adds
  information the pairwise form discards. Use [Paternity](paternity.md) or
  [General pedigrees](pedigrees.md).
