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
`P(G)` is the Hardy-Weinberg frequency at `θ = 0`.

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
siblings?"). Both multiply over the loci both profiles and the table share, and
skip the rest without complaint, unlike `randomMatchProbability`. No shared
locus throws. `kinshipLRAt` throws if the alternative gives the observed pair
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

## Limits

- `θ = 0` throughout; see [Theta](../concepts/theta.md).
- Non-inbred pairs only. An inbred pedigree needs the nine condensed Jacquard
  coefficients.
- A third person (the mother in a trio, a second sibling, an untyped link) adds
  information the pairwise form discards. Use [Paternity](paternity.md) or
  [General pedigrees](pedigrees.md).
