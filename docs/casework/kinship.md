# Kinship

Two people share alleles for two reasons: coincidence, and common ancestry. A
kinship likelihood ratio separates them.

## IBD coefficients

For a pair of non-inbred people, an entire relationship reduces to three
numbers - the probability that at a given locus they share 0, 1 or 2 alleles
**identical by descent**:

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

A parent and child always share exactly one allele by descent - hence
`(0, 1, 0)`. Full siblings inherit each parental allele independently, giving
the quarter/half/quarter split.

Other relationships:

```jennifer
def dfc as kinship.Relationship init kinship.relationship(
    "double first cousins", 0.5625, 0.375, 0.0625);
```

The three must be in `[0, 1]` and sum to 1.

### Three relationships that look identical

`HALF_SIBLINGS`, `GRANDPARENT` and `AVUNCULAR` have the same coefficients, and
therefore produce **exactly the same likelihood ratio**. Autosomal STRs cannot
tell a half-brother from an uncle from a grandfather. That is a property of the
markers, not a shortcut in this code - the three constants exist so the
hypothesis can be named honestly in a report, not because the numbers differ.

Separating them needs [lineage markers](lineage.md), a denser panel, or the
extra structure a [pedigree](pedigrees.md) supplies.

## The joint probability

At one locus, for genotypes `Ga` and `Gb`:

```
P(Ga, Gb) = k0 · P(Ga) P(Gb)
          + k1 · T1(Ga, Gb)
          + k2 · [Ga == Gb] · P(Ga)
```

The first term is independence - no shared allele. The third is identity -
both shared, so the genotypes must be equal. The middle term `T1` is the
interesting one: one of B's alleles is a copy of one of A's, and the other is
drawn from the population.

```jennifer
kinship.jointGenotypeProbability($db, $ga, $gb, kinship.FULL_SIBLINGS);
```

`T1` is symmetric in its arguments, so the LR does not depend on which person
you call A - and the test suite asserts that rather than assuming it.

## Likelihood ratios

```jennifer
def lr as float init kinship.kinshipLR($a, $b, kinship.FULL_SIBLINGS, $db);
```

tests the relationship against being unrelated, multiplying over every locus
the two profiles and the database share. Above 1 favours the relationship;
below 1 favours unrelatedness.

To test one relationship against another rather than against unrelatedness:

```jennifer
def lr as float init kinship.kinshipLRAgainst($a, $b,
    kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, $db);
```

That is the form most real questions take. "Are they full siblings or half
siblings" is a genuine question; "are they full siblings or complete strangers"
usually is not.

Per locus:

```jennifer
kinship.kinshipLRAt($db, $ga, $gb, $rel, $alt);
```

Which loci contributed:

```jennifer
kinship.pairLoci($db, $a, $b);
```

Unlike `randomMatchProbability`, kinship uses the intersection without
complaint - two people's panels genuinely differ, and there is no other
sensible reading. Print `pairLoci` if the count matters.

## Numbers worth recognising

With allele frequency `p`, at one locus:

| Case | LR against unrelated |
| --- | --- |
| Parent/child, `aa` and `ab` | `1 / (2 pa)` |
| Full sibs, both `aa` | `(1 + 1/p)² / 4` |
| Full sibs, both `ab` | `1/4 + (1 + pa + pb) / (8 pa pb)` |
| Identical, genotype `g` | `1 / P(g)` |
| Parent/child sharing no allele | `0` |

At `p = 0.1`, two people both homozygous `10/10` give a sibling LR of
`(1 + 10)²/4 = 30.25`. The test suite checks exactly these closed forms.

Note the last row. Parent and child *must* share an allele at every locus, so a
locus where they share none drives the LR to zero - a single-locus exclusion.
Full siblings need not, so the same observation only counts against them.

## Theta is not applied

Kinship uses the plain product rule, `θ = 0`. The sub-population correction
does not carry into the `T1` term as a simple substitution, and applying it to
some terms and not others would produce a number that answers no question. See
[Theta](../concepts/theta.md) and [Limitations](../project/limitations.md).

## Inbreeding

The three IBD coefficients describe **non-inbred** pairs. An inbred pedigree
needs the nine condensed Jacquard coefficients, which this pairwise form does
not carry. A consanguineous pedigree is also rejected by
[`pedigree.j`](pedigrees.md) rather than mis-summed.

## When pairwise is not enough

As soon as a third person is involved - a mother in a paternity trio, a second
sibling, an untyped link between two typed people - the pairwise form throws
away information. Those cases go to [Paternity](paternity.md) and
[General pedigrees](pedigrees.md), both of which condition on everyone at once.

## Related

- [Paternity](paternity.md) - the trio special case
- [General pedigrees](pedigrees.md) - three or more relatives
- [`kinship` reference](../reference/kinship.md)
