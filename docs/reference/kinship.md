# `kinship`

```jennifer
import "@jennifer/forensicgenetics/kinship.j" as kinship;
```

Relationship likelihood ratios and paternity. Builds on `forensicgenetics`, so its
signatures take `forensics.Profile`, `forensics.Genotype` and
`forensics.FrequencyDb`. Errors are `Error{kind: "kinship"}`.

Uses the plain product rule throughout - `θ = 0`. See
[Theta](../concepts/theta.md).

## Relationship constants

| Name | k0 | k1 | k2 |
| --- | --- | --- | --- |
| `UNRELATED` | 1 | 0 | 0 |
| `PARENT_CHILD` | 0 | 1 | 0 |
| `FULL_SIBLINGS` | 0.25 | 0.5 | 0.25 |
| `HALF_SIBLINGS` | 0.5 | 0.5 | 0 |
| `GRANDPARENT` | 0.5 | 0.5 | 0 |
| `AVUNCULAR` | 0.5 | 0.5 | 0 |
| `FIRST_COUSINS` | 0.75 | 0.25 | 0 |
| `IDENTICAL` | 0 | 0 | 1 |

`HALF_SIBLINGS`, `GRANDPARENT` and `AVUNCULAR` are IBD-identical and give the
same likelihood ratio; the three names exist so a hypothesis can be stated
honestly.

## Mutation constants

| Name | rate | decay |
| --- | --- | --- |
| `NO_MUTATION` | 0.0 | 0.0 |
| `DEFAULT_MUTATION` | 0.002 | 0.1 |

`DEFAULT_MUTATION` is a placeholder in the right order of magnitude, not a
published per-locus figure. Override with `withLocusRate`.

## Types

```jennifer
def struct Relationship {
    name as string,
    k0 as float,      # P(share no allele IBD)
    k1 as float,      # P(share exactly one)
    k2 as float       # P(share both)
};

def struct MutationModel {
    rate as float,                        # per-locus, per-meiosis
    decay as float,                       # geometric decay per extra step
    perLocus as map of string to float    # locus -> rate, overriding `rate`
};
```

## Relationships

| Signature | Returns |
| --- | --- |
| `relationship(name as string, k0 as float, k1 as float, k2 as float)` | `Relationship` - coefficients must be in `[0, 1]` and sum to 1 |

## Transmission and mutation

| Signature | Returns |
| --- | --- |
| `transmissionProbability(g as forensics.Genotype, allele as string)` | `float` - 1.0, 0.5 or 0.0 |
| `mutationModel(rate as float, decay as float)` | `MutationModel` - both in `[0, 1)` |
| `withLocusRate(model, locus as string, rate as float)` | `MutationModel` - a copy with one locus overridden |
| `mutationProbability(model, locus as string, source as string, target as string)` | `float` |

`mutationProbability` is `1 - rate` when `source == target`, and
`rate · (1 - decay) · decay^(k-1) / 2` for a `k`-step change. A microvariant
difference below one repeat counts as one step. A changed allele with no
numeric repeat count throws.

## Joint probabilities and kinship LRs

| Signature | Returns |
| --- | --- |
| `jointGenotypeProbability(db, ga, gb, rel)` | `float` - `k0·P(Ga)P(Gb) + k1·T1 + k2·[Ga==Gb]·P(Ga)` |
| `kinshipLRAt(db, ga, gb, rel, alt)` | `float` - one locus, `rel` over `alt` |
| `pairLoci(db, a, b)` | `list of string` - the loci both profiles and the database share |
| `kinshipLR(a, b, rel, db)` | `float` - `rel` against `UNRELATED`, across loci |
| `kinshipLRAgainst(a, b, rel, alt, db)` | `float` - `rel` against `alt`, across loci |

`jointGenotypeProbability` throws if the two genotypes are at different loci;
`kinshipLRAt` throws if the alternative gives the observed pair probability
zero. `kinshipLR` and `kinshipLRAgainst` throw if the pair shares no covered
locus.

## Paternity

| Signature | Returns |
| --- | --- |
| `trioLoci(db, mother, child, father)` | `list of string` - loci all three and the database share |
| `paternityIndexAt(db, mother, child, father, model)` | `float` - `X / Y` at one locus |
| `paternityIndex(mother, child, allegedFather, db)` | `float` - CPI under `DEFAULT_MUTATION` |
| `paternityIndexWith(mother, child, allegedFather, db, model)` | `float` - CPI under an explicit model |
| `inconsistentLoci(mother, child, allegedFather, db)` | `list of string` - loci with a zero strict PI |
| `probabilityOfPaternity(cpi as float)` | `float` - `CPI / (CPI + 1)` |
| `probabilityOfPaternityPrior(cpi as float, prior as float)` | `float` - prior in `(0, 1)` |

`paternityIndexAt` and everything built on it throw on a **maternal
exclusion** - a child carrying neither of the mother's alleles - rather than
returning zero. The mutation model applies to the alleged father's transmission
only. See [Paternity](../casework/paternity.md).

## Errors

`Error{kind: "kinship"}` for: IBD coefficients out of range or not summing to
1; a rate or decay outside `[0, 1)`; mismatched loci in a genotype pair; a zero
denominator in a likelihood ratio; no shared locus for a pair or a trio; a
maternal exclusion; a negative CPI; a prior outside `(0, 1)`; and a mutation
step between alleles with no repeat count.
