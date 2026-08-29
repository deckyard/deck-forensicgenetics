# `forensicgenetics`

```jennifer
import "@jennifer/forensicgenetics/" as forensics;
```

The entry module: profiles, reference frequencies, match probability and
mixtures. Errors are `Error{kind: "forensicgenetics"}`.

## Constants

| Name | Type | Value | Meaning |
| --- | --- | --- | --- |
| `THETA_NONE` | `float` | `0.0` | no sub-population correction |
| `THETA_GENERAL` | `float` | `0.01` | NRC II, broad population |
| `THETA_ISOLATED` | `float` | `0.03` | NRC II, isolated population |
| `MIN_ALLELE_COUNT` | `float` | `5.0` | the numerator of the `5/(2N)` floor |
| `MAX_FREQUENCY_SUM` | `float` | `1.05` | the largest allele-frequency sum a locus may declare |

## Types

```jennifer
def struct LocusFrequencies {
    locus as string,
    n as int,                             # individuals typed
    alleles as map of string to float     # allele label -> published frequency
};

def struct FrequencyDb {
    name as string,
    loci as map of string to LocusFrequencies
};

def struct Genotype {
    locus as string,
    a as string,                          # canonically ordered pair
    b as string
};

def struct Profile {
    id as string,
    loci as map of string to Genotype
};

def struct Mixture {
    id as string,
    loci as map of string to list of string
};
```

## Allele labels

| Signature | Returns |
| --- | --- |
| `alleleValue(allele as string)` | `float` - the repeat count; throws if not numeric |
| `isNumericAllele(allele as string)` | `bool` - whether `alleleValue` would succeed |

## Frequency database

The generic, vendor-independent type. This module reads no published format;
a reader deck converts one into a `FrequencyDb` (see
[`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter)).

| Signature | Returns |
| --- | --- |
| `db(name as string)` | `FrequencyDb` - empty |
| `withLocus(fdb, locus as string, n as int, freqs as map of string to float)` | `FrequencyDb` - a copy with the locus added or replaced |
| `loci(fdb)` | `list of string` |
| `hasLocus(fdb, locus as string)` | `bool` |
| `sampleSize(fdb, locus as string)` | `int` - N |
| `minFrequency(fdb, locus as string)` | `float` - `5/(2N)`, or `0.0` when N is unknown |
| `alleleFrequency(fdb, locus as string, allele as string)` | `float` - floored at `5/(2N)`, capped at 1 |
| `combineFrequencies(dbs as list of FrequencyDb, name as string)` | `FrequencyDb` - several tables pooled into one, over their shared loci |

`withLocus` rejects an empty locus name, a negative `n`, a frequency outside
`(0, 1]`, or allele frequencies summing past `MAX_FREQUENCY_SUM`. A locus may
legitimately list only some of its alleles, so a sum *below* 1 is fine; a sum
meaningfully above it is a misread table, and accepting one yields genotype
"probabilities" greater than 1. Every lookup throws with a message naming the
missing locus.

`combineFrequencies` throws on an empty list, when the inputs share no locus,
or when a shared locus has no sample size anywhere; it keeps only the loci
every input covers, so check `loci` on the result to see what survived. Pooling
describes one merged population, and merging populations whose frequencies
differ biases the match probability *downwards* - the Wahlund effect. For an
uncertain population use `randomMatchProbabilityAcross` instead; see
[Limitations](../project/limitations.md).

## Genotypes and profiles

| Signature | Returns |
| --- | --- |
| `genotype(locus as string, a as string, b as string)` | `Genotype` - pair put in canonical order |
| `isHomozygous(g)` | `bool` |
| `alleles(g)` | `list of string` - one label for a homozygote, two otherwise |
| `profile(id as string, genotypes as list of Genotype)` | `Profile` - a repeated locus throws |
| `withGenotype(p, g)` | `Profile` - a copy, replacing any genotype at that locus |
| `profileLoci(p)` | `list of string` |
| `hasGenotype(p, locus as string)` | `bool` |
| `genotypeAt(p, locus as string)` | `Genotype` - throws if untyped there |
| `commonLoci(fdb, p)` | `list of string` - intersection, in profile order |
| `restrict(p, keep as list of string)` | `Profile` - only the named loci |

## Match probability

| Signature | Returns |
| --- | --- |
| `genotypeFrequency(fdb, g, theta as float)` | `float` - NRC II rec. 4.2, capped at 1 |
| `randomMatchProbability(fdb, p, theta as float)` | `float` - the product across loci |
| `log10RandomMatchProbability(fdb, p, theta as float)` | `float` - negative; cannot underflow |
| `matchLikelihoodRatio(fdb, p, theta as float)` | `float` - `1 / RMP` |
| `randomMatchProbabilityAcross(dbs as list of FrequencyDb, weights as list of float, p, theta as float)` | `float` - the RMP averaged over candidate populations, weighted by `weights` |
| `matchLikelihoodRatioAcross(dbs as list of FrequencyDb, weights as list of float, p, theta as float)` | `float` - `1 / randomMatchProbabilityAcross` |

`theta` must be in `[0, 1)`. All five panel functions throw on an empty
profile or on a locus the database does not cover - they do not silently drop
loci. See [Match probability](../casework/match-probability.md).

The `Across` pair is for a donor whose population is uncertain: `weights` are
the prior probabilities that the donor came from each population, one per
entry of `dbs`, and must be non-negative and sum to 1. Every table must cover
every locus in the profile, and `theta` is applied within each population.
They throw on an empty list, a length mismatch, a negative weight, or weights
that do not sum to 1; `matchLikelihoodRatioAcross` also throws if the weighted
probability is zero. Averaging the answer this way is not the same as pooling
the tables with `combineFrequencies`, which biases it - see
[Reference frequency data](../concepts/frequency-data.md).

## Mixtures

| Signature | Returns |
| --- | --- |
| `mixture(id as string, observed as map of string to list of string)` | `Mixture` - duplicates collapsed; an empty locus throws |
| `mixtureLoci(m)` | `list of string` |
| `inclusionProbability(fdb, locus as string, observed as list of string)` | `float` - `(Σ p)²`, capped at 1 |
| `exclusionProbability(fdb, locus as string, observed as list of string)` | `float` - `1 - PI` |
| `cpi(fdb, m)` | `float` - the product across loci |
| `cpe(fdb, m)` | `float` - `1 - CPI` |
| `isIncluded(p, m)` | `bool` - non-exclusion of one profile; throws if they share no locus |

None of these take a theta; see [Theta](../concepts/theta.md).

## Errors

`Error{kind: "forensicgenetics"}` for: an unknown locus, an unknown frequency
format, an empty or duplicated locus name, a frequency outside `(0, 1]`, a
locus whose frequencies sum past `MAX_FREQUENCY_SUM`, an empty allele label, a
theta outside `[0, 1)`, an empty profile or mixture locus, a profile and
mixture sharing no locus, a non-numeric allele passed to `alleleValue`, and a
locus with neither a frequency nor a sample size to floor with.
