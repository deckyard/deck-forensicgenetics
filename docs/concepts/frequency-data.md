# Reference frequency data

Every match probability is conditional on a reference population, and that
population arrives as a `FrequencyDb`: per locus, the allele frequencies and
`n`, the number of individuals typed (2N alleles). The struct fields are in the
[Cheatsheet](../reference/cheatsheet.md#forensicgeneticsfrequencydb).

## The 5/(2N) floor

A reference sample of N people has seen 2N alleles at a locus, so an allele
rarer than about 1/(2N) is often absent from the table. NRC II floors every
allele frequency at

```
5 / (2N)
```

for that locus's own N: 0.0025 at N = 1000, 0.025 at N = 100. The floor is
applied on read, so the stored table keeps its published numbers:

```jennifer
forensicgenetics.minFrequency($db, "TH01");              # 5 / (2N)
forensicgenetics.alleleFrequency($db, "TH01", "9.3");    # max(published, floor), capped at 1
```

An unlisted allele is scored at the floor; a listed allele below the floor is
raised to it. A locus with `n = 0` has no floor, and `alleleFrequency` throws
for an allele it does not list. `MIN_ALLELE_COUNT` is the 5.

## Building a table

```jennifer
def db as forensicgenetics.FrequencyDb init forensicgenetics.db("Example population");
$db = forensicgenetics.withLocus($db, "TH01", 1000,
    {"6": 0.25, "7": 0.15, "8": 0.30, "9.3": 0.30});
$db = forensicgenetics.withLocus($db, "D3S1358", 1000,
    {"15": 0.25, "16": 0.20, "17": 0.30, "18": 0.25});
```

`withLocus` returns a copy and throws on an empty locus name, a negative `n`,
a frequency outside `(0, 1]` (leave an unobserved allele out; the floor covers
it), or frequencies summing past `MAX_FREQUENCY_SUM` (1.05). A sum below 1 is
accepted, since a table may list only part of a locus.

Each call copies the whole database; see
[Limitations](../project/limitations.md#building-a-database-or-profile-one-entry-at-a-time-is-quadratic)
for the cost at hundreds of loci.

## Published tables

The core reads no file format. A reader deck imports the core and produces a
`FrequencyDb`:

| Deck | Converts |
| --- | --- |
| [`@deckyard/strider`](https://github.com/deckyard/deck-strider-frequency-converter) | STRidER-style XML |

```jennifer
import "@deckyard/strider/" as strider;

def db as forensicgenetics.FrequencyDb init strider.load($xml, "STRidER AT");
```

The deck bundles no reference data. The table in `examples/lib/exampledb.j` is
synthetic. Real autosomal STR frequencies for European populations come from
[STRidER](https://strider.online/), the ENFSI reference database, whose
[formulae page](https://strider.online/formulae) is the specification these
estimators follow.

## Which dataset

STRidER publishes per-country, per-continent and whole-database tables, and
they give different answers. Name the dataset when you build the table; the
label travels as `FrequencyDb.name`.

```jennifer
forensicgenetics.loci($db);                  # locus names, insertion order
forensicgenetics.hasLocus($db, "TH01");      # bool
forensicgenetics.sampleSize($db, "TH01");    # N
```

A locus the table lacks throws everywhere, with a message naming the locus and
the table.

### When the population is uncertain

For "the donor was German or Austrian", use the law of total probability:

```jennifer
forensicgenetics.randomMatchProbabilityAcross([$de, $at], [0.5, 0.5], $suspect,
    forensicgenetics.THETA_GENERAL);
forensicgenetics.matchLikelihoodRatioAcross([$de, $at], [0.5, 0.5], $suspect,
    forensicgenetics.THETA_GENERAL);
```

It computes the match probability within each population and averages the
results, weighted by the prior that the donor came from each. The weights must
be non-negative and sum to 1 (within 0.001).

`combineFrequencies` pools tables instead:

```jennifer
forensicgenetics.combineFrequencies([$de, $at], "DE+AT");
```

It sums observed allele counts (`p * 2N`) over the loci every input carries, and
describes one merged population. Pooling averages allele frequencies where the
`Across` form averages genotype probabilities. A merged population whose parts
differ has fewer heterozygotes than Hardy-Weinberg predicts (the Wahlund
effect), so the RMP from a pooled table is biased downwards and overstates the
evidence.

Two populations, one locus, equal sample sizes:

| | allele 6 | allele 9.3 |
| --- | --- | --- |
| population A | 0.8 | 0.2 |
| population B | 0.2 | 0.8 |
| pooled | 0.5 | 0.5 |

For a `6,6` homozygote at theta 0:

| Method | Match probability |
| --- | --- |
| pooled table, `0.5²` | 0.25 |
| across, `0.5(0.8²) + 0.5(0.2²)` | 0.34 |

The pooled answer overstates this match by a third at one locus, and the gap
compounds across a profile. Use `combineFrequencies` only when the merged
group is itself the reference population.

A conservative alternative is to compute in each candidate population and
report the largest RMP.
