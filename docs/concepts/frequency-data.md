# Reference frequency data

Every match probability in this deck is a statement about a **reference
population**, and that population arrives as a table of allele frequencies. The
table is not an implementation detail: change it and every number changes.

## What a FrequencyDb holds

```jennifer
def struct LocusFrequencies {
    locus as string,
    n as int,                            # individuals typed at this locus
    alleles as map of string to float    # allele label -> published frequency
};

def struct FrequencyDb {
    name as string,
    loci as map of string to LocusFrequencies
};
```

The `n` is the part people leave out, and it is the part that matters most
after the frequencies themselves.

## The 5/2N floor

A reference database of a thousand people has seen two thousand alleles at each
locus. Anything rarer than about one in two thousand it has probably not seen
at all - and an allele it has not seen has a published frequency of zero, which
would make the whole profile impossible.

NRC II's answer is a floor: no allele is ever scored rarer than

```
5 / (2N)
```

for that locus's own N. Five over the number of alleles observed. With N = 1000
that is 0.0025; with N = 100 it is 0.025, ten times larger, because a small
reference sample supports much weaker claims about rarity.

The deck applies it on read, so the stored table keeps its published numbers:

```jennifer
forensics.minFrequency($db, "TH01");              # 5 / (2N)
forensics.alleleFrequency($db, "TH01", "9.3");    # max(published, floor), capped at 1
```

An allele the table never observed comes back at the floor rather than at zero.
An observed allele rarer than the floor is raised to it. This is why
`LocusFrequencies` carries `n`: without it there is no floor to apply, and
`alleleFrequency` throws rather than return zero for an unobserved allele.

`MIN_ALLELE_COUNT` is the `5`, exported so it is visible rather than buried.

## Building a table by hand

```jennifer
def db as forensics.FrequencyDb init forensics.db("Example population");
$db = forensics.withLocus($db, "TH01", 1000,
    {"6": 0.25, "7": 0.15, "8": 0.30, "9.3": 0.30});
$db = forensics.withLocus($db, "D3S1358", 1000,
    {"15": 0.25, "16": 0.20, "17": 0.30, "18": 0.25});
```

`withLocus` returns a copy and validates as it goes: an empty locus name, a
negative `n`, or a frequency outside `(0, 1]` is an error. A frequency of
exactly zero is rejected - an allele that was not observed simply should not be
in the table, and the floor will handle it.

It also rejects a locus whose frequencies sum past `MAX_FREQUENCY_SUM` (1.05).
Listing only *some* of a locus's alleles is legitimate, so a sum below 1 is
fine; a sum well above it means the table was misread, and accepting one yields
genotype "probabilities" greater than 1.

## Getting a published table in

The core reads **no published format at all**. It defines `FrequencyDb`, and a
**reader deck** converts a particular source into one:

| Deck | Converts |
| --- | --- |
| [`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter) | STRidER-style XML |

```jennifer
import "@mplx/strider/" as strider;

def db as forensics.FrequencyDb init strider.load($xml, "STRidER AT");
```

Readers import this deck, never the reverse. That is what keeps the estimators
independent of any vendor's file format, and what lets a reader for another
source slot in without touching them.

Or build one directly, which is what the tests and examples do:

```jennifer
def db as forensics.FrequencyDb init forensics.withLocus(
    forensics.db("example"), "TH01", 1000,
    {"6": 0.25, "7": 0.15, "8": 0.30, "9.3": 0.30});
```

## Which dataset?

A match probability is not a property of a profile alone. It is a property of a
profile **and** a reference population, and STRidER publishes per-country,
per-continent and whole-database tables that give different answers. So the
dataset is named, not assumed - the label lands in `FrequencyDb.name` and
travels with the table.

### When the population is uncertain

Suppose the donor could be German or Austrian. Two operations exist, they give
different numbers, and only one answers that question.

**The law of total probability** - the right tool:

```jennifer
forensics.randomMatchProbabilityAcross([$de, $at], [0.5, 0.5], $suspect,
    forensics.THETA_GENERAL);
forensics.matchLikelihoodRatioAcross([$de, $at], [0.5, 0.5], $suspect,
    forensics.THETA_GENERAL);
```

It computes the match probability *within* each population and averages the
results, weighted by how likely the donor is to come from each. The weights are
priors and must sum to 1.

**Pooling** - one merged population:

```jennifer
forensics.combineFrequencies([$de, $at], "DE+AT");
```

This sums observed allele counts weighted by sample size, over the loci every
input carries. It describes a single merged population, which is the honest use
for it.

### Why the two differ

Pooling averages *allele frequencies*; the across form averages *genotype
probabilities*. The gap between them is the **Wahlund effect**: merging
populations whose frequencies differ leaves the pool short of heterozygotes
relative to Hardy-Weinberg. Since the match probability assumes Hardy-Weinberg,
a pooled table biases it downwards - overstating the evidence.

Two populations, one locus, equal sample sizes:

| | allele 6 | allele 9.3 |
| --- | --- | --- |
| population A | 0.8 | 0.2 |
| population B | 0.2 | 0.8 |
| pooled | 0.5 | 0.5 |

For a `6,6` homozygote with no theta correction:

| Method | Match probability |
| --- | --- |
| pooled table, `0.5²` | **0.25** |
| across, `0.5(0.8²) + 0.5(0.2²)` | **0.34** |

The pooled table overstates that match by a third, and on a full profile the
discrepancy compounds across every locus.

### The conservative alternative

Often you need not choose. Compute in each candidate population and report the
largest; if the conclusion holds across every plausible population, the choice
stops being a point of dispute.

## Inspecting a table

```jennifer
forensics.loci($db);                  # locus names, insertion order
forensics.hasLocus($db, "TH01");      # bool
forensics.sampleSize($db, "TH01");    # N
```

A locus the table does not carry is an error everywhere, with a message naming
the missing locus and the table - not the runtime's generic missing-key
complaint.

## Where real data comes from

The deck bundles no reference data. `testdata/frequencies.xml` contains
**synthetic** numbers, generated only to exercise the reader and the examples;
it is labelled as such in the file itself.

Real autosomal STR frequencies for European populations come from
[STRidER](https://strider.online/), the ENFSI reference database, whose
[formulae page](https://strider.online/formulae) is also the specification the
estimators here follow. Other populations and panels have their own published
tables. Whichever you use, the population it describes is the population your
match probability is about - say which one in anything you write down.
