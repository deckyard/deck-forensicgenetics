# Population statistics

A published frequency table comes with a forensic-efficiency table: per locus,
how informative it is for identification and for parentage, and whether the
sample it was estimated from fits Hardy-Weinberg proportions. `popstats.j`
computes both.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/popstats.j" as popstats;
```

## The table for a frequency database

```jennifer
def table as list of popstats.LocusStats init popstats.statsTable($db);
def all as popstats.CombinedStats init popstats.combinedStats($table);
```

One `LocusStats` per locus:

| Field | Meaning |
| --- | --- |
| `expectedHeterozygosity` | chance two random alleles differ, `1 - Σ p²` |
| `matchProbability` | chance two unrelated people share a genotype |
| `powerOfDiscrimination` | `1 - matchProbability` |
| `pic` | polymorphism information content (Botstein et al. 1980) |
| `typicalPaternityIndex` | `1 / (2 Σ p²)` |
| `powerOfExclusionTrio` | share of random men excluded, mother typed |
| `powerOfExclusionDuo` | share of random men excluded, mother untyped |

All assume Hardy-Weinberg proportions. Each has its own function, e.g.
`popstats.powerOfDiscrimination($db, "TH01")`.

`matchProbability` averages over all genotypes. It is not the random match
probability of a particular profile; that is
`forensicgenetics.randomMatchProbability`.

### Which power of exclusion

Both are the exact Hardy-Weinberg expressions of Jamieson & Taylor (1997).
Use the trio value when the mother is typed, the duo value when only child and
alleged father are. A duo excludes a man only when he shares no allele with
the child, so it is always lower. `observedPowerOfExclusion` is the Brenner &
Morris (1990) approximation from observed heterozygosity, which PowerStats
reports; it is not interchangeable with either.

### Combined values

`combinedStats` multiplies over loci assuming independence. A long panel's
match probability underflows; read `log10MatchProbability`. The same holds for
the non-exclusion products (`log10NonExclusionTrio`, `log10NonExclusionDuo`).
For per-locus powers from elsewhere, `combinedPower` and
`log10CombinedComplement` do the same arithmetic on a `list of float`.

### Small samples

`expectedHeterozygosity` is the plug-in estimate. `unbiasedHeterozygosity`
applies Nei's (1978) `2N / (2N - 1)` correction with the locus's `n`. Report it
when the table states a sample size.

## Statistics of a genotype sample

A sample is a `list of forensicgenetics.Genotype` at one locus:

```jennifer
def ho as float init popstats.observedHeterozygosity($sample);
def freqs as map of string to float init popstats.estimateFrequencies($sample);
def own as forensicgenetics.FrequencyDb init popstats.withEstimatedLocus(
    forensicgenetics.db("my sample"), $sample);
def he as float init popstats.unbiasedHeterozygosity($own, "TH01");
```

`withEstimatedLocus` stores the counted frequencies with `n` = sample size, so
every frequency statistic and the core's match probabilities apply to the
sample. Compare `ho` with `he`: a deficit of heterozygotes suggests null
alleles, allele dropout, or a structured sample.

## Hardy-Weinberg tests

```jennifer
def p2 as float init popstats.hweExactTest($sample);
def pk as float init popstats.hwePermutationTest($sample, popstats.DEFAULT_PERMUTATIONS, 1);
```

`hweExactTest` is the exact test of Wigginton et al. (2005) and accepts at most
two alleles. STR loci have more; `hwePermutationTest` shuffles the sample's
alleles, re-pairs them, and counts tables no more likely than the observed one
(Guo & Thompson 1992). The p-value is `(r + 1) / (B + 1)`, with standard error
about `sqrt(p (1 - p) / B)`. Report the permutation count and seed; the same
seed reproduces the value.

Testing many loci at 0.05 produces false rejections; correct for the number of
tests (e.g. Bonferroni) before calling a locus out of equilibrium.

See `examples/popstats.j` for a full table and a sample.

See also [Population structure and linkage](population-structure.md) for
F-statistics across populations and tests of independence between loci.
