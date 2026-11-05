# Population structure and linkage

Two checks on the samples behind a frequency database: how far its
sub-populations differ (the θ that the NRC II correction stands in for), and
whether its loci are independent, which every product over loci assumes.
Both are in `popstats.j`.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/popstats.j" as popstats;
```

## F-statistics

Weir & Cockerham (1984) estimate three correlations from genotype samples of
several populations:

| Estimate | Function | Measures |
| --- | --- | --- |
| θ (F_ST) | `popstats.fst` | coancestry of two alleles within a population, relative to the whole |
| F (F_IT) | `popstats.fit` | correlation of the two alleles in an individual, relative to the whole |
| f (F_IS) | `popstats.fis` | the same within populations: inbreeding, or a heterozygote deficit |

They are ratios of three variance components, `a` (between populations),
`b` (between individuals within populations) and `c` (within individuals).
`fStatistics` computes the components of one locus; the ratio functions take
them:

```jennifer
def stats as popstats.FStatistics init popstats.fStatistics([$sampleA, $sampleB, $sampleC]);
def theta as float init popstats.fst($stats);
```

Each sample is a `list of forensicgenetics.Genotype` at one locus, one per
person. From profiles grouped by population, `fStatisticsTable` does every
locus at once, leaving out people untyped at a locus:

```jennifer
def table as list of popstats.FStatistics init popstats.fStatisticsTable($populations, ["D3S1358", "TH01", "FGA"]);
def all as popstats.FStatistics init popstats.combineFStatistics($table);
def thetaAll as float init popstats.fst($all);
```

The multilocus value is the ratio of the components summed over loci, Weir &
Cockerham's estimator, not the mean of per-locus ratios.

### Reading the values

- The components are unbiased, the ratio nearly so. When populations do not
  differ, θ scatters around 0 and is often slightly negative: for identical
  samples of `n` people it is exactly `-(1 - Ho / (2 He)) / (n - 1)`, about
  `-1 / (2 (n - 1))` in Hardy-Weinberg proportions. Read a small negative
  value as 0, not as an error.
- A θ estimated between sub-populations of one database is one input to the
  choice of `forensicgenetics.THETA_GENERAL` / `THETA_ISOLATED`, not a
  replacement: NRC II's θ is a coancestry allowance for the unknown
  sub-population of a suspect, and an estimate from a few groups carries
  sampling error across loci.
- `fis` well above 0 at one locus only suggests a null allele or allele
  drop-out in the typing; across all loci, inbreeding or a Wahlund effect
  from pooling unlike groups.

### One population

`inbreedingCoefficient(genotypes)` is f for a single sample, the same
estimator with one population: `fis(fStatistics([genotypes]))`. With one
population `fst` and `fit` throw.

## Linkage disequilibrium

`ldPermutationTest` tests whether genotypes at one locus are independent of
genotypes at another in the same people. Gametic phase is not needed and
Hardy-Weinberg proportions are not assumed, which suits unphased STR data.

```jennifer
def p as float init popstats.ldPermutationTest($th01, $fga, popstats.DEFAULT_PERMUTATIONS, 1);
```

The two lists must hold the same people in the same order. The p-value is a
seeded Monte Carlo estimate, like `hwePermutationTest`: the same seed gives
the same answer, and its standard error is about `sqrt(p (1 - p) / B)`.

For a whole database, `ldPermutationTests` runs every pair of loci:

```jennifer
def pairs as list of popstats.LdTest init popstats.ldPermutationTests($profiles, 2000, 7);
```

### Many pairs, many tests

A 20-locus set has 190 pairs. Under independence about 5% of them fall below
0.05 by chance, so a handful of small p-values is expected. Judge them
against a correction: Bonferroni rejects at `0.05 / 190`. A pair that stays
significant, or a locus that is significant against many others, points to
physical linkage, a typing artefact, or a sample pooled from different
populations.

Cost is about pairs × permutations × people: budget the permutation count on
a large set, and raise it only for pairs near the threshold.

See [Formulae](../project/formulae.md#f-statistics) and the
[popstats reference](../reference/popstats.md).
