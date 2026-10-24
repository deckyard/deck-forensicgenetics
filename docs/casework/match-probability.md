# Match probability

A stain matches a suspect's single-source profile. The random match
probability (RMP) is the probability that an unrelated person from the
reference population carries the same profile.

## One locus

```jennifer
def f as float init forensicgenetics.genotypeFrequency($db, $g, forensicgenetics.THETA_GENERAL);
```

NRC II recommendation 4.2 with the theta correction, on allele frequencies
floored at `5/(2N)`, capped at 1. See [Theta](../concepts/theta.md) and
[Reference frequency data](../concepts/frequency-data.md).

## Across a panel

```jennifer
def rmp as float init forensicgenetics.randomMatchProbability($db, $p, $theta);
def lr as float init forensicgenetics.matchLikelihoodRatio($db, $p, $theta);
```

`randomMatchProbability` is the product of `genotypeFrequency` over every locus
of the profile. `matchLikelihoodRatio` is `1 / RMP`: "the stain came from this
person" against "it came from an unrelated person from the reference
population". For a relative as the alternative, the RMP is the wrong number;
see [Kinship](kinship.md).

## Long panels

Twenty loci give an RMP near 10⁻²⁵, which a float holds. Longer panels can
underflow. `log10RandomMatchProbability` sums logarithms instead:

```jennifer
def logRmp as float init forensicgenetics.log10RandomMatchProbability($db, $p, $theta);
io.printf("about 1 in 10^%f|prec=1\n", -$logRmp);
```

## Missing loci are refused, not skipped

`randomMatchProbability` throws on a profile locus the table does not cover:

```
frequency database 'Example' has no locus 'SE33'
```

Dropping the locus would raise the RMP without any sign in the output. To
score the shared loci only, narrow the profile and report which loci remain:

```jennifer
def narrowed as forensicgenetics.Profile init forensicgenetics.restrict($p,
    forensicgenetics.commonLoci($db, $p));
def rmp as float init forensicgenetics.randomMatchProbability($db, $narrowed, $theta);
def used as list of string init forensicgenetics.profileLoci($narrowed);
```

An empty profile throws as well.

## Worked example

`examples/matchprobability.j`, synthetic table, `θ = 0.01`:

```
locus       genotype    frequency
D3S1358      15/16        0.033724
vWA          17/17        0.005966
FGA          21/24        0.028526
D8S1179      12/14        0.002559
D21S11       29/31.2      0.014348
D18S51       14/17        0.009885
TH01         6/9.3        0.030870
D16S539      11/12        0.025092

log10 RMP = -14.792   ->  about 1 in 6.2 x 10^14
```

No single locus is decisive; the answer is the product of eight ordinary
frequencies.

## What to report with the number

- the reference population: `$db.name` and the table's source
- theta: 0.01 and 0.03 differ by about 1.7 orders of magnitude on this profile
- the loci used, especially after `restrict`
- that the alternative is an unrelated person, and a kinship LR if relatives
  are possible

See also [Mixtures](mixtures.md) for more than one contributor.
