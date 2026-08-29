# Match probability

The single-source question: a stain matches a suspect's profile, so how
surprising would that match be if the stain had come from someone else?

## One locus

```jennifer
def f as float init forensics.genotypeFrequency($db, $g, forensics.THETA_GENERAL);
```

The NRC II recommendation 4.2 genotype frequency - `pa²` and `2 pa pb`
corrected for co-ancestry. The formulae and the choice of theta are in
[Theta and the sub-population correction](../concepts/theta.md); the `5/(2N)`
floor applied to `pa` and `pb` on the way in is in
[Reference frequency data](../concepts/frequency-data.md).

## Across a panel

```jennifer
def rmp as float init forensics.randomMatchProbability($db, $p, $theta);
def lr as float init forensics.matchLikelihoodRatio($db, $p, $theta);
```

`randomMatchProbability` multiplies `genotypeFrequency` over every locus of the
profile. `matchLikelihoodRatio` is `1 / RMP` - the weight of evidence for *the
stain came from this person* against *it came from an unrelated person from the
reference population*.

Note what the denominator hypothesis says: **unrelated**. A relative is not an
unrelated person, and for a brother the right number is not the RMP at all -
see [Kinship](kinship.md).

## Long panels: use the logarithm

Twenty loci put an RMP somewhere around 10⁻²⁵, which a float carries fine.
Fifty or a hundred loci do not. `log10RandomMatchProbability` sums logarithms
instead of multiplying probabilities:

```jennifer
def logRmp as float init forensics.log10RandomMatchProbability($db, $p, $theta);
io.printf("about 1 in 10^%f|prec=1\n", -$logRmp);
```

It is also the form a report usually quotes - "one in 10¹⁸" rather than a
number with eighteen leading zeros.

## Missing loci are refused, not skipped

`randomMatchProbability` throws if the profile carries a locus the frequency
table does not cover:

```
frequency database 'Example' has no locus 'SE33'
```

This is deliberate. Dropping the locus quietly would raise the RMP - weaker
evidence - without anything in the output saying a locus had been discarded.
That is the failure mode where a wrong number looks exactly like a right one.

When narrowing the panel really is what you want, say so:

```jennifer
def narrowed as forensics.Profile init forensics.restrict($p,
    forensics.commonLoci($db, $p));
def rmp as float init forensics.randomMatchProbability($db, $narrowed, $theta);
```

Now the loci that contributed are `forensics.profileLoci($narrowed)`, and you
can print them.

An empty profile is an error too.

## A worked example

From `examples/matchprobability.j`, against the bundled synthetic table at
`θ = 0.01`:

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

Two things to read off it. The homozygote `17/17` at vWA contributes far more
than any heterozygote - squaring a frequency is what drives a match
probability. And the whole answer is a product of eight ordinary-looking
numbers; no single locus is decisive, which is exactly why a panel is used.

## Reporting checklist

A match probability is meaningless without the assumptions that produced it.
Whatever you print, print these beside it:

- **which reference population** - `$db.name`, and where the table came from
- **which theta** - 0.01 and 0.03 differ by two orders of magnitude here
- **which loci contributed** - `forensics.profileLoci($p)`, especially if the
  profile was narrowed
- **that the alternative is an unrelated person** - and, if relatives are in
  the frame, a kinship LR as well

## Related

- [Mixtures](mixtures.md) - when the stain has more than one contributor
- [Kinship](kinship.md) - when the alternative suspect is a relative
- [Theta](../concepts/theta.md) - the correction and its cost
- [`forensicgenetics` reference](../reference/forensicgenetics.md) - every signature
