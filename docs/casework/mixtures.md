# Mixtures

A mixed stain has no genotype. All that can be read off it is the set of
alleles present at each locus, with no reliable way to say which contributor
brought which - and, often, no reliable count of contributors at all.

The deck answers the modest question that needs no such assumptions: **what
fraction of the population could not be excluded as a contributor?**

## Building a mixture

```jennifer
def stain as forensics.Mixture init forensics.mixture("stain", {
    "D3S1358": ["15", "16", "17"],
    "vWA": ["16", "17"],
    "FGA": ["21", "22", "24"],
    "TH01": ["6", "7", "9.3"]
});
```

Repeated alleles are collapsed; a locus with an empty list is an error.

## Probability of inclusion

At one locus, a person cannot be excluded if *both* their alleles are among
those the stain shows. Under Hardy-Weinberg the chance of that for a random
person is the square of the summed frequencies:

```
PI = (Σ observed allele frequencies)²
PE = 1 - PI
```

```jennifer
forensics.inclusionProbability($db, "D3S1358", ["15", "16", "17"]);
forensics.exclusionProbability($db, "D3S1358", ["15", "16", "17"]);
```

The sum is capped at 1.0, which the `5/(2N)` floor can otherwise push it past
at a locus showing many alleles.

Combined across loci:

```jennifer
def cpi as float init forensics.cpi($db, $stain);
def cpe as float init forensics.cpe($db, $stain);
```

`CPI` is the product of the per-locus `PI`; `CPE = 1 - CPI` is the number a
mixture report usually quotes.

## Is a given profile included?

```jennifer
forensics.isIncluded($suspect, $stain);   # bool
```

This is the non-exclusion test that CPI quantifies, applied to one person. It
skips loci the mixture does not cover and loci where the profile is untyped,
and returns false as soon as one locus shows the profile carries an allele the
stain does not.

If the profile and the mixture share no locus at all it throws rather than
returning `true`. Having no data is not the same as not being excluded, and a
vacuous "included" is exactly the kind of answer that gets misread.

## A worked example

From `examples/matchprobability.j`:

```
mixture over 4 loci
  D3S1358    PI = 0.1572  PE = 0.8428
  vWA        PI = 0.0079  PE = 0.9921
  FGA        PI = 0.0636  PE = 0.9364
  TH01       PI = 0.1342  PE = 0.8658

CPI = 0.000011  - 1 in 94808 could not be excluded
CPE = 0.999989  - the share of the population excluded
```

vWA carries most of the weight because the stain shows only two alleles there;
D3S1358 shows three, so more of the population survives it. A mixture with many
alleles at every locus excludes almost nobody, and its CPI approaches 1.

## Theta is not applied

`inclusionProbability` and `cpi` take no theta. The probability of inclusion
counts how much of the population cannot be excluded - a different question
from a co-ancestry-corrected match probability, and one where the NRC II
sampling formula does not apply. Mixing the two would produce a number that
answers neither question. See [Theta](../concepts/theta.md).

## What CPI does not say

CPI/CPE is a deliberately weak statistic, and its weakness is the source of
most misreporting:

- **It does not say the suspect is a contributor.** It says what fraction of
  the population is not excluded. Everyone in that fraction is equally
  not-excluded.
- **It does not use the suspect's profile.** `cpi` never sees it. Two suspects
  with different profiles, both included, get the same CPI.
- **It ignores peak heights and contributor counts.** A minor contributor whose
  alleles are barely above threshold counts exactly as much as the major one.
- **It degrades badly with dropout.** If a contributor's allele failed to
  amplify, the stain does not show it, and a true contributor is excluded by a
  locus that was never observed properly. CPI has no way to express that, which
  is why guidance generally restricts it to complete, unambiguous mixtures.

## What would be stronger

A likelihood ratio that states both hypotheses explicitly - "the suspect and
one unknown" against "two unknowns" - and models drop-in, dropout and peak
heights. That is full **mixture deconvolution**, the EuroForMix space, and it
is a larger effort layered on this base rather than part of it. It is not in
this deck; see [Limitations](../project/limitations.md).

Where the alternative contributor might be a **relative** rather than a random
person, [Kinship](kinship.md) is the relevant machinery instead.

## Related

- [Profiles and genotypes](../concepts/profiles.md) - the `Mixture` type
- [Match probability](match-probability.md) - the single-source case
- [`forensicgenetics` reference](../reference/forensicgenetics.md)
