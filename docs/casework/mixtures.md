# Mixtures

A mixed stain has no genotype, only the alleles present at each locus. The
deck answers the question that needs no assumption about contributors: what
fraction of the population cannot be excluded?

```jennifer
def stain as forensicgenetics.Mixture init forensicgenetics.mixture("stain", {
    "D3S1358": ["15", "16", "17"],
    "vWA": ["16", "17"],
    "FGA": ["21", "22", "24"],
    "TH01": ["6", "7", "9.3"]
});
```

Repeated alleles collapse; a locus with an empty list throws.

## Probability of inclusion

A person is not excluded at a locus if both their alleles are among those
observed. Under Hardy-Weinberg:

```
PI = (Σ observed allele frequencies)²     sum capped at 1 before squaring
PE = 1 - PI
```

The frequencies are floored at `5/(2N)`, which is why the cap can matter at a
locus showing many alleles.

```jennifer
forensicgenetics.inclusionProbability($db, "D3S1358", ["15", "16", "17"]);
forensicgenetics.exclusionProbability($db, "D3S1358", ["15", "16", "17"]);

def cpi as float init forensicgenetics.cpi($db, $stain);   # product of PI over loci
def cpe as float init forensicgenetics.cpe($db, $stain);   # 1 - CPI
```

`cpi` throws on a mixture with no loci and on a locus the table lacks.

No theta: PI counts who cannot be excluded, which is not a co-ancestry-corrected
match probability. See [Theta](../concepts/theta.md).

## Is a profile included?

```jennifer
forensicgenetics.isIncluded($suspect, $stain);   # bool
```

The non-exclusion test for one person. Loci only one side carries are skipped;
the first locus where the profile has an allele the stain lacks returns
`false`. A profile and mixture with no shared locus throw instead of returning
`true`.

## Worked example

`examples/matchprobability.j`:

```
mixture over 4 loci
  D3S1358    PI = 0.1572  PE = 0.8428
  vWA        PI = 0.0079  PE = 0.9921
  FGA        PI = 0.0636  PE = 0.9364
  TH01       PI = 0.1342  PE = 0.8658

CPI = 0.000011  - 1 in 94808 could not be excluded
CPE = 0.999989  - the share of the population excluded
```

vWA shows two alleles and excludes the most; a stain showing many alleles at
every locus excludes almost nobody, and its CPI approaches 1.

## What CPI does not say

- It does not say the suspect contributed, only how much of the population is
  not excluded.
- It never sees the suspect's profile: two included suspects get the same CPI.
- It ignores peak heights and the number of contributors.
- It fails under dropout. A true contributor whose allele did not amplify is
  excluded at that locus, which is why guidance restricts CPI to complete,
  unambiguous mixtures.

A likelihood ratio with explicit hypotheses ("suspect and one unknown" against
"two unknowns") that models drop-in, dropout and peak heights is full mixture
deconvolution (the EuroForMix approach). It is not in this deck; see
[Limitations](../project/limitations.md).
