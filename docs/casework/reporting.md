# Reporting a match

A match probability's allele frequencies are estimates from a database of a
few hundred or thousand people. `reporting.j` implements three published
allowances for that sampling error, and the ENFSI verbal scale for the
likelihood ratio.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/reporting.j" as reporting;
```

## The size-bias correction

The suspect and the offender share the profile; before taking frequencies,
add both to the database. A heterozygote allele seen `x` times among `2N`
becomes `(x + 2) / (2N + 4)`, a homozygote's `(x + 4) / (2N + 4)`:

```jennifer
def lr as float init reporting.sizeBiasLikelihoodRatio($db, $suspect,
    forensicgenetics.THETA_GENERAL);
def logRmp as float init reporting.log10SizeBiasRandomMatchProbability($db, $suspect,
    forensicgenetics.THETA_GENERAL);
```

The correction replaces the `5/(2N)` floor: an allele the table lacks scores
`2 / (2N + 4)`. It vanishes as `N` grows. `sizeBiasAlleleFrequency` shows one
allele's adjusted frequency.

## The Bayesian estimator

FSR-GUI-0012 section 11.4 takes each allele's posterior mean under a Dirichlet
prior as its frequency: `(x + 1/k) / (2N + 1)` with `k` the alleles the table
lists (11.4.1), or `(x + 1/(k + 2)) / (2N + 1)` to allow for unobserved
alleles (11.4.2):

```jennifer
def lr as float init reporting.bayesianLikelihoodRatio($db, $suspect,
    forensicgenetics.THETA_GENERAL, reporting.PRIOR_ONE_OVER_K_PLUS_TWO);
```

An allele the table lacks scores `(1/(k + 2)) / (2N + 1)`, smaller than the
size-bias `2 / (2N + 4)`. `bayesianFrequency` shows one allele's estimate;
`PRIOR_UNIFORM` gives `(x + 1) / (2N + k)`. Unlike the size-bias correction,
the estimate does not depend on the profile being scored.

## A posterior bound on the RMP

The Bayesian alternative treats each locus's frequencies as unknown, with a
Dirichlet posterior given the database counts, and reports a quantile of the
RMP's posterior instead of a point value. Curran, Buckleton, Triggs & Weir
(2002) preferred it among the methods they compared. A 99% one-sided bound:

```jennifer
def bound as float init reporting.likelihoodRatioLowerBound($db, $suspect,
    forensicgenetics.THETA_GENERAL, reporting.PRIOR_ONE_OVER_K, 0.99,
    reporting.DEFAULT_DRAWS, 1);
```

The LR exceeds `bound` with posterior probability 0.99. For several quantiles
from one simulation, or for a profile whose LR passes float range:

```jennifer
def q as list of float init reporting.log10RmpQuantiles($db, $suspect,
    forensicgenetics.THETA_GENERAL, reporting.PRIOR_UNIFORM, [0.5, 0.99],
    reporting.DEFAULT_DRAWS, 1);
```

The prior is the caller's choice: `PRIOR_UNIFORM` (one pseudo-count per
allele), `PRIOR_ONE_OVER_K` (`1/k`; Curran & Buckleton 2011 compare the
two) or `PRIOR_ONE_OVER_K_PLUS_TWO` (`1/(k + 2)`, FSR-GUI-0012 11.4.2). Here
`PRIOR_ONE_OVER_K`'s `k` also counts profile alleles the table lacks, so the
prior stays a proper Dirichlet of mass 1; its posterior mean matches 11.4.1
only when every profile allele is listed. The simulation is seeded; report
the seed and the number of draws with the bound. The quantile's Monte Carlo error is about `sqrt(q (1 - q) / draws)`
on the probability scale, so a 99% bound at 10,000 draws rests on about 100
draws in the tail.

## Verbal equivalents

```jennifer
def v as reporting.VerbalEquivalent init reporting.verbalEquivalent($lr);
io.printf("%s\n", $v.statement);
```

The scale is the illustration in the ENFSI Guideline for Evaluative Reporting
in Forensic Science (2015), guidance note 4:

| LR | Verbal equivalent |
| --- | --- |
| 1 | The forensic findings do not support one proposition over the other. |
| 2 - 10 | ... provide weak support for the first proposition relative to the alternative. |
| 10 - 100 | ... provide moderate support for the first proposition rather than the alternative |
| 100 - 1000 | ... moderately strong support ... |
| 1000 - 10,000 | ... strong support ... |
| 10,000 - 1,000,000 | ... very strong support ... |
| 1,000,000 and above | ... extremely strong support ... |

An LR below 1 reads `1/LR` on the same scale as support for the alternative.
The guideline does not say which end of a range is inclusive and calls the
ranges indicative; the deck closes the upper end (`(10, 100]`, mirrored
`[0.01, 0.1)`) and puts LRs in `(1, 2)` in the weak band. `VerbalEquivalent`
carries the band's bounds and which ends belong to it, so a report can state
them. For an LR held as a logarithm (`mixture.log10MixtureLR`) use
`verbalEquivalentLog10`.

The guideline asks that a laboratory use one scale for all its evidence types,
and that the alternative proposition be stated with the conclusion: weak
support for the first proposition is not support for the alternative.

## Worked example

`examples/reporting.j` scores a five-locus profile on the synthetic example
table: the plug-in LR, the size-bias LR, the Bayesian LR, the 99% bound, and
their verbal equivalents.
