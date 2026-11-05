# Mixtures

A mixed stain has no genotype, only the alleles present at each locus. The
core answers the question that needs no assumption about contributors: what
fraction of the population cannot be excluded? `mixture.j` answers the
question that does: how much more probable is the stain if the suspect
contributed? See
[Likelihood ratio with drop-out and drop-in](#likelihood-ratio-with-drop-out-and-drop-in).

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

A likelihood ratio with explicit hypotheses uses the suspect's profile and
models drop-out; see the next section.

## Likelihood ratio with drop-out and drop-in

`mixture.j` computes the semi-continuous (qualitative) LR of LRmix (Gill et
al. 2007; Haned et al. 2012): allele presence, a drop-out probability per
contributor, a drop-in probability, θ. Peak heights are not used.

```jennifer
import "@deckyard/forensicgenetics/mixture.j" as mixture;

# Hp: victim and suspect. Hd: victim and one unknown.
def hp as mixture.Hypothesis init mixture.hypothesis(
    "Hp", [$victim, $suspect], [0.0, 0.1], 0, 0.1);
def hd as mixture.Hypothesis init mixture.hypothesis(
    "Hd", [$victim], [0.0], 1, 0.1);

def lr as float init mixture.mixtureLR($db, $stain, $hp, $hd, 0.05, forensicgenetics.THETA_GENERAL);
def log10Lr as float init mixture.log10MixtureLR($db, $stain, $hp, $hd, 0.05, 0.01);
def perLocus as map of string to float init mixture.locusLRs($db, $stain, $hp, $hd, 0.05, 0.01);
```

`hypothesis(name, known, knownDropout, unknowns, unknownDropout)`: the known
contributors with one drop-out probability each, and a number of unknowns
sharing one. The last two arguments of the LR functions are the drop-in
probability `c` and θ.

### The model

At one locus, with `E` the observed alleles:

- an allele a contributor carries is lost only if every copy drops out:
  probability `d` per copy, so `d²` for a homozygote and `d1 · d2` for an
  allele two contributors share;
- an observed allele nobody carries is a drop-in, `c · p_a` each; no drop-in
  at all is `1 - c`;
- the unknowns' genotypes are summed over, weighted by the Balding-Nichols
  sampling formula conditioned on every known contributor of either
  hypothesis (`mixture.conditioningProfiles(hp, hd)`).

`LR = Π over loci of P(E | Hp) / P(E | Hd)`. Formulae and sources:
[Formulae](../project/formulae.md#mixtures-semi-continuous-likelihood-ratio).

### Checks the model passes

- Suspect alone against one unknown, `d = 0`, `c = 0`: the LR equals
  `forensicgenetics.matchLikelihoodRatio`, at any θ.
- Two unknowns, `d = 0`, `c = 0`, evidence `{a, b}`:
  `(pa + pb)⁴ - pa⁴ - pb⁴`.
- No allele seen, `d` equal for suspect and unknown: LR 1.

### Choosing parameters

- **Number of contributors.** `mixture.minimumContributors(stain)` is the
  floor, `ceil(max alleles at a locus / 2)`. Drop-out can hide contributors;
  compare LRs over plausible counts.
- **Drop-out.** `d = 0` makes a missing suspect allele an exclusion: Hp gives
  probability 0 and the LR is 0 (`log10MixtureLR` throws). LRmix practice is
  to report over a range of `d`, or a conservative value from a drop-out
  sensitivity plot.
- **Drop-in.** `c = 0` makes an allele no contributor carries impossible; if
  Hd also cannot explain a locus, `locusLR` throws.
- **Cost.** Each unknown multiplies the work by the number of genotypes at the
  locus. At most `mixture.MAX_UNKNOWNS` (3) unknowns and
  `mixture.MAX_LOCUS_ALLELES` (10) alleles per locus; a locus over either
  throws. Leave a locus out with `mixture.restrictMixture(stain, keep)`.

### Worked example

`examples/mixture.j` (synthetic table), victim + suspect against victim +
unknown, `c = 0.05`, θ = 0.01:

```
stain needs at least 2 contributors

drop-out   log10 LR
0.05     1.408
0.10     1.672
0.30     2.052

per locus, d = 0.1
  D3S1358    LR = 0.70
  vWA        LR = 37.19
  TH01       LR = 14.15
  D16S539    LR = 0.13
```

D3S1358 lacks the suspect's 18: a drop-out under Hp, so the locus favours Hd.
D16S539 shows a 14 nobody known carries: under Hp it must be a drop-in, under
Hd the unknown can carry it. At `d = 0` the suspect is excluded at D3S1358.

### What the model does not do

- No peak heights, stutter or degradation: those need a continuous model
  (EuroForMix, STRmix).
- One drop-out probability per contributor, the same at every locus.
- Relatives as alternative contributors are not modelled.

See [Limitations](../project/limitations.md#mixtures-stop-at-the-semi-continuous-lr).
