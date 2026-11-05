# Mutation models and silent alleles

`kinship.j` has one mutation model: an unbounded stepwise ladder, applied to
the alleged father only, the same at every locus unless overridden. The
[`mutation`](../reference/mutation.md) module adds the models Familias offers,
separate paternal and maternal models, published per-locus rates, and silent
(null) alleles. It answers the same trio, duo and pairwise questions, and
with no mutation and no silent allele it returns exactly what `kinship.j`
returns.

## From a model to a matrix

A `MutationSpec` names a model; at a locus it becomes a matrix over that
locus's alleles, `M[i][j]` = P(allele `i` is transmitted as `j`), every row
summing to 1:

```jennifer
def la as mutation.LocusAlleles init mutation.locusAlleles($db, "D8S1179", [], 0.0);
def m as mutation.MutationMatrix init mutation.mutationMatrix(mutation.stepwiseModel(0.002, 0.1), $la);
mutation.transitionProbability($m, "16", "15");
mutation.overallRate($m, $la);      # Σ π_i (1 - M_ii)
mutation.isStationary($m, $la);     # π M = π ?
```

| Model | Mutation goes to | Notes |
| --- | --- | --- |
| `equalModel(rate)` | every other allele equally | not stationary unless frequencies are equal |
| `proportionalModel(rate)` | allele `j` in proportion to `p_j` | stationary by construction; ignores repeat length |
| `stepwiseModel(rate, range)` | nearby repeat counts, `range` per extra step | microvariants keep to their own class |
| `extendedStepwiseModel(rate, range, rate2)` | as Stepwise, plus `rate2` between classes (15 -> 15.2) | Familias 3's model |

The matrix is normalised over the alleles the table lists (plus any typed
allele it lacks), as Familias does. `kinship.mutationProbability` is not: its
ladder runs on past the table's alleles, so the two stepwise models give
different numbers at the same rate.

## Stationary or not

A model is stationary when one generation of mutation leaves the allele
frequencies unchanged, `π M = π`. Equal and Stepwise generally are not, which
inflates or deflates rare alleles in untyped generations (Simonsson & Mostad
2016). `stationaryModel` applies their PM transformation: each row off the
diagonal is scaled so that `π` becomes stationary and the expected rate stays
the same.

```jennifer
def spec as mutation.MutationSpec init mutation.stationaryModel(mutation.extendedStepwiseModel(0.002, 0.1, 0.0001));
```

PM does not always exist. A rare allele (a typed allele at the `5/(2N)` floor,
say) can need a scaled diagonal below zero, and a Stepwise matrix with
microvariants and `rate2` 0 has no unique stationary distribution; both throw,
as does a locus with more than `MAX_STATIONARY_ALLELES` (200) alleles.
In a trio or duo every transmitting parent is typed, so stationarity matters
only through the changed matrix entries.

## Sex-specific models

```jennifer
def pm as mutation.ParentageModel init mutation.parentageModel(
    mutation.stepwiseModel(0.002, 0.1),       # paternal
    mutation.stepwiseModel(0.0005, 0.1));     # maternal
def cpi as float init mutation.paternityIndex($mother, $child, $allegedFather, $db, $pm, 0.0);
mutation.paternityIndexAt($db, $motherGt, $childGt, $fatherGt, $pm, 0.0);
```

The maternal model acts on the mother's transmission in both `X` and `Y`, so
a child carrying neither of the mother's alleles no longer throws a maternal
exclusion once the maternal rate is above 0: the index weighs a maternal
mutation instead. With `parentageModel(spec, mutation.NO_MUTATION)` the
mother is Mendelian, as in `kinship.j`. The random man of the alternative does
not mutate. `duoPaternityIndex` uses the paternal model only.

## Published per-locus rates

`aabb2003Counts` is the AABB 2003 summary of apparent mutations in paternity
testing at 18 loci, as NIST STRBase published it (source and date in
[Formulae and sources](../project/formulae.md#published-mutation-counts)).
`aabb2003Rates` turns it into rates:

| `which` | Rate | Caveat |
| --- | --- | --- |
| `"total"` | all mutations / all meioses (STRBase's column) | sex-averaged |
| `"paternal"` | mutations assigned to the father / paternal meioses | leaves out the unassigned ones: a lower bound |
| `"maternal"` | the same for the mother | as above |

```jennifer
def rates as map of string to float init mutation.aabb2003Rates("paternal");
def spec as mutation.MutationSpec init mutation.withLocusRates(mutation.stepwiseModel(0.002, 0.1), $rates);
def km as kinship.MutationModel init mutation.kinshipModelWithRates(kinship.DEFAULT_MUTATION, $rates);
```

Loci not in the table keep the blanket rate. The names follow
`panels.normalizeLocus` (`vWA`, `Penta D`); rates for other loci, or from
another source, go in with `withLocusRate` / `withLocusRates` (or
`kinship.withLocusRate`).

## Silent alleles

A primer-binding-site variant can stop an allele from amplifying. A person
showing a single peak `a` may then be `a/silent`, and a father and child who
look like opposite homozygotes (`a` and `b`) need not exclude. Give the locus a
silent frequency `q`:

```jennifer
def pm as mutation.ParentageModel init mutation.withSilentFrequency(mutation.MENDELIAN, "TH01", 0.001);
mutation.duoPaternityIndexAt($db, $childGt, $fatherGt, $pm, 0.0);
mutation.kinshipLRAt($db, $ga, $gb, kinship.FULL_SIBLINGS, kinship.UNRELATED, 0.001, 0.0);
mutation.kinshipLR($a, $b, kinship.FULL_SIBLINGS, kinship.UNRELATED, $db, {"TH01": 0.001}, 0.0);
```

The silent allele is one more allele with frequency `q`; the table's
frequencies are read as frequencies among detectable alleles and scaled by
`1 - q`. Every single-peak genotype is summed over `a/a` and `a/silent`
(Familias' treatment). For opposite homozygotes, child `a` and alleged father
`b`, the duo index is `q / ((πa + 2q)(πb + 2q))`; in a trio with a
heterozygous mother `a/c`, `q / ((πb + 2q)(πa + q))`. Both tend to 0 with
`q`, and at `q = 0` every function equals its `kinship.j` original. The
silent allele does not mutate.

## Worked example

`examples/mutation.j` runs the trio of `examples/paternity.j` (synthetic
frequencies) through each model at a paternal rate of 0.002:

```
D8S1179  mother 12/13  child 13/15  alleged father 16/16

paternal model                PI at D8S1179   CPI
Equal                         0.0060827       127.5
Proportional                  0.0023358       48.9
Stepwise (range .1)           0.0437956       923.5
Extended stepwise             0.0437956       923.4
Extended stepwise, PM         0.0182257       384.8
kinship.j stepwise ladder     0.0218978

Extended stepwise             stationary false  expected rate 0.002000
Extended stepwise, PM         stationary true  expected rate 0.002000

AABB 2003 paternal rate at D8S1179 = 0.00159;  CPI = 734.3

duo at TH01, child 7, alleged father 9.3 (single peaks)
silent frequency 0.000  PI 0.00000
silent frequency 0.001  PI 0.10893
silent frequency 0.010  PI 0.77990
log10 of the parent/child LR with q = .01: -0.108
```

The one-step change 16 -> 15 is what Stepwise favours; Equal spreads the same
rate over all nine D8S1179 alleles, and Proportional sends it mostly to
common alleles (15 is at 0.0411). The normalised Stepwise row puts about
twice the ladder's weight on the one-step change, because the ladder halves
`rate` between the two directions and loses mass past the table's ends. PM
rescales each row off the diagonal by its `d_i`; here the 16 row by 0.42, and
the CPI falls by a factor of 2.4 on the same expected rate. The model is a
case decision; report it.
