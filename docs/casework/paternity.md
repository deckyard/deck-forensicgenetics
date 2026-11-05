# Paternity

With mother, child and alleged father all typed, conditioning on the mother
tells you which of the child's alleles the father had to supply. That is more
informative than a pairwise father/child LR, so paternity has its own trio
calculation.

## The paternity index

At one locus:

```
X  = P(child | mother, alleged father)
Y  = P(child | mother, random man from the reference population)
PI = X / Y
```

`X` uses the alleged father's transmission probabilities, with mutation; `Y`
uses population allele frequencies. The full formula is in
[Formulae and sources](../project/formulae.md#paternity-the-trio-index).

```jennifer
kinship.paternityIndexAt($db, $motherGt, $childGt, $fatherGt, kinship.NO_MUTATION);

def cpi as float init kinship.paternityIndex($mother, $child, $allegedFather, $db);
kinship.trioLoci($db, $mother, $child, $allegedFather);   # the loci used
```

The combined paternity index (CPI) is the product over the loci all three
profiles carry. A locus one of them lacks is skipped; a locus all three carry
but the table does not cover throws, because dropping it could hide an
exclusion - the one locus where the alleged father cannot explain the child
is exactly the one that must not vanish. Narrow the profiles with
`forensicgenetics.restrict` to leave a locus out on purpose. No shared locus
throws. `inconsistentLoci` follows the same rule.

## Probability of paternity

```jennifer
def w as float init kinship.probabilityOfPaternity($cpi);              # prior 0.5
def w01 as float init kinship.probabilityOfPaternityPrior($cpi, 0.1);  # prior 0.1
```

`W = CPI / (CPI + 1)` at prior 0.5, or `W = CPI·π / (CPI·π + (1-π))` for a
prior `π` in `(0, 1)`. `W` restates the CPI and adds no evidence; a 99.8%
rests on assuming the man was as likely as not the father before any DNA was
examined. Report both numbers and the prior.

## Maternal exclusion throws

If the child carries neither of the mother's alleles at a locus, `Y` is zero
and no index exists. That points to a sample or pedigree error (a mislabelled
tube, the wrong mother), so the call throws:

```
locus 'D8S1179': the mother cannot have transmitted either of the child's
alleles - this is a maternal exclusion, not a paternity index
```

## Mutation

STR loci mutate at roughly one or two per thousand meioses, almost always by
one repeat. Over a fifteen-locus panel a true father will occasionally fail to
explain one locus. Without a mutation model that locus gives `X = 0` and a
false exclusion.

```jennifer
def bad as list of string init kinship.inconsistentLoci($mother, $child,
    $allegedFather, $db);
```

lists the Mendelian-inconsistent loci. One is the classic single
inconsistency; several mean exclusion, and no mutation model should rescue
that.

The stepwise model (formula in
[Formulae and sources](../project/formulae.md#stepwise-mutation-model)) gives
`1 - rate` for no change and `rate · (1 - decay) · decay^(k-1) / 2` for a
`k`-step change in one direction. A microvariant difference below one repeat
(15 to 15.2) is one step. It applies to the alleged father's transmission
only; with `NO_MUTATION` the index is the strict Mendelian one.

```jennifer
kinship.mutationProbability(kinship.DEFAULT_MUTATION, "D8S1179", "16", "15");
```

| Model | rate | decay | Use |
| --- | --- | --- | --- |
| `NO_MUTATION` | 0.0 | 0.0 | strict Mendelian; one inconsistency excludes |
| `DEFAULT_MUTATION` | 0.002 | 0.1 | used by `paternityIndex`; a placeholder, not a published rate |

A rescued trio's CPI is proportional to the rate assumed at the inconsistent
locus, and observed rates vary about tenfold between loci. Set published
per-locus rates:

```jennifer
def model as kinship.MutationModel init kinship.DEFAULT_MUTATION;
$model = kinship.withLocusRate($model, "D8S1179", 0.0016);
$model = kinship.withLocusRate($model, "FGA", 0.0028);

def cpi as float init kinship.paternityIndexWith($mother, $child,
    $allegedFather, $db, $model);
```

Published per-locus rates (AABB 2003, via NIST STRBase), Familias' Equal,
Proportional and Extended stepwise models, a maternal model and silent
alleles are in [Mutation models and silent alleles](mutation-models.md).
`mutation.kinshipModelWithRates` puts the published rates into a
`kinship.MutationModel`.

## Worked example

`examples/paternity.j`. The alleged father is 16/16 at D8S1179 and the child
needs a 15:

```
locus        mother  child    father     PI strict  PI mutation
D3S1358      15/16    15/17    17/18      3.7679 3.7638
vWA          16/17    17/19    19/19      4.6751 4.6657
FGA          21/22    22/24    24/25      6.0901 6.0834
D8S1179      12/13    13/15    16/16      0.0000 0.0219
D21S11       29/30    30/31.2  31.2/32.2  2.2604 2.2579
D18S51       14/15    15/17    17/18      5.9880 5.9814
TH01         6/7      7/9.3    9.3/10     6.4516 6.4445
D16S539      11/12    12/13    13/14      2.2655 2.2630

Mendelian-inconsistent loci: ["D8S1179"]

strict Mendelian CPI  = 0.000000
  -> excluded: no allowance is made for mutation

CPI with mutation     = 460.77
W (prior 0.5)         = 0.99783441
W (prior 0.1)         = 0.98084163

theta  trio CPI  duo PI
0.00     460.77    1.89
0.01     345.16    1.31
0.03     214.83    0.71
```

Each consistent locus loses about 0.1% to the `(1 - rate)` factor. The
inconsistent locus contributes 0.0219 instead of 0. The prior moves `W` from
99.78% to 98.08% on the same evidence.
At `θ = 0.03` the trio CPI falls by half. The duo index is far smaller:
without the mother, either of the child's alleles may be paternal, and
D8S1179 still needs a mutation.

## With theta

```jennifer
def cpi as float init kinship.paternityIndexTheta($mother, $child, $allegedFather,
    $db, forensicgenetics.THETA_GENERAL);
def cpi2 as float init kinship.paternityIndexWithTheta($mother, $child,
    $allegedFather, $db, $model, 0.01);
kinship.paternityIndexAtTheta($db, $motherGt, $childGt, $fatherGt, $model, 0.01);
```

`X` is unchanged: under paternity the child's paternal allele is the alleged
father's, transmitted with the mutation model, and no population draw is
involved. In `Y` the random man's allele `c` comes from the same
sub-population as the typed mother and alleged father, so it is drawn given
their four alleles:

```
F(c) = (m θ + (1-θ) pc) / (1 + 3θ)      m = copies of c among the mother's and alleged father's alleles
```

(Balding & Nichols 1995; Buckleton et al. 2016, ch. 11.) Mother, child and
alleged father all `aa` give `(1 + 3θ) / (4θ + (1-θ) pa)`. Theta and the
mutation model combine without interaction: one acts only on `X`, the other
only on `Y`. At `theta = 0` each function returns exactly its theta-0
counterpart. The random man is assumed not to mutate, as in the theta-0 index.

## Duos

When the mother is not typed, the duo index weighs the alleged father and
child alone:

```jennifer
def pi as float init kinship.duoPaternityIndex($child, $allegedFather, $db,
    kinship.DEFAULT_MUTATION, 0.01);
kinship.duoPaternityIndexAt($db, $childGt, $fatherGt, kinship.NO_MUTATION, 0.0);
```

`X` takes one of the child's alleles from the alleged father (with mutation)
and draws the other from the population given his two; `Y` draws both given
his two and then the first. With `NO_MUTATION` it equals
`kinshipLRTheta(child, father, PARENT_CHILD, db, theta)`; theta 0 is allowed
and gives the product-rule duo index. It multiplies over the loci both carry,
with the same rule for uncovered loci as the trio. A duo has no maternal
exclusion; it carries less information than a trio, so prefer the trio
whenever the mother is typed. Formulae in
[Formulae and sources](../project/formulae.md#kinship-and-paternity-with-theta).

When the alleged father cannot be typed, the case goes through his relatives:
see [General pedigrees](pedigrees.md).
