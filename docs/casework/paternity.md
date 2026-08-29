# Paternity

A paternity case is not a pairwise comparison. Father, mother and child are all
available, and conditioning on the mother's genotype is strictly more
informative than treating father and child as an anonymous pair - the mother
tells you which of the child's alleles the father had to supply.

So paternity gets its own trio calculation.

## The paternity index

At one locus, the paternity index is a ratio of two probabilities of the
child's genotype:

```
X = P(child | mother, alleged father)
Y = P(child | mother, a random man from the reference population)
PI = X / Y
```

`X` uses the alleged father's actual transmission probabilities; `Y` uses the
population allele frequencies in their place.

```jennifer
kinship.paternityIndexAt($db, $motherGt, $childGt, $fatherGt, kinship.NO_MUTATION);
```

Across the panel, the **combined paternity index**:

```jennifer
def cpi as float init kinship.paternityIndex($mother, $child, $allegedFather, $db);
```

The loci used are those all three profiles and the database share:

```jennifer
kinship.trioLoci($db, $mother, $child, $allegedFather);
```

## Probability of paternity

```jennifer
def w as float init kinship.probabilityOfPaternity($cpi);            # prior 0.5
def w as float init kinship.probabilityOfPaternityPrior($cpi, 0.1);  # explicit prior
```

`W = CPI / (CPI + 1)` under the conventional neutral prior, or
`W = CPI·π / (CPI·π + (1-π))` for a stated prior `π`.

`W` is a **restatement** of the CPI, not additional evidence, and its
plausible-looking 99.8% hides an assumption that the man was as likely as not to
be the father before any DNA was examined. Quote both numbers, and say which
prior the second one assumes.

## A maternal exclusion is an error

If the child carries neither of the mother's alleles at a locus, `Y` is zero
and no paternity index exists. That is a sample-handling or pedigree problem -
a mislabelled tube, the wrong mother - not a paternity index of zero, so the
deck throws rather than returning a number:

```
locus 'D8S1179': the mother cannot have transmitted either of the child's
alleles - this is a maternal exclusion, not a paternity index
```

## Mutation: surviving one inconsistent locus

STR loci mutate, at roughly one or two per thousand meioses, and almost always
by a single repeat unit. Over a fifteen-locus panel that is common enough to
expect occasionally - so a true father will sometimes fail to explain one
locus.

Without a mutation model, that single locus gives `X = 0`, drives the CPI to
zero, and excludes him. **That is a false exclusion**, and avoiding it is what
the mutation model is for.

```jennifer
def bad as list of string init kinship.inconsistentLoci($mother, $child,
    $allegedFather, $db);
```

lists the loci that are Mendelian-inconsistent. One is the classic single
inconsistency; several mean exclusion, and no mutation model should rescue that.

### The model

The probability that an allele is transmitted as a different one, `k` repeat
steps away in a given direction:

```
P(k steps, one direction) = rate · (1 - decay) · decay^(k-1) / 2
P(no change)              = 1 - rate
```

The rates over all non-zero changes sum to `rate`, and one-step changes
dominate. A microvariant difference smaller than a full repeat - 15 to 15.2 -
counts as one step.

```jennifer
kinship.mutationProbability(kinship.DEFAULT_MUTATION, "D8S1179", "16", "15");
```

Two constants:

- **`NO_MUTATION`** - strict Mendelian transmission. One inconsistency
  excludes.
- **`DEFAULT_MUTATION`** - `rate = 0.002`, `decay = 0.1`. This is a
  **placeholder** in the right order of magnitude, not a published figure.

`paternityIndex` uses `DEFAULT_MUTATION`; `paternityIndexWith` takes an
explicit model.

### Use real rates

Observed mutation rates vary over roughly an order of magnitude between loci,
and a rescued trio's CPI is directly proportional to the rate assumed at the
inconsistent locus. Substitute published per-locus rates:

```jennifer
def model as kinship.MutationModel init kinship.DEFAULT_MUTATION;
$model = kinship.withLocusRate($model, "D8S1179", 0.0016);
$model = kinship.withLocusRate($model, "FGA", 0.0028);

def cpi as float init kinship.paternityIndexWith($mother, $child,
    $allegedFather, $db, $model);
```

The mutation model applies to the **alleged father's transmission only** - it
is his transmission that is in question. With `NO_MUTATION` the calculation
reduces exactly to the Mendelian one.

## A worked example

From `examples/paternity.j`. The alleged father is 16/16 at D8S1179 and the
child needs a 15:

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
```

Three things to read off it. Every consistent locus loses about 0.1% to the
`(1 - rate)` discount - mutation is not free even where nothing mutated. The
inconsistent locus contributes 0.0219 instead of 0, which costs roughly two
orders of magnitude but does not exclude. And the choice of prior moves `W`
from 99.78% to 98.08% on identical evidence.

## When the father is unavailable

Everything above needs the alleged father's profile. When he is dead or
missing, the case rests on his relatives - his parents, his other children -
and that is a pedigree question, not a trio one. See
[General pedigrees](pedigrees.md), where the same example is worked through
his parents instead.

## Related

- [Kinship](kinship.md) - the general pairwise machinery
- [General pedigrees](pedigrees.md) - when the father is unavailable
- [`kinship` reference](../reference/kinship.md)
