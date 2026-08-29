# General pedigrees

[Kinship](kinship.md) answers questions about a pair. That approach breaks down
the moment a third person is involved - and in practice a third person almost
always is.
A missing person is identified through two siblings and a niece. An alleged
father is dead, and only his parents can be typed. A disaster victim is matched
to a family group rather than to an individual.

Those need the likelihood of the **whole pedigree**.

## The problem, and peeling

The likelihood of a pedigree at one locus is a sum over every genotype every
untyped member could have. Taken literally that is exponential in the number of
untyped people - a handful of them and it stops being computable.

**Elston-Stewart peeling** makes it linear. The pedigree is a tree of
individuals and nuclear families; each branch is summarised - *peeled* - into a
vector over the genotypes of the single person connecting it to the rest. The
sums move inward instead of nesting. This module peels by passing those vectors
along the tree's edges.

Two consequences you will meet:

- **Only simple pedigrees.** Peeling needs the individual/family graph to be a
  forest. Loops are rejected, not mis-summed. See below.
- **Allele lumping.** Alleles no member's genotype shows are pooled into one
  class carrying their combined frequency. Because untyped members are
  unobserved those alleles are exchangeable, so the pooling is **exact**, not
   an approximation - and this keeps a twenty-allele locus down to a
  handful of genotype classes.

## Building a pedigree

```jennifer
def ped as pedigree.Pedigree init pedigree.pedigree("case 42");
$ped = pedigree.withFounder($ped, "GF", "m");
$ped = pedigree.withFounder($ped, "GM", "f");
$ped = pedigree.withChild($ped, "FA", "m", "GF", "GM");
```

A **founder** has no parents in the pedigree; their genotype comes from the
population. Everyone else names **both** parents - a member with only one
recorded parent is rejected, because the missing parent is still a source of
alleles and leaving them out changes the answer.

`sex` is `"m"`, `"f"` or `""`. It is recorded but not used: the autosomal
calculation is symmetric in the two parental roles.

Inspecting:

```jennifer
pedigree.members($ped);            # ids, insertion order
pedigree.isFounder($ped, "GF");    # bool
pedigree.childrenOf($ped, "GF");   # ids
pedigree.personOf($ped, "FA");     # the Person
```

## Validation

```jennifer
pedigree.validate($ped);      # throws, naming the problem
pedigree.isSimple($ped);      # bool
```

`validate` checks that every named parent exists, that nobody is their own
ancestor, and that the individual/family graph is a forest. Every likelihood
function calls it first.

### Loops

A consanguineous mating closes a loop:

```jennifer
# Two full siblings have a child together.
$p = pedigree.withChild($p, "A", "m", "GF", "GM");
$p = pedigree.withChild($p, "B", "f", "GF", "GM");
$p = pedigree.withChild($p, "C", "m", "A", "B");    # loop
```

```
pedigree 'loop' contains a loop (a consanguineous mating or a repeated
marriage); Elston-Stewart peeling needs a simple pedigree
```

Handling loops needs cutset conditioning, which is not implemented. Rejecting
them is the honest alternative to returning a number that is quietly wrong.

A half-sibling pedigree - one father, two mothers - is **not** a loop, and the
check does not over-reject it.

## Likelihoods

At one locus, with a genotype map for the typed members:

```jennifer
def l as float init pedigree.locusLikelihood($ped, $genotypes, $db, "TH01");
```

Members absent from `genotypes` are untyped and summed out; ids that are not
pedigree members are ignored.

Across a panel, from profiles keyed by member id:

```jennifer
def loci as list of string init pedigree.pedigreeLoci($ped, $profiles, $db);
def l as float init pedigree.pedigreeLikelihood($ped, $profiles, $db);
def logL as float init pedigree.log10PedigreeLikelihood($ped, $profiles, $db);
```

A locus counts if the database covers it and at least one pedigree member is
typed there. `pedigreeLikelihood` underflows on a long panel - prefer the
log form, or take a ratio.

## The likelihood ratio

State both hypotheses as pedigrees over the same people and divide:

```jennifer
def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
```

The ratio is taken locus by locus, so it stays in range where either likelihood
alone would underflow.

Both pedigrees must have the **same typed members**, and `pedigreeLR` refuses
otherwise - scoring two hypotheses against different data makes their ratio
meaningless. Untyped members may differ freely, and that is exactly how a
hypothesis introduces an unknown alternative father.

```jennifer
pedigree.typedMembers($ped, $profiles);   # who actually contributes
```

## A worked example: deficient paternity

The alleged father is unavailable, so the case is made through his parents.
Four people are typed - the two grandparents, the mother and the child - and
the man himself never is.

```
  H1                              H2
    GF x GM                         GF   GM   (unrelated bystanders)
       |                                 UN x MO
      FA x MO                                |
           |                                CH
          CH
```

```jennifer
def h1 as pedigree.Pedigree init pedigree.pedigree("son of GF x GM");
$h1 = pedigree.withFounder($h1, "GF", "m");
$h1 = pedigree.withFounder($h1, "GM", "f");
$h1 = pedigree.withFounder($h1, "MO", "f");
$h1 = pedigree.withChild($h1, "FA", "m", "GF", "GM");
$h1 = pedigree.withChild($h1, "CH", "m", "FA", "MO");

def h2 as pedigree.Pedigree init pedigree.pedigree("unrelated man");
$h2 = pedigree.withFounder($h2, "GF", "m");
$h2 = pedigree.withFounder($h2, "GM", "f");
$h2 = pedigree.withFounder($h2, "MO", "f");
$h2 = pedigree.withFounder($h2, "UN", "m");
$h2 = pedigree.withChild($h2, "CH", "m", "UN", "MO");

def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
```

`examples/pedigree.j` prints it locus by locus:

```
locus       L(H1)          L(H2)          LR
D3S1358      0.0000210529 0.0000111749 1.8839
vWA          0.0000030558 0.0000026145 1.1688
FGA          0.0000002152 0.0000000707 3.0451
D8S1179      0.0000001059 0.0000000174 6.0827
D21S11       0.0000000211 0.0000000187 1.1302
D18S51       0.0000026876 0.0000008977 2.9940
TH01         0.0000013055 0.0000004047 3.2258
D16S539      0.0000003625 0.0000003200 1.1328

combined LR = 504.29

for contrast, the pairwise grandparent LR (GF vs CH alone) = 53.67
```

The last line is the point of the module. Asking the pairwise question - *is
the child a grandchild of GF?* - throws GM away entirely, and with her the
constraint that the two grandparents must **jointly** explain the child's
paternal allele. Using both at once is worth roughly ten times as much.

## Correctness

Peeling is easy to get subtly wrong, so the test overlay checks it against
`kinship.j` rather than against itself. A pedigree LR built by peeling must
reproduce the closed-form result that module computes an entirely different
way, and the suite asserts that for parent/child, full siblings, half siblings,
grandparent, avuncular, and the trio paternity index - all to the last digit.

Agreement between two independent derivations is much stronger evidence than
either matching a number captured from itself. See
[Testing and contributing](../project/testing.md).

## Related

- [Kinship](kinship.md) - the pairwise form, and where it runs out
- [Paternity](paternity.md) - the trio, when the father is available
- [`pedigree` reference](../reference/pedigree.md)
