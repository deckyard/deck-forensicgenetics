# General pedigrees

[Kinship](kinship.md) handles a pair. A missing person identified through two
siblings and a niece, or an alleged father known only through his parents,
needs the likelihood of the whole pedigree.

## Peeling

The likelihood at one locus sums over every genotype of every untyped member,
which is exponential if done literally. Elston-Stewart peeling makes it linear
in the number of people: each branch of the individual/family tree is reduced
to a vector over the genotypes of the person joining it to the rest, and these
vectors are passed along the tree's edges. The recursion is in
[Formulae and sources](../project/formulae.md#pedigree-likelihood).

- Simple pedigrees only: the individual/family graph must be a forest. Loops
  throw (see [Loops](#loops)).
- Alleles no typed member shows are lumped into one class with their combined
  frequency. Untyped members' alleles are exchangeable, so this is exact, and
  it keeps a twenty-allele locus to a few genotype classes.
- `θ = 0`.

## Building a pedigree

```jennifer
def ped as pedigree.Pedigree init pedigree.pedigree("case 42");
$ped = pedigree.withFounder($ped, "GF", "m");
$ped = pedigree.withFounder($ped, "GM", "f");
$ped = pedigree.withChild($ped, "FA", "m", "GF", "GM");
```

A founder's genotype prior is the population frequency. Everyone else names
both parents: a member with one parent recorded throws, because the missing
parent still contributes alleles. Parents may be added after the child.

`sex` is `"m"`, `"f"` or `""`. It is recorded and not used: the autosomal
calculation treats both parents alike.

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

`validate` checks that every named parent is a member, that nobody is their
own ancestor, and that the individual/family graph is a forest. Every
likelihood function calls it first.

### Loops

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

Loops need cutset conditioning, which is not implemented. A half-sibling
pedigree (one father, two mothers) is not a loop.

## Likelihoods

```jennifer
def l as float init pedigree.locusLikelihood($ped, $genotypes, $db, "TH01");

def loci as list of string init pedigree.pedigreeLoci($ped, $profiles, $db);
def lp as float init pedigree.pedigreeLikelihood($ped, $profiles, $db);
def logL as float init pedigree.log10PedigreeLikelihood($ped, $profiles, $db);
```

`$genotypes` maps member id to `Genotype` at one locus; `$profiles` maps member
id to `Profile`. Members without an entry are untyped and summed out; ids that
are not members are ignored. A locus is used if the table covers it and at
least one member is typed there. `pedigreeLikelihood` underflows on a long
panel; use the log form or a ratio.

## The likelihood ratio

```jennifer
def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
pedigree.typedMembers($h1, $profiles);   # who contributes data
```

The ratio is taken locus by locus, so it stays in range where either
likelihood alone would underflow. Both pedigrees must type the same members,
or the call throws; untyped members may differ, which is how a hypothesis adds
an unknown alternative father. A denominator likelihood of zero at any locus
throws.

## Worked example: deficient paternity

The alleged father is unavailable. The grandparents, the mother and the child
are typed.

```
  H1                              H2
    GF x GM                         GF   GM   (unrelated bystanders)
       |                                 UN x MO
      FA x MO                                |
           |                                CH
          CH
```

The pedigrees are built as in the [Quick start](../getting-started/quick-start.md#5-a-pedigree).
`examples/pedigree.j` prints:

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

The pairwise question (is CH a grandchild of GF?) ignores GM and the
constraint that both grandparents must jointly explain the child's paternal
allele. Using both is worth about ten times as much here.

## Correctness

The test overlay checks peeling against `kinship.j`'s closed forms, which share
no code with it: parent/child, full siblings, half siblings, grandparent,
avuncular and the trio paternity index, to a relative tolerance of 1e-9. See
[Testing and contributing](../project/testing.md#cross-validation-between-modules).
