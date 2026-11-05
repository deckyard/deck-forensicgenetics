# X-chromosomal STRs

X-STRs help where autosomal markers cannot separate the hypotheses, typically
deficiency cases that involve a female: are two women paternal half sisters,
or is this woman the paternal grandmother of that girl? A man has one X,
from his mother, and passes it whole to every daughter and to no son. So:

- a father and daughter share an X allele IBD at every locus, as do paternal
  half sisters and a paternal grandmother with her granddaughter;
- full sisters always share their paternal X;
- father and son, and paternal half brothers, share nothing on the X.

## Profiles carry their sex

```jennifer
import "@deckyard/forensicgenetics/xstr.j" as xstr;

def man as xstr.XProfile init xstr.maleProfile("AF", {"DXS7132": "13", "HPRTB": "12"});
def woman as xstr.XProfile init xstr.femaleProfile("C", [
    forensicgenetics.genotype("DXS7132", "13", "14"),
    forensicgenetics.genotype("HPRTB", "12", "13")
]);
```

A male genotype is stored with `a == b`; it is hemizygous, not homozygous,
and every function reads it by the profile's `sex`.

## Match probability

```jennifer
def model as xstr.XFrequencies init xstr.independentModel($db);
def pm as float init xstr.matchProbability($model, $man);     # Π p
def pf as float init xstr.matchProbability($model, $woman);   # Π 2pq or p²
```

`$db` is an ordinary `FrequencyDb` of X-STR allele frequencies, floored at
`5/(2N)` as in the core. `log10MatchProbability` sums logarithms. Theta is 0.

## Kinship

A relationship is a set of X-IBD coefficients for an ordered pair of sexes
(Pinto, Gusmão & Amorim 2011). Two females have `(k0, k1, k2)` as on an
autosome; a pair including a male has only `k0` and `k1`.

| Constant | A - B | (k0, k1, k2) |
| --- | --- | --- |
| `FATHER_DAUGHTER` | M - F | (0, 1, 0) |
| `FATHER_SON` | M - M | (1, 0, 0) |
| `MOTHER_SON` | F - M | (0, 1, 0) |
| `MOTHER_DAUGHTER` | F - F | (0, 1, 0) |
| `FULL_SISTERS` | F - F | (0, 1/2, 1/2) |
| `FULL_BROTHERS` | M - M | (1/2, 1/2, 0) |
| `BROTHER_SISTER` | M - F | (1/2, 1/2, 0) |
| `PATERNAL_HALF_SISTERS` | F - F | (0, 1, 0) |
| `MATERNAL_HALF_SISTERS` | F - F | (1/2, 1/2, 0) |
| `MATERNAL_HALF_BROTHERS` | M - M | (1/2, 1/2, 0) |
| `PATERNAL_GRANDMOTHER_GRANDDAUGHTER` | F - F | (0, 1, 0) |
| `MATERNAL_GRANDMOTHER_GRANDDAUGHTER` | F - F | (1/2, 1/2, 0) |
| `MATERNAL_GRANDFATHER_GRANDDAUGHTER` | M - F | (1/2, 1/2, 0) |
| `PATERNAL_AUNT_NIECE` | F - F | (1/2, 1/2, 0) |
| `PATERNAL_UNCLE_NIECE` | M - F | (1/2, 1/2, 0) |
| `MATERNAL_AUNT_NIECE` | F - F | (1/4, 3/4, 0) |

Other pedigrees: `xRelationship(name, sexA, sexB, k0, k1, k2)`.
`coancestry(rel)` gives the X coancestry coefficient; the tests check every
constant against Gomes et al. (2020), Table 5.

```jennifer
def lr as float init xstr.kinshipLR($a, $b, xstr.PATERNAL_HALF_SISTERS, $model);
def vs as float init xstr.kinshipLRAgainst(
    $a, $b, xstr.PATERNAL_HALF_SISTERS, xstr.MATERNAL_HALF_SISTERS, $model);
```

Autosomally, paternal and maternal half sisters have the same coefficients
and the same LR. On the X the paternal pair must share an allele at every
locus, so a single locus without one excludes it (LR 0, no mutation is
modelled in `kinshipLR`).

## Paternity of a daughter

```jennifer
import "@deckyard/forensicgenetics/kinship.j" as kinship;

def trio as float init xstr.daughterTrioIndex(
    $mother, $daughter, $allegedFather, $model, kinship.NO_MUTATION);
def duo as float init xstr.daughterDuoIndex(
    $daughter, $allegedFather, $model, kinship.DEFAULT_MUTATION);
```

The father's transmission is his one allele, so with a mother typed, a
daughter's paternal allele is fixed wherever the mother does not carry
both of the daughter's alleles, and the index there is `1/p` of that
allele. Without mutation the duo index equals `kinshipLR` for
`FATHER_DAUGHTER`. The mutation model is `kinship.j`'s stepwise model,
applied to the father only.

## Linked loci

X-STRs sit on one chromosome, and kits cluster them. The Investigator Argus
X-12 kit types twelve X-STRs in four linkage groups of three, which the
manufacturer handles as haplotypes:

```jennifer
def groups as map of string to list of string init xstr.argusX12Groups();
# LG1 DXS8378, DXS10135, DXS10148; LG2 DXS7132, DXS10074, DXS10079;
# LG3 HPRTB, DXS10101, DXS10103;   LG4 DXS7423, DXS10134, DXS10146
def lg1 as xstr.LinkageGroup init xstr.linkageGroup("LG1", $groups["LG1"], 0.001, $haplotypes);
def linked as xstr.XFrequencies init xstr.withLinkageGroup($model, $lg1);
```

`$haplotypes` maps `haplotypeKey(alleles)` (alleles in the group's locus
order, joined by `-`) to a frequency from a population study; the deck
bundles none. `minFrequency` is the score for a haplotype the table lacks;
the caller chooses it, no floor is imposed.

Two models, and what each assumes:

| Model | Population | Inheritance |
| --- | --- | --- |
| `independentModel(db)` | loci in linkage equilibrium: haplotype frequency = product of allele frequencies | each locus inherited independently |
| `withLinkageGroup` | haplotype frequencies as tabled (linkage disequilibrium allowed) | a group inherited as one block, no recombination inside it; groups independent |

The independent model needs only allele frequencies, which are widely
published, but misstates both the population (linkage disequilibrium in a
cluster) and the inheritance (a cluster's alleles travel together). The
grouped model needs haplotype frequencies, which need far larger samples: a
three-locus group with ten alleles per locus has a thousand possible
haplotypes, so many observed haplotypes are rare or unseen and the floor
matters. With one meiosis between the two people (father-daughter,
mother-daughter, mother-son) and haplotypes in equilibrium the two models
agree; across several meioses (sisters, grandparents) they do not. In the
tests, two full brothers sharing a two-locus haplotype at `.1 × .1` give
50.5 grouped and 30.25 independent.

A profile typed at part of a group throws: type the whole group, or score
those loci with a model without the group. Female phase is unknown and is
summed over: every pair of haplotypes her genotypes allow, `2^(h-1)` pairs
for `h` heterozygous loci, at most `MAX_GROUP_HETEROZYGOUS_LOCI` (12) per
group.

`examples/xstr.j` runs a half-sister question, a daughter trio and duo, and
the same pair with LG1 as a haplotype, on synthetic numbers. See
[Limitations](../project/limitations.md#linked-x-strs-are-scored-without-recombination).
