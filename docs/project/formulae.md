# Formulae and sources

Notation: `pa`, `pb` are allele frequencies after the floor; `θ` is the
co-ancestry coefficient; `N` is a sample or database size.

## Minimum allele frequency

```
pmin = 5 / (2N)
```

Applied on read by `alleleFrequency`, per locus, with that locus's own `N`. An
unlisted allele is scored at `pmin`; a listed allele below `pmin` is raised to
it; the result is capped at 1. No floor when `N = 0`.

*NRC II (National Research Council, "The Evaluation of Forensic DNA Evidence",
1996), the minimum-allele-frequency recommendation.*

## Single-locus genotype frequency

NRC II recommendation 4.2, the Balding-Nichols sampling formula, conditioned on
one copy of the genotype already observed:

```
homozygote  aa:  [2θ + (1-θ)pa] · [3θ + (1-θ)pa] / [(1+θ)(1+2θ)]

heterozygote ab: 2 · [θ + (1-θ)pa] · [θ + (1-θ)pb] / [(1+θ)(1+2θ)]
```

At `θ = 0` these reduce to Hardy-Weinberg: `pa²` and `2 pa pb`. `θ` is in
`[0, 1)`; the result is capped at 1.

*NRC II recommendation 4.2; Balding & Nichols (1994).*

## Random match probability and the single-source LR

```
RMP = Π over loci of the genotype frequency
LR  = 1 / RMP
```

`log10 RMP` sums `log10` of the per-locus frequencies.

For an uncertain population with priors `w_i` (non-negative, summing to 1):

```
RMP_across = Σ over populations i of w_i · RMP_i
```

## Mixtures: inclusion and exclusion

```
PI  = (min(1, Σ over observed alleles at the locus of p))²
PE  = 1 - PI
CPI = Π over loci of PI
CPE = 1 - CPI
```

No `θ`: PI counts who cannot be excluded.

## Pooling frequency tables

```
p_pooled(a) = Σ_i p_i(a) · 2N_i / Σ_i 2N_i        N_pooled = Σ_i N_i
```

over the loci every input carries.

## Kinship: joint genotype probability

For a non-inbred pair with IBD coefficients `(k0, k1, k2)` (values in
[Kinship](../casework/kinship.md#ibd-coefficients)):

```
P(Ga, Gb) = k0 · P(Ga) P(Gb)
          + k1 · T1(Ga, Gb)
          + k2 · [Ga == Gb] · P(Ga)
```

`P(G)` is the Hardy-Weinberg genotype frequency (`θ = 0`). The IBD-1 term
conditions on one of B's alleles being a copy of one of A's:

```
T1(Ga, Gb) = P(Ga) · Σ over i in {Ga.a, Ga.b} of ½ · Q(Gb, i)

Q(G, i) = P(genotype G | one allele fixed at i, the other drawn from the population)
        = pb            if G is (b, b) and i == b
        = pb2           if G is (b1, b2) and i == b1
        = pb1           if G is (b1, b2) and i == b2
        = 0             otherwise
```

`T1` is symmetric, and `Σ over all Gb of P(Ga, Gb) = P(Ga)`; the test overlay
asserts both.

```
kinship LR = P(Ga, Gb | rel) / P(Ga, Gb | alt)
```

Closed forms the tests check:

```
parent/child, aa × ab        1 / (2 pa)
full sibs, both aa           (1 + 1/p)² / 4
full sibs, both ab           1/4 + (1 + pa + pb) / (8 pa pb)
identical, genotype g        1 / P(g)
```

*Standard IBD decomposition; Thompson, "Statistical Inference from Genetic
Data on Pedigrees"; Evett & Weir, "Interpreting DNA Evidence".*

## Paternity: the trio index

With `T(G, x)` the Mendelian transmission probability (the count of `x` in `G`,
over 2) and `C = (c1, c2)` the child's genotype:

```
P(C | mother M, father F) = T(M,c1)·F(c2) + T(M,c2)·F(c1)     if c1 ≠ c2
                          = T(M,c1)·F(c1)                      if c1 = c2

X = P(C | M, alleged father)     with F(x) = mutation-aware transmission
Y = P(C | M, random man)         with F(x) = p(x)

PI  = X / Y
CPI = Π over loci of PI
```

`Y = 0` is a maternal exclusion and throws.

## Probability of paternity

```
W = CPI / (CPI + 1)                                 prior 0.5
W = CPI·π / (CPI·π + (1 - π))                       prior π in (0, 1)
```

## Stepwise mutation model

For allele `s` transmitted as `t`, a change of `k` repeat steps in one
direction:

```
P(t | s) = 1 - rate                                   if t = s
         = rate · (1 - decay) · decay^(k-1) / 2       otherwise

k = max(1, round(|value(s) - value(t)|))
F(x) = ½ · P(x | f1) + ½ · P(x | f2)                 alleged father (f1, f2)
```

The non-zero changes in both directions sum to `rate`. A microvariant
difference below one repeat counts as `k = 1`. Applied to the alleged father's
transmission only. `rate` and `decay` are in `[0, 1)`.

## Pedigree likelihood

At one locus, over all genotype assignments `g` to the pedigree's members:

```
L = Σ over g of  Π over founders i of P(gi)
               · Π over non-founders i of P(gi | g_father(i), g_mother(i))
               · Π over i of pen(i, gi)

P(gc | gf, gm) = T(gf, x)·T(gm, y) + T(gf, y)·T(gm, x)   for gc = (x, y), x ≠ y
               = T(gf, x)·T(gm, x)                        for gc = (x, x)

pen(i, g) = 1                    if i is untyped
          = [g == observed(i)]   otherwise
```

Evaluated by Elston-Stewart peeling in its message-passing form. The
individual/family graph of a simple pedigree is a tree; each edge carries a
vector over the joining person's genotypes:

```
m(i → F)(g) = pen(i, g) · prior(i, g) · Π over F' ≠ F incident to i of m(F' → i)(g)

m(F → i)(g), i a child of F, parents f and m:
  = Σ over gf, gm of  m(f → F)(gf) · m(m → F)(gm) · P(g | gf, gm)
                    · Π over siblings s of [ Σ over gs of P(gs | gf, gm) · m(s → F)(gs) ]

m(F → i)(g), i a parent of F with spouse s:
  = Σ over gs of  m(s → F)(gs)
                · Π over children c of [ Σ over gc of P(gc | g, gs) · m(c → F)(gc) ]

L = Σ over g of pen(r, g) · prior(r, g) · Π over F incident to r of m(F → r)(g)
```

`prior(i, g)` is the Hardy-Weinberg genotype frequency for a founder and 1
otherwise. Any member may be the root `r`; a pedigree in several components is
peeled once per component and the results multiplied.

Alleles no typed member shows are pooled into one class with their combined
frequency; untyped members' alleles are exchangeable, so this is exact. If the
observed alleles' floored frequencies sum to 1 or more, they are rescaled to
sum to 1 and there is no lumped class.

*Elston & Stewart (1971); the message-passing form is equivalent to the
classical anterior/posterior recursion.*

## Lineage markers

```
point estimate   k / N
```

The reported estimate is the Clopper-Pearson exact-binomial one-sided upper
bound: the largest `p` for which at most `k` successes in `N` trials still has
probability at least `1 - c`; equivalently the `c` quantile of
`Beta(k+1, N-k)`. Computed as the upper end of the two-sided Clopper-Pearson
interval at level `2c - 1`:

```jennifer
stats.proportionCi(k, n, 2c - 1, "clopper-pearson").upper
```

`1 - α^(1/N)` at `k = 0` (the "rule of three" at 95%), with `α = 1 - c`;
exactly 1 at `k = N`.

```
lineage LR = 1 / bound
```

*Clopper & Pearson (1934).*

## Alignment

Global Needleman-Wunsch with match `+1`, mismatch `-1`, gap `-2`. Ties in the
traceback prefer the diagonal, then a deletion, then an insertion, so an
equal-scoring run is reported as substitutions. Indels are then slid to the 3'
end of any homopolymer run.

*Needleman & Wunsch (1970). EMPOP's phylogenetic alignment rules are not
applied; see [Limitations](limitations.md).*

## Primary sources

- **NRC II**: National Research Council, *The Evaluation of Forensic DNA
  Evidence*, National Academies Press, 1996. Recommendations 4.1 and 4.2, and
  the minimum-allele-frequency rule.
- **STRidER**: the ENFSI STR reference database, <https://strider.online/>,
  and its [formulae page](https://strider.online/formulae), the intended
  validation target for these estimators.
- **EMPOP**: <https://empop.online/>, mtDNA haplotypes and alignment
  conventions.
- **YHRD**: <https://yhrd.org/>, Y-STR haplotypes.
- Balding & Nichols (1994); Clopper & Pearson (1934); Elston & Stewart (1971);
  Needleman & Wunsch (1970); Evett & Weir, *Interpreting DNA Evidence* (1998);
  Thompson, *Statistical Inference from Genetic Data on Pedigrees* (2000).
