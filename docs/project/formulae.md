# Formulae and sources

Every estimator the deck implements, in one place, with where it comes from.
Notation: `pa`, `pb` are allele frequencies after the floor; `θ` is the
co-ancestry coefficient; `N` is a sample or database size.

## Minimum allele frequency

```
pmin = 5 / (2N)
```

Applied on read by `alleleFrequency`, per locus, using that locus's own `N`. An
allele the table never observed is scored at `pmin` rather than zero; an
observed allele rarer than `pmin` is raised to it; the result is capped at 1.

*NRC II (National Research Council, "The Evaluation of Forensic DNA Evidence",
1996), the minimum-allele-frequency recommendation.*

## Single-locus genotype frequency

NRC II recommendation 4.2 - the Balding-Nichols sampling formula, conditioned
on one copy of the genotype already being observed:

```
homozygote  aa:  [2θ + (1-θ)pa] · [3θ + (1-θ)pa] / [(1+θ)(1+2θ)]

heterozygote ab: 2 · [θ + (1-θ)pa] · [θ + (1-θ)pb] / [(1+θ)(1+2θ)]
```

At `θ = 0` these reduce to Hardy-Weinberg: `pa²` and `2 pa pb`.

*NRC II recommendation 4.2; Balding & Nichols (1994).*

## Random match probability and the single-source LR

```
RMP = Π over loci of the genotype frequency
LR  = 1 / RMP
```

`log10 RMP` sums `log10` of the per-locus frequencies instead.

## Mixtures: inclusion and exclusion

```
PI  = (Σ over observed alleles at the locus of p)²      capped at 1
PE  = 1 - PI
CPI = Π over loci of PI
CPE = 1 - CPI
```

No `θ` - PI counts who cannot be excluded rather than estimating a
co-ancestry-corrected match probability.

## Kinship: joint genotype probability

For a non-inbred pair with IBD coefficients `(k0, k1, k2)`:

```
P(Ga, Gb) = k0 · P(Ga) P(Gb)
          + k1 · T1(Ga, Gb)
          + k2 · [Ga == Gb] · P(Ga)
```

where `P(G)` is the Hardy-Weinberg genotype frequency (`θ = 0`) and the IBD-1
term conditions on one of B's alleles being a copy of one of A's:

```
T1(Ga, Gb) = P(Ga) · Σ over i in {Ga.a, Ga.b} of ½ · Q(Gb, i)

Q(G, i) = P(genotype G | one allele fixed at i, the other drawn from the population)
        = pb            if G is (b, b) and i == b
        = pb2           if G is (b1, b2) and i == b1
        = pb1           if G is (b1, b2) and i == b2
        = 0             otherwise
```

`T1` is symmetric in its arguments, and `Σ over all Gb of P(Ga, Gb) = P(Ga)` -
both are asserted in the test overlay.

```
kinship LR = P(Ga, Gb | rel) / P(Ga, Gb | alt)
```

*Standard IBD decomposition; Thompson, "Statistical Inference from Genetic
Data on Pedigrees"; Evett & Weir, "Interpreting DNA Evidence".*

### IBD coefficients

| Relationship | k0 | k1 | k2 |
| --- | --- | --- | --- |
| unrelated | 1 | 0 | 0 |
| parent/child | 0 | 1 | 0 |
| full siblings | 1/4 | 1/2 | 1/4 |
| half siblings, grandparent, avuncular | 1/2 | 1/2 | 0 |
| first cousins | 3/4 | 1/4 | 0 |
| identical | 0 | 0 | 1 |

### Closed forms the tests check

```
parent/child, aa × ab        1 / (2 pa)
full sibs, both aa           (1 + 1/p)² / 4
full sibs, both ab           1/4 + (1 + pa + pb) / (8 pa pb)
identical, genotype g        1 / P(g)
```

## Paternity: the trio index

With `T(G, x)` the Mendelian transmission probability - the count of `x` in `G`
over 2 - and `C = (c1, c2)` the child's genotype:

```
P(C | mother M, father F) = T(M,c1)·F(c2) + T(M,c2)·F(c1)     if c1 ≠ c2
                          = T(M,c1)·F(c1)                      if c1 = c2

X = P(C | M, alleged father)     with F(x) = mutation-aware transmission
Y = P(C | M, random man)         with F(x) = p(x)

PI  = X / Y
CPI = Π over loci of PI
```

`Y = 0` is a maternal exclusion, not an index of zero, and throws.

## Probability of paternity

```
W = CPI / (CPI + 1)                                 prior 0.5
W = CPI·π / (CPI·π + (1 - π))                       prior π
```

A restatement of the CPI under a stated prior, not additional evidence.

## Stepwise mutation model

For a change of `k` repeat steps in a given direction:

```
P(k steps, one direction) = rate · (1 - decay) · decay^(k-1) / 2
P(no change)              = 1 - rate
```

The rates over all non-zero changes sum to `rate`, so the model is a proper
distribution over the allele ladder. A microvariant difference below one full
repeat counts as `k = 1`. Applied to the alleged father's transmission only.

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

Evaluated by **Elston-Stewart peeling** in its message-passing form. The
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
otherwise. Any member may serve as the root `r`; a pedigree in several
components is peeled once per component and the results multiplied.

**Allele lumping.** Alleles no member's genotype shows are pooled into one
class carrying their combined frequency. Because untyped members are
unobserved, those alleles are exchangeable and the pooling is exact.

*Elston & Stewart (1971); the message-passing formulation is equivalent to the
classical anterior/posterior recursion.*

## Lineage markers

```
point estimate   k / N
```

The conservative estimate is the **Clopper-Pearson exact-binomial one-sided
upper bound**: the largest `p` for which observing at most `k` successes in `N`
trials still has probability at least `1 - c`. Equivalently the `c` quantile of
`Beta(k+1, N-k)`. Computed as the upper end of a two-sided Clopper-Pearson
interval at level `2c - 1`:

```jennifer
stats.proportionCi(k, n, 2c - 1, "clopper-pearson").upper
```

Special cases: `1 - α^(1/N)` at `k = 0` (the "rule of three" at 95%), and
exactly 1 at `k = N`.

```
lineage LR = 1 / bound
```

*Clopper & Pearson (1934).*

## Alignment

Global Needleman-Wunsch with match `+1`, mismatch `-1`, gap `-2`; ties in the
traceback prefer the diagonal, so an equal-scoring run is reported as
substitutions rather than as an indel pair. Indels are then slid to the 3' end
of any homopolymer run, which is what makes an extra base in a C-stretch
number from the run's 3' end.

*Needleman & Wunsch (1970). Note that EMPOP's phylogenetic alignment rules are
**not** applied - see [Limitations](limitations.md).*

## Primary sources

- **NRC II** - National Research Council, *The Evaluation of Forensic DNA
  Evidence*, National Academies Press, 1996. Recommendations 4.1 and 4.2, and
  the minimum-allele-frequency rule.
- **STRidER** - the ENFSI STR reference database, <https://strider.online/>,
  and its [formulae page](https://strider.online/formulae), the intended
  validation target for these estimators.
- **EMPOP** - <https://empop.online/>, mtDNA haplotypes and the alignment
  conventions.
- **YHRD** - <https://yhrd.org/>, Y-STR haplotypes.
- Balding & Nichols (1994); Clopper & Pearson (1934); Elston & Stewart (1971);
  Needleman & Wunsch (1970); Evett & Weir, *Interpreting DNA Evidence* (1998);
  Thompson, *Statistical Inference from Genetic Data on Pedigrees* (2000).
