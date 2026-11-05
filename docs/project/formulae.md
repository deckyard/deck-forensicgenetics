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

## Mixtures: semi-continuous likelihood ratio

At one locus, `E` the observed alleles, `K` the known contributors with
drop-out `d_k`, `u` unknowns with drop-out `d_u` and genotypes `G`:

```
P(E | H) = Σ over G of P(E | K, G) · P(G | typed)

P(E | K, G) = Π over alleles a carried by someone of
                  (1 - D_a)   if a ∈ E
                  D_a         if a ∉ E
            · Π over a ∈ E carried by nobody of c · p_a      (drop-in)
            · (1 - c)   if no allele dropped in

D_a = Π over contributors i of d_i^(copies of a in i)        homozygote: d²

LR = Π over loci of P(E | Hp) / P(E | Hd)
```

`P(G | typed)` draws the unknowns' alleles one at a time by the
Balding-Nichols sampling formula, after the `n` alleles of the typed profiles
(every known contributor of either hypothesis) and earlier draws:

```
P(next = a) = (n_a θ + (1-θ) p_a) / (1 + (n-1) θ)
```

with a factor 2 per heterozygous unknown. At `θ = 0` this is
`Π P(G_i)` under Hardy-Weinberg. Alleles in neither `E` nor a typed profile
form one class of frequency `1 - Σ p` (exact: none of them is observed).
Unknowns share `d_u`, so combinations are enumerated as multisets weighted
by their number of orderings.

Identity the tests check: one known `S`, `d = 0`, `c = 0`, against one
unknown gives `1 / genotypeFrequency(S, θ)` (NRC II 4.2).

*Curran, Gill & Bill (2005), "Interpretation of repeat measurement DNA
evidence allowing for multiple contributors and population substructure",
Forensic Sci. Int. 148; Gill et al. (2007), "LoComatioN", Forensic Sci. Int.; Haned, Slooten & Gill (2012), "Exploratory data analysis for the
interpretation of low template DNA mixtures", Forensic Sci. Int. Genet. 6;
Gill & Haned (2013), "A new methodological framework to interpret complex DNA
profiles using likelihood ratios", Forensic Sci. Int. Genet. 7; Curran et al.
(1999), "Interpreting DNA mixtures in structured populations", J. Forensic
Sci. 44; Balding & Nichols (1994).*

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

## Kinship and paternity with theta

The `...Theta` functions of `kinship.j` and the duo index. Every allele not
identical by descent is drawn from one sub-population in sequence by the
sampling formula:

```
D(x | S) = (m θ + (1-θ) px) / (1 + (n-1) θ)       n = |S|, m = copies of x in S
Pθ(x1, ..., xr) = Π over i of D(xi | x1, ..., x(i-1))
```

`D(x | ∅) = px`. For genotypes `Ga = (a1, a2)`, `Gb = (b1, b2)` and `o(G)` = 2
for a heterozygote, 1 for a homozygote:

```
P(Ga, Gb) = k0 · o(Ga) o(Gb) Pθ(a1, a2, b1, b2)
          + k1 · o(Ga) Pθ(a1, a2) · Σ over i in {a1, a2} of ½ · Qθ(Gb, i)
          + k2 · [Ga == Gb] · o(Ga) Pθ(a1, a2)

Qθ(G, i) = Q(G, i) with the free allele x scored D(x | a1, a2), not px
```

IBD 1 draws three alleles: B's IBD allele is a copy, not a draw. The terms
are exchangeable, so `P(Ga, Gb)` is symmetric; at `θ = 0` each reduces to
the product-rule term above. Closed forms the tests check, with
`u = θ + (1-θ) pa`, `v = θ + (1-θ) pb`:

```
parent/child, aa × aa        (1 + 2θ) / (3θ + (1-θ) pa)
parent/child, ab × ab        (1 + 2θ)(u + v) / (4 u v)
parent/child, ab × ac        (1 + 2θ) / (4 u)
identical, genotype g        1 / (NRC II 4.2 frequency of g)
full sibs                    1/4 + 1/2 · LR(parent/child) + 1/4 · LR(identical)
```

Trio, mother `M = (m1, m2)`, alleged father `(f1, f2)`: `X` as in the theta-0
index; `Y` with

```
F(x) = D(x | m1, m2, f1, f2) = (m θ + (1-θ) px) / (1 + 3θ)
```

Mother, child and alleged father all `aa`: `PI = (1 + 3θ) / (4θ + (1-θ) pa)`.
Mutation acts only on `X` and theta only on `Y`, so the two combine without
interaction; the random man does not mutate.

Duo, child `C = (c1, c2)`, alleged father `f = (f1, f2)`, `F` the
mutation-aware transmission:

```
X = F(c1) D(c2 | f) + F(c2) D(c1 | f)              c1 ≠ c2
  = F(c1) D(c1 | f)                                c1 = c2
Y = o(C) · D(c1 | f) · D(c2 | f, c1)
```

With no mutation the duo index equals the parent/child kinship LR at the same
`θ`.

*Balding & Nichols (1994) for the sampling formula; Balding & Nichols (1995)
for paternity; Fung & Hu, "Statistical DNA Forensics" (2008), ch. 4
(parentage) and ch. 5 (kinship), and Buckleton, Bright & Taylor, "Forensic
DNA Evidence Interpretation", 2nd ed. (2016), ch. 4 (relatedness) and ch. 11
(parentage), for the kinship and paternity forms with θ.*

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
observed alleles' floored frequencies sum to exactly 1 there is no lumped
class. If they sum past 1 (a small `N` with every allele observed) the locus
throws: kept as they are, an untyped founder's prior would sum to more than 1
and the ratio would depend on how many untyped people each hypothesis names;
rescaled to 1, a floored allele would fall below its `5/(2N)` floor and the
module would stop agreeing with `kinship.j`. Leave the locus out with
`restrict`, or use a larger table.

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

## Population statistics

`popstats.j`. `p_i` are the listed frequencies rescaled to sum to 1 (no
`5/(2N)` floor); `Sk = Σ p_i^k`; `N` individuals, `2N` alleles. Frequency
statistics assume Hardy-Weinberg proportions.

```
He        = 1 - S2                                          Nei (1973)
He_unb    = 2N / (2N - 1) · (1 - S2)                        Nei (1978), eq. 2
MP        = Σ over genotypes of P(g)² = 2 S2² - S4          Fisher (1951); Jones (1972)
PD        = 1 - MP
PIC       = 1 - S2 - Σ_{i<j} 2 p_i² p_j² = 1 - S2 - (S2² - S4)
                                                            Botstein et al. (1980), p. 320
TPI       = 1 / (2 S2)                                      Brenner & Morris (1990), V.C.3
PE_trio   = 1 - 2 S2 + S3 + 2 S4 - 3 S5 - 2 S2² + 3 S2 S3   Jamieson & Taylor (1997), eq. 1a
PE_duo    = 1 - 4 S2 + 2 S2² + 4 S3 - 3 S4                  Jamieson & Taylor (1997), eq. 2a
```

`PE_trio` is the exclusion probability with the mother typed, `PE_duo` with
only child and alleged father typed. Both were checked against direct
enumeration over all mother, child and man genotypes; two alleles at .5 give
3/16 and 1/8.

Across independent loci:

```
MP_comb  = Π MP          log10 MP_comb = Σ log10 MP
PD_comb  = 1 - Π MP
PE_comb  = 1 - Π (1 - PE)              log10 (1 - PE_comb) = Σ log10 (1 - PE)
TPI_comb = Π TPI
```

From a sample of `n` people with `H` the heterozygote share and `h = 1 - H`
(Brenner & Morris use the letters the other way round):

```
p̂_i      = copies of allele i / 2n                          gene counting
MP_obs   = Σ over observed genotypes of x_g²                 Fisher (1951)
TPI_obs  = 1 / (2h)                                          Brenner & Morris (1990), V.C.3
PE_obs   = H² (1 - 2 H h²)                                   Brenner & Morris (1990), V.B.2
```

### Hardy-Weinberg tests

Conditional on the allele counts `n_i`, a genotype table `n_ij` has the
Levene-Haldane probability

```
P(n_ij | n_i) = n! · Π n_i! · 2^H / ((2n)! · Π_{i<=j} n_ij!)
```

with `H` the number of heterozygotes (Levene 1949; Haldane 1954; Guo &
Thompson 1992). The p-value is the total probability of tables no more
likely than the observed one; probabilities within a relative `1e-7` are
ties.

- `hweExactTest`, biallelic: the heterozygote-count distribution from its
  mode by the recurrences of Wigginton et al. (2005), eq. 2:

  ```
  P(h - 2) = P(h) · h (h - 1) / (4 (hom_r + 1)(hom_c + 1))
  P(h + 2) = P(h) · 4 hom_r hom_c / ((h + 2)(h + 1))
  ```

  then normalised.
- `hwePermutationTest`, any number of alleles: the `2n` alleles are shuffled
  (Fisher-Yates) and re-paired `B` times, which samples tables from the
  conditional distribution above (Guo & Thompson 1992, Monte Carlo test). With
  `r` shuffled tables no more likely than the observed one,
  `p = (r + 1) / (B + 1)` (Phipson & Smyth 2010). The generator is Park-Miller,
  `x ← 48271 x mod (2^31 - 1)` (Park & Miller 1988; Park, Miller & Stockmeyer 1993), with rejection so each
  draw is exactly uniform.

### F-statistics

Weir & Cockerham (1984), *Evolution* 38(6):1358-1370. Per allele, `r`
population samples, `n_i` individuals in sample `i`, `p_i` the allele's
frequency and `h_i` the share of individuals heterozygous for it:

```
n̄  = Σ n_i / r                    n_c = (r n̄ - Σ n_i² / (r n̄)) / (r - 1)
p̄  = Σ n_i p_i / (r n̄)            h̄   = Σ n_i h_i / (r n̄)
s² = Σ n_i (p_i - p̄)² / ((r - 1) n̄)
a  = n̄ / n_c · [s² - (p̄ (1 - p̄) - (r - 1) s² / r - h̄ / 4) / (n̄ - 1)]
b  = n̄ / (n̄ - 1) · [p̄ (1 - p̄) - (r - 1) s² / r - (2 n̄ - 1) h̄ / (4 n̄)]
c  = h̄ / 2

θ = Σ a / Σ (a + b + c)     F = Σ (a + b) / Σ (a + b + c)     f = Σ b / Σ (b + c)
```

Sums run over alleles for one locus, over alleles and loci for the multilocus
estimate. With one sample (`r = 1`) the `s²` terms vanish and only `f` is
defined:

```
f = 1 - 2 (n - 1) Ho / (2 n He - Ho)        He = 1 - Σ p̂², Ho observed
```

Verified (5 November 2026) against the Genepop documentation (R package
`genepop` 1.2.17, `inst/doc/GenepopS.tex`, section "ANOVA estimators: single-
and multilocus definitions"): its worked example `Fmulti.txt` (3 samples, 3
loci) gives per-locus and 1984-multilocus FIS, FST and FIT that
`fStatistics` reproduces to the 4 places printed (asserted in
`src/popstats_test.j`), and its single-group form
`FIS = (MSI - MSG) / (MSI + MSG)` equals `b / (b + c)` at `r = 1`. The
components also agree to rounding with hierfstat's ANOVA implementation
(`R/nwc.R`, `sigma² = MSG, (MSI - MSG)/2, (MSP - MSI)/(2 n_c)`) on random
data. Genepop's own multilocus weighting and Weir's (1996) `S1, S2, S3`
weighting differ slightly from the 1984 sums and are not implemented.

### Linkage disequilibrium

Genotypic disequilibrium between two loci, phase unknown (the test of
Genepop's option 2; Raymond & Rousset 1995, *J. Hered.* 86(3):248-249). The
`n` individuals form a contingency table, rows the genotypes at the first
locus, columns those at the second. Under independence, given the margins,

```
P(n_ij | r_i, c_j) = Π r_i! · Π c_j! / (n! · Π n_ij!)
```

`ldPermutationTest` shuffles the second locus's genotypes among individuals
(Fisher-Yates, the Park-Miller generator above), which draws tables from
exactly this distribution. Zaykin, Zhivotovsky & Weir (1995, *Genetica*
96(1-2):169-178) permute alleles and score the conditional probability of the
multi-locus array, which tests within- and between-locus association
together; here single-locus genotypes are permuted whole (conditioning on
them), so only between-locus association is tested. With `r` of `B` shuffled tables no more probable than the
observed one (ties within a relative `1e-7`), `p = (r + 1) / (B + 1)`
(Phipson & Smyth 2010). Genepop evaluates the same conditional distribution
by Markov chain and by default uses the G statistic; this is its probability
test.

## Kit panels

Locus lists only; no estimator. Autosomal STRs per panel, non-autosomal
markers excluded (`panels.Panel.otherMarkers`):

- **CODIS13**: Budowle B, Moretti TR, Niezgoda SJ, Brown BL (1998), "CODIS and
  PCR-based short tandem repeat loci: law enforcement tools", *Proc. 2nd
  European Symposium on Human Identification*, Promega, 73-88.
- **CODIS20**: Hares DR (2015), "Selection and implementation of expanded CODIS
  core loci in the United States", *Forensic Sci. Int. Genet.* 17:33-34.
- **ESS**: Council Resolution 2009/C 296/01 of 30 November 2009, *OJ* C 296,
  Annex 1; Gill P et al. (2006), *Forensic Sci. Int.* 156:242-244.
- **Kits**: the manufacturers' manuals - Thermo Fisher GlobalFiler, NGM SElect
  and Identifiler Plus user guides; Promega PowerPlex Fusion, Fusion 6C and
  ESX 17 technical manuals; QIAGEN Investigator ESSplex SE QS and IDplex Plus
  handbooks.

## Database search

Allele sharing at a locus is identity by state: the size of the multiset
intersection of the two genotypes,

```
IBS(Ga, Gb) = |{a1, a2} ∩ {b1, b2}|  (multisets),  IBS ∈ {0, 1, 2}
```

Direct-search verdicts are counting rules, not statistics (see
[search](../reference/search.md)). A familial search ranks by the pairwise
kinship LR of [Kinship: joint genotype probability](#kinship-joint-genotype-probability),
`rel` against unrelated, evaluated by `kinship.kinshipLRAt` per locus and
multiplied over the loci both profiles carry.

*Bieber FR, Brenner CH, Lazer D (2006), "Finding criminals through DNA of
their relatives", Science 312:1315-1316, for LR-ranked familial searching.*

## X-chromosomal STRs

`xstr.j`, θ = 0. A male has one allele (or haplotype) `x`, a female an
unordered pair `G`. With `P(G)` the Hardy-Weinberg genotype frequency and
`Q(G, i)` as in [Kinship](#kinship-joint-genotype-probability):

```
male match probability        px
female match probability      P(G)

male A, male B                k0 · px py + k1 · [x == y] · px
male A, female B              px · (k0 · P(Gb) + k1 · Q(Gb, x))
female A, female B            k0 · P(Ga) P(Gb) + k1 · T1(Ga, Gb) + k2 · [Ga == Gb] · P(Ga)

LR = P(A, B | rel) / P(A, B | alt)
```

`k2 = 0` for a pair with a male. X coancestry, the probability that an allele
drawn from each is IBD: `k2/2 + k1/4` (two females), `k1` (two males),
`k1/2` (mixed). The relationship constants are derived from X transmission (a
father passes his X to daughters only, whole; a mother passes one of her two)
and checked against the coancestry values of Gomes et al. (2020), Table 5.

Daughter `C = (c1, c2)`, mother `M`, alleged father's allele `f`, `F(x)` the
probability he passes `x` (the stepwise `P(x | f)` of `kinship.j`, 1 or 0
without mutation):

```
trio   X = T(M, c1) F(c2) + T(M, c2) F(c1)        Y = T(M, c1) p(c2) + T(M, c2) p(c1)
duo    X = F(c1) p(c2) + F(c2) p(c1)              Y = P(C)
       (one term each when c1 == c2)              PI = X / Y
```

`T(M, x)` is the mother's Mendelian transmission probability. Without
mutation the duo index equals the father-daughter LR.

Linkage groups: a group's alleles form one haplotype, scored from the
group's table (`max(listed, minFrequency)`). A female's probability at a group
sums over the unordered haplotype pairs her genotypes allow (`2^(h-1)` for `h`
heterozygous loci); a pairwise term sums over both persons' pairs. In the trio
the mother's transmission averages over her pairs weighted by
`P(pair) / Σ P(pair)`. Mutation within a group is the product of per-locus
probabilities. A group is inherited without recombination; groups and
ungrouped loci multiply.

*Pinto N, Gusmão L, Amorim A (2011), "X-chromosome markers in kinship
testing: a generalisation of the IBD approach identifying situations where
their contribution is crucial", Forensic Sci. Int. Genet. 5(1):27-32, for the
X-chromosomal IBD approach; Gomes I, Pinto N, Antão-Sousa S, Gomes V, Gusmão L,
Amorim A (2020), "Twenty years later: a comprehensive review of the X
chromosome use in forensic genetics", Front. Genet. 11:926, Table 5 (X
coancestry by sex), read 2026-11-05. Argus X-12 linkage groups: QIAGEN,
Investigator Argus X-12 QS Handbook, 07/2022, Table 1, read 2026-11-05.
Models with recombination and linkage disequilibrium, not implemented:
Tillmar AO, Egeland T, Lindblom B, Holmlund G, Mostad P (2011), Forensic Sci.
Int. Genet. 5(5):506-511; Kling D, Dell'Amico B, Tillmar AO (2015),
"FamLinkX", Forensic Sci. Int. Genet. 17:1-7.*
## Sampling uncertainty of a match probability

Counts are read back from the table: `x_a = p_a · 2N` (unrounded; 0 for an
unlisted allele). No `5/(2N)` floor.

Size-bias correction: the suspect's and the offender's genotype are added to
the database before the frequency is taken:

```
heterozygote ab:  p'_a = (x_a + 2) / (2N + 4)      p'_b = (x_b + 2) / (2N + 4)
homozygote  aa:   p'_a = (x_a + 4) / (2N + 4)
```

`sizeBiasGenotypeFrequency` evaluates NRC II 4.2 at `p'`; at `θ = 0` that is
`2 p'_a p'_b` or `p'_a²`. The combination with `θ` is the deck's (the
sources read state the frequencies only).

*Size-bias ("pseudocounting") method of Balding & Nichols (1994), "DNA profile
match probability calculation: how to allow for population stratification,
relatedness, database selection and single bands", Forensic Sci. Int.
64:125-140, as stated in Forensic Science Regulator, "Guidance: The use of
allele frequency databases to assign a likelihood ratio (LR) in human DNA
cases", FSR-GUI-0012 Issue 1 (February 2026), section 11.3,
<https://assets.publishing.service.gov.uk/media/5f5253a68fa8f54298e94927/FSR-G-213_DNA_profile_interpretation_Issue_2_Final.pdf>
(read 2026-11-05).
The original paper was not read.*

Posterior bound (the highest-posterior-density approach): at each locus

```
(p_1, ..., p_k, p_rest) ~ Dirichlet(x_1 + α, ..., x_k + α, r)

α = 1 (PRIOR_UNIFORM)  or  α = 1/k (PRIOR_ONE_OVER_K)
k = listed alleles + profile alleles the table lacks
r = 2N - Σ x_i, the count a partial table leaves (no prior; dropped if 0)
```

Only the profile's alleles enter the RMP, so the other classes are summed
into one (the Dirichlet aggregation property; exact). Per draw:
`G_j ~ Gamma(α_j, 1)` independently, `p_j = G_j / Σ G`, the genotype
frequency is NRC II 4.2 at `p`, and the RMP the product over independent
loci. The reported value is the nearest-rank quantile of `draws` simulated
`log10 RMP`; `likelihoodRatioLowerBound` is `10^-(level quantile)`.

Gamma deviates, `α ≥ 1` (Marsaglia & Tsang 2000):

```
d = α - 1/3,  c = 1 / sqrt(9d)
repeat: x ~ N(0, 1), v = (1 + c x)³ (retry if 1 + c x ≤ 0), u ~ U(0, 1)
until u < 1 - 0.0331 x⁴  or  ln u < x²/2 + d (1 - v + ln v)
return d v
```

and for `α < 1`, `Gamma(α + 1) · U^(1/α)`, taken in logarithms so a small `α`
cannot underflow. Normal deviates by the polar method: `v1, v2 ~ U(-1, 1)`,
`s = v1² + v2²` (retry unless `0 < s < 1`), `v1 m` and `v2 m` with
`m = sqrt(-2 ln s / s)`. Uniforms from Park-Miller (multiplier 48271), as
`popstats.hwePermutationTest`.

*Curran, Buckleton, Triggs & Weir (2002), "Assessing uncertainty in DNA
evidence caused by sampling effects", Science & Justice 42(1):29-37 (abstract
read: they prefer the highest-posterior-density method, with non-zero θ);
Curran & Buckleton (2011), "An investigation into the performance of methods
for adjusting for sampling uncertainty in DNA likelihood ratio calculations",
Forensic Sci. Int. Genet. 5(5):512-516 (abstract read: uniform against 1/k
prior); FSR-GUI-0012 section 11.4 (the 1/k prior); Cowell (2015),
arXiv:1509.08361, section 3 (the Dirichlet posterior `Dir(α_i + m_i)` from
database counts). Marsaglia & Tsang (2000), "A simple method for generating
gamma variables", ACM Trans. Math. Softw. 26(3):363-372; Marsaglia & Bray
(1964), "A convenient method for generating normal variables", SIAM Rev.
6(3):260-264. The full texts of Curran et al. (2002), Marsaglia & Tsang and
Marsaglia & Bray were not read; the two sampling algorithms were checked
against the GNU Scientific Library's implementations (`randist/gamma.c`,
which cites Marsaglia & Tsang, and `randist/gauss.c`).*

## Verbal equivalents of a likelihood ratio

| LR | Label (`VerbalEquivalent.label`) |
| --- | --- |
| `= 1` | no support |
| `(1, 10]` | weak support |
| `(10, 100]` | moderate support |
| `(100, 1000]` | moderately strong support |
| `(1000, 10⁴]` | strong support |
| `(10⁴, 10⁶]` | very strong support |
| `> 10⁶` | extremely strong support |

Below 1 the same labels apply to `1/LR` as support for the alternative:
`[0.1, 1)` weak, ..., `< 10⁻⁶` extremely strong. `statement` is the
guideline's first phrasing ("The forensic findings provide moderate support
for the first proposition rather than the alternative."), with the table's
elided "The forensic findings" written out and, below 1, the two propositions
exchanged.

The guideline's table lists the ranges as "2 - 10", "10 - 100", ...,
"1,000,000 and above" without saying which end is inclusive and calls them
indicative; upper-inclusive bands follow the guideline's speaker-recognition
example scale (p. 64: `1 < LR ≤ 10`, `0.1 ≤ LR < 1`, ...), and `(1, 2)` is put
in the weak band. Both are the deck's reading.

*ENFSI, "ENFSI Guideline for Evaluative Reporting in Forensic Science",
approved version 3.0 (2015), guidance note 4, table on p. 17,
<https://enfsi.eu/wp-content/uploads/2016/09/m1_guideline.pdf> (read
2026-11-05). The guideline cites Nordgaard, Ansell, Drotz & Jaeger (2012),
"Scale of conclusions for the value of evidence", Law, Probability and Risk
11:1-24.*
## Familias mutation models

`mutation.j`. At a locus with detectable alleles `1..n`, frequencies `π_i`
(floored table values times `1 - q`) and an optional silent allele `s`
(frequency `q`), a model is a matrix `M[i][j]` = P(`i` transmitted as `j`):

```
Equal          M[i][j] = rate / (n - 1)                     i ≠ j
               M[i][i] = 1 - rate

Proportional   p_j = π_j / S,  S = Σ_detectable π
               alpha = rate / S / Σ_j p_j (1 - p_j)
               M[i][j] = alpha · p_j                        i ≠ j
               M[i][i] = 1 - alpha + alpha · p_i

Stepwise       class(v) = round(10 (v - round(v)))         microvariant class
               C_i = alleles with class(v_j) = class(v_i), j ≠ i
               O_i = alleles of other classes
               M[i][j] = rate · range^k_ij / Σ_{l in C_i} range^k_il    j in C_i, k = |v_i - v_j|
               M[i][j] = rate2 / |O_i|                                  j in O_i
               M[i][i] = 1 - rate          if O_i is empty
                       = 1 - rate2         if C_i is empty
                       = 1 - rate - rate2  otherwise

Silent         M[s][s] = 1,  M[s][j] = M[i][s] = 0
```

Stepwise with `rate2 = 0` is Familias' "Stepwise"; with `rate2 > 0` its
"Extended stepwise". Every row sums to 1; an entry outside `[0, 1]` throws.
Proportional's expected rate `Σ π_i (1 - M[i][i])` is `rate` and it is
stationary.

**PM stationarisation.** With `p` the detectable frequencies normalised to 1,
`v` the stationary distribution of `M` (`v M = v`, `Σ v = 1`, by Gaussian
elimination), `R = 1 - Σ p_i M[i][i]` and `c = 1 - Σ v_i M[i][i]`:

```
d_i = R v_i / (c p_i)
S[i][j] = d_i M[i][j]                  i ≠ j
S[i][i] = 1 - d_i (1 - M[i][i])
```

Then `p S = p` and `Σ p_i (1 - S[i][i]) = R`. It throws when some
`d_i (1 - M[i][i]) > 1` or `v` is not unique.

**Parentage.** For a parent with genotype `(g1, g2)` and matrix `M`,
`T(x) = ½ M[g1][x] + ½ M[g2][x]`. Trio, child `(c1, c2)`:

```
P(C | M, F) = T_M(c1) T_F(c2) + T_M(c2) T_F(c1)        (one term if c1 = c2)
X uses T_F of the alleged father; Y uses T_F(x) = D(x | mother's, alleged father's alleles)
```

`T_M` uses the maternal matrix in both `X` and `Y`; the random man does not
mutate. The duo is the `kinship.j` duo with `T_F` from the paternal matrix.
`D` is the sampling formula of the theta section over the locus frequencies.

*Egeland, Kling & Mostad, "Relationship Inference with Familias and R",
Academic Press (2016); model definitions as implemented in the Familias R
package 2.6.4 (`FamiliasLocus`), cross-checked with pedmut 0.9.1
(`mutationMatrix`); Simonsson & Mostad (2016) for stationary models, PM as in
Familias' `stabilize` and pedmut's `makeStationary`.*

## Silent alleles

A single-peak genotype `a` is `{a/a, a/s}`; a two-peak genotype is itself.
Each likelihood sums over the true genotypes behind the typed ones:

```
P(typed) = Σ over true genotypes G consistent with the typed ones of P(G)
```

with the trio and duo weighting each parental pair by its joint probability
(product rule, or the sampling formula with `θ`), and the pairwise kinship
summing `P(Ga, Gb)` of the kinship section over the true pairs. Closed forms
the tests check (`θ = 0`, `πx = (1-q) px`):

```
duo, child a, alleged father b                 PI = q / ((πa + 2q)(πb + 2q))
trio, mother a/c, child a, alleged father b    PI = q / ((πb + 2q)(πa + q))
Σ over phenotypes B (and B = s/s) of P(A = a, B) = πa² + 2 πa q
```

At `q = 0` every function equals its `kinship.j` original.

*Familias' treatment: a silent allele is the last allele of the frequency
vector, matching a single typed peak with one copy (C++ `pers::execute`, R
`FamiliasLocus`).*

## Published mutation counts

`mutation.aabb2003Counts`: "Apparent Mutations Observed at STR Loci in the
Course of Paternity Testing", data from the AABB 2003 Annual Report (a
compilation of several years, 44 laboratories), as published by NIST STRBase
at <https://strbase.nist.gov/mutation.htm>. Read on 5 November 2026 from the
Internet Archive copy of 18 December 2020
(<https://web.archive.org/web/20201218004320/https://strbase.nist.gov/mutation.htm>);
the live STRBase site no longer serves the page. For each of 18 loci it gives
maternal and paternal mutations / meioses, mutations from either parent, and
totals; every row's parts add up to its totals (asserted). Rates:

```
total     = total mutations / total meioses
paternal  = paternal mutations / paternal meioses     (unassigned excluded)
maternal  = maternal mutations / maternal meioses     (unassigned excluded)
```

Each total rate rounds to STRBase's printed percentage (asserted).

## Primary sources

- **NRC II**: National Research Council, *The Evaluation of Forensic DNA
  Evidence*, National Academies Press, 1996. Recommendations 4.1 and 4.2
  (Equations 4.10 are Balding & Nichols 1994); theta 0.01, or 0.03 for
  isolated populations (ch. 4); rare alleles rebinned to a count of at least
  five (ch. 4, 5), which the deck applies as the 5/(2N) floor.
- **STRidER**: the ENFSI STR reference database, <https://strider.online/>,
  and its [formulae page](https://strider.online/formulae), the intended
  validation target for these estimators.
- **EMPOP**: <https://empop.online/>, mtDNA haplotypes and alignment
  conventions.
- **YHRD**: <https://yhrd.org/>, Y-STR haplotypes.
- **LRmix**: Haned et al. (2012); Gill & Haned (2013); Curran, Gill & Bill
  (2005); Gill et al. (2007).
- Balding & Nichols (1994, 1995); Fung & Hu, *Statistical DNA Forensics*
  (2008); Buckleton, Bright & Taylor, *Forensic DNA Evidence Interpretation*,
  2nd ed. (2016); Clopper & Pearson (1934); Elston & Stewart (1971);
  Needleman & Wunsch (1970); Evett & Weir, *Interpreting DNA Evidence* (1998);
  Thompson, *Statistical Inference from Genetic Data on Pedigrees* (2000).
- Population statistics: Botstein, White, Skolnick & Davis (1980), *Am. J.
  Hum. Genet.* 32:314-331; Brenner & Morris (1990), *Proceedings for the
  International Symposium on Human Identification 1989*, Promega, 21-53; Fisher (1951),
  *Heredity*; Guo & Thompson (1992), *Biometrics*; Haldane (1954), *J.
  Genet.*; Jamieson & Taylor (1997), *Anim. Genet.*; Jones (1972), *J.
  Forensic Sci. Soc.*; Levene (1949), *Ann. Math. Stat.*; Nei (1973), *PNAS*; Nei
  (1978), *Genetics*; Park & Miller (1988), *Commun. ACM* 31:1192-1201;
  Park, Miller & Stockmeyer (1993), *Commun. ACM* 36(7):108-110; Phipson & Smyth
  (2010), *Stat. Appl. Genet. Mol. Biol.*; Wigginton, Cutler & Abecasis
  (2005), *Am. J. Hum. Genet.*
- Population structure and linkage: Weir & Cockerham (1984), *Evolution*
  38(6):1358-1370; Zaykin, Zhivotovsky & Weir (1995), *Genetica*
  96(1-2):169-178; Raymond & Rousset (1995), *J. Hered.* 86(3):248-249; the
  Genepop documentation (Rousset, R package `genepop` 1.2.17, CRAN) for the
  worked F-statistics example and the genotypic-disequilibrium test.
- Reporting: Forensic Science Regulator FSR-GUI-0012 (2026), sections 11.3
  and 11.4; Curran, Buckleton, Triggs & Weir (2002), *Sci. Justice*
  42(1):29-37; Curran & Buckleton (2011), *Forensic Sci. Int. Genet.*
  5(5):512-516; Cowell (2015), arXiv:1509.08361; Marsaglia & Tsang (2000),
  *ACM Trans. Math. Softw.* 26(3):363-372; Marsaglia & Bray (1964), *SIAM
  Rev.* 6(3):260-264; ENFSI Guideline for Evaluative Reporting in Forensic
  Science, approved version 3.0 (2015).
- Mutation models and silent alleles: Egeland, Kling & Mostad, *Relationship
  Inference with Familias and R*, Academic Press (2016), ISBN
  978-0-12-802402-7; the Familias R package 2.6.4 (CRAN) and pedmut 0.9.1
  (CRAN); Simonsson & Mostad (2016), "Stationary mutation models", *Forensic
  Sci. Int. Genet.* 23:217-225, doi:10.1016/j.fsigen.2016.04.005; NIST
  STRBase, "Mutations Observed at STR Loci" (AABB 2003 annual report data).
