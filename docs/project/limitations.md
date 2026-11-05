# Limitations

## Not a validated forensic tool

The deck implements published estimators and is tested against them, but it has
not been through the validation any jurisdiction requires of software used in
casework. Treat it as a library for research, teaching and cross-checking, and
validate independently before it informs a real conclusion.

### Not for legal or court use

**Do not use this deck for a legal proceeding.** It is not intended for, and is
not recommended for, any calculation that informs a prosecution, a defence, a
paternity ruling, an immigration decision, a disaster-victim identification, or
any other legal or official determination. Jurisdictions require validated,
accredited software for that purpose, operated under an accredited process;
this deck is neither. If a number from here would ever reach a court, reproduce
it with validated software and let that result be the one you rely on.

### No warranty, no liability

This software is provided **as is**, without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose, and non-infringement. **You use it entirely
at your own risk.**

The author accepts no responsibility and no liability for any claim, damage,
loss, or other consequence arising from this software or from any number it
produces, however that number was reached and whoever relied on it. Verifying
that a result is correct, appropriate, and lawful for your purpose is your
responsibility alone.

This restates - and does not narrow, extend, or replace - the warranty and
liability disclaimers of the
[LGPL-3.0-only](https://www.gnu.org/licenses/lgpl-3.0.html) licence this deck
is released under. Where the two differ, the licence governs.

Everything below is a specific, technical limitation. None of it narrows the
two statements above.

## Data

### The example frequency table is synthetic

`examples/lib/exampledb.j` builds a table of made-up numbers. Every figure an
example prints illustrates the arithmetic, not a real population. Use a real
table, for example from [STRidER](https://strider.online/), for anything you
intend to rely on.

### Pooling populations is not an uncertain population

`combineFrequencies` merges tables into one. When the inputs' frequencies
differ, the merged table violates Hardy-Weinberg (the Wahlund effect): it
scores homozygous loci too rare, overstating the evidence, and heterozygous
loci too common, understating it, so the error is not a conservative bound in
either direction. For "the donor was German or Austrian" use
`randomMatchProbabilityAcross`. In the test suite's example the two differ by a
third at one homozygous locus, and the gap compounds across a profile. See
[Reference frequency data](../concepts/frequency-data.md#when-the-population-is-uncertain).

### Reading a published table is a separate deck

The core reads no file format.
[`@deckyard/strider`](https://github.com/deckyard/deck-strider-frequency-converter)
converts STRidER-style XML and has its own limitations, chief among them that
it has not been validated against a real STRidER export.

### No reference sequence is bundled

`lineage.align` takes the rCRS (NC_012920.1) region you sequenced as an
argument.

## Statistical scope

### No theta in pedigree calculations

[`pedigree.j`](../casework/pedigrees.md) uses `θ = 0`. The sub-population
correction needs the sampling formula over every founder allele of the
pedigree jointly, which peeling's per-person messages do not carry; applying θ
to founder priors alone would be wrong. Pairwise kinship and trio / duo
paternity with θ are in [`kinship.j`](../casework/kinship.md) (the
`...Theta` functions). See [Theta](../concepts/theta.md).

The original `kinship.j` functions remain at `θ = 0` for compatibility; pass
θ through the `...Theta` variants.

### Kinship is for outbred pairs

The three IBD coefficients describe non-inbred individuals. An inbred pedigree
needs the nine condensed Jacquard coefficients.

### Peeling handles simple pedigrees only

A consanguineous mating, or any other loop in the individual/family graph, is
rejected by `pedigree.validate`. Loops need cutset conditioning, which is not
implemented. A half-sibling pedigree is not a loop.

### Mixtures stop at the semi-continuous LR

`mixture.j` uses allele presence, drop-out and drop-in (LRmix). Peak heights,
stutter and degradation need a continuous model (EuroForMix, STRmix), which is
not implemented. Drop-out is one probability per contributor, the same at
every locus; homozygote drop-out is the `d²` approximation, not the
Balding & Buckleton (2009) two-parameter form. Relatives as alternative
contributors are not modelled. CPI / CPE remain for the
no-assumption question and fail under drop-out. See
[Mixtures](../casework/mixtures.md).

### EMPOP's phylogenetic alignment rules are not applied

`lineage.align` returns a plain optimal alignment with each indel run slid
3' as a block, through homopolymers and tandem repeats alike. EMPOP places
indels where the mtDNA phylogeny prefers, which can differ for a complex
indel. Check such a designation against EMPOP before quoting it.

An insertion before the first reference base is numbered after position
`refStart - 1`; with `refStart = 1` that is position 0 (`0.1G`), which
`parseDifference` accepts for insertions only.

### The default mutation rates are placeholders

`DEFAULT_MUTATION` is `rate = 0.002`, `decay = 0.1`: the right order of
magnitude for paternal autosomal STR mutation, not a published figure. When a
mutation model rescues a trio from a single inconsistency, the CPI is
proportional to the rate at that locus, and observed per-locus rates vary
about tenfold. Set published rates with `kinship.withLocusRate`, or take the
AABB 2003 rates bundled in [`mutation.j`](../casework/mutation-models.md)
(`mutation.aabb2003Rates`, `mutation.kinshipModelWithRates`).

### Mutation models stop at trios and duos

`mutation.j`'s Familias models act in trio and duo paternity only. Pairwise
kinship has no mutation term, and `pedigree.j` has no mutation at all; a
pedigree with an untyped generation between typed relatives (where a
non-stationary model would distort frequencies) is not covered. The bundled
rates are one compilation (AABB 2003, 18 loci); its per-sex rates leave out
mutations of unassigned origin, and no allele- or length-specific rates are
modelled. The PM stationarisation fails for some models (a rare allele, or
Stepwise with microvariants and `rate2` 0); the other Simonsson & Mostad
methods (DP, RM) and Familias' `MaxStabilizedMutrate` bound are not
implemented. Familias' "Custom" matrices cannot be supplied.

### Silent alleles have one frequency per locus

A silent allele is a per-locus frequency the caller supplies; nothing
estimates it, and it does not depend on which allele is silenced or on the
kit. It never mutates. A person showing no peak at a locus (silent/silent)
cannot be entered; leave the locus out of that profile.

### Population statistics assume Hardy-Weinberg and independent loci

`popstats.j`'s frequency statistics (He, MP, PD, PIC, TPI, PE) are expected
values under Hardy-Weinberg proportions, and the combined values multiply over
loci. Neither holds exactly in a structured population. The `observed...`
functions use the sample's genotype shares instead. Frequencies are used
unfloored and rescaled to sum to 1; a partial table (sum below 0.95) is
refused.

### The multi-allelic HWE test is Monte Carlo

`hweExactTest` is exact for two alleles only, and returns 0.0 when the
p-value is below float range (about 1e-308). For STR loci,
`hwePermutationTest` estimates the exact p-value by permutation, with standard
error about `sqrt(p (1 - p) / B)`. Full enumeration and Guo & Thompson's
Markov chain are not implemented. No multiple-testing correction is applied.

### F-statistics are the 1984 estimator only; the LD test is Monte Carlo

`fStatistics` implements Weir & Cockerham (1984) for unlinked autosomal loci;
there is no hierarchical (region / population) analysis, no bootstrap or
jackknife interval over loci, no allele-size statistic (R_ST) and no
alternative multilocus weighting (Weir 1996, Genepop). An individual missing at
a locus is dropped from that locus only. `ldPermutationTest` estimates the
exact p-value by permutation (the probability statistic; Genepop's default G
statistic and its Markov chain are not implemented), within one population
only; there is no combination across populations (Fisher's method) and no
multiple-testing correction in `ldPermutationTests`. Haplotype frequencies and
`D` / `r²` are not estimated.

### X-STRs use theta 0

`xstr.j` has no sub-population correction. Applying the sampling formula to
X-STRs needs the haploid male draws and the haplotype alleles of a linkage
group in one sequence; that extension is not implemented. Report the X-STR
LR as a theta-0 figure.

### Linked X-STRs are scored without recombination

A `LinkageGroup` is inherited as one block, with no recombination between its
loci, and groups are multiplied as if unlinked. Neither holds exactly: the
Argus X-12 groups span a few centimorgans within, and the groups lie on one
chromosome. The error is nil for a single meiosis from a father (his X passes
whole) and grows with the number of female meioses between the two people.
Models with recombination fractions and linkage disequilibrium (Tillmar et al.
2011; FamLinkX, Kling, Dell'Amico & Tillmar 2015) are not implemented. The
independent model ignores linkage disequilibrium and linkage altogether. See
[X-chromosomal STRs](../casework/x-chromosome.md#linked-loci).

### X-STR haplotype tables and floors are the caller's

No X-STR allele or haplotype frequencies are bundled, and `linkageGroup`
imposes no floor: an unlisted haplotype scores at the caller's
`minFrequency`, or throws at 0. Haplotype tables need large samples, and a
rare-haplotype floor can dominate the LR; report the table and the floor.
A female's phase within a group has `2^(h-1)` resolutions for `h`
heterozygous loci; `xstr.MAX_GROUP_HETEROZYGOUS_LOCI` (12) caps `h`, where
one group costs about 0.15 s in a pairwise LR and 0.2 s in a trio index.
Real groups (Argus X-12: three loci) are far below it.

X-STR kinship is pairwise (and trio / duo for a daughter); pedigrees with
more typed people, X-STR mutation in `kinshipLR`, mutation in the mother, and
X-STR power-of-exclusion statistics are not implemented.

### Sampling-uncertainty allowances are single-source only

`reporting.j` applies the size-bias correction, the Bayesian estimator and the
Dirichlet posterior to the single-source RMP and its LR. Mixture, kinship and pedigree LRs get
neither. The posterior treats loci as independent and reads counts as
`p · 2N` from the published frequencies, so a table rounded to few digits
gives slightly wrong counts. Combining the size-bias frequencies with θ
through NRC II 4.2 is the deck's choice; the sources read state the
frequencies only.

### The verbal scale's band edges are a reading

The ENFSI (2015) table gives ranges ("10 - 100") without inclusive ends and
starts its first band at 2. `reporting.verbalEquivalent` closes each band at
the top and treats `(1, 2)` as weak support; `VerbalEquivalent` states the
bounds it used. The guideline presents the scale as an illustration; a
laboratory with its own scale should not use this one.

### Search verdicts are counting rules

`search.search` classes candidates by counts of identical, drop-out and
mismatching loci against thresholds the caller sets. It does not weigh rare
against common alleles, model drop-in, or treat a mismatch at a high-mutation
locus differently. Report a candidate with an RMP or LR over its loci.
Drop-out is recognised only as a query homozygote against a candidate
heterozygote carrying that allele.

### Familial search ranks, it does not identify

`search.familialSearch` uses `kinship.kinshipLR` (theta 0). In a large
database unrelated people reach high LRs by chance, so the ranking is a
shortlist; half siblings, grandparents and avuncular relatives are
indistinguishable. In `examples/search.j` a planted full sibling typed at
eight loci ranks 19th of 200. There is no `log10` variant: a candidate whose
LR passes the float maximum throws a `"search"` error naming it; score that
pair with `kinship.log10KinshipLR`.

### Built-in panels cover a fixed list of kits

`panels.j` ships eleven panels whose locus lists are documented by the
manufacturer or the standard; other kits need `panels.customPanel`.
`normalizeLocus` knows the aliases listed in [panels](../reference/panels.md)
and passes any other name through unchanged.

### Long panels underflow or overflow

`randomMatchProbability` and `pedigreeLikelihood` return plain products. Use
`log10RandomMatchProbability`, `log10PedigreeLikelihood`, or `pedigreeLR`,
which takes the ratio locus by locus. `mixture.mixtureLR` multiplies locus
LRs and throws past the float maximum; `mixture.log10MixtureLR` sums their
logarithms.

A ratio past the float maximum (about `1.8e308`) throws the module's error
instead of returning infinity: `matchLikelihoodRatio` (RMP below `1e-308`; the
log10 LR is `-log10RandomMatchProbability`), `matchLikelihoodRatioAcross`, the
`kinship.kinshipLR...` functions (`log10KinshipLR...`), `pedigree.pedigreeLR`
(`log10PedigreeLR`), `mixture.mixtureLR` (`log10MixtureLR`),
`search.familialSearch` (`kinship.log10KinshipLR` for the pair), and the trio
and duo CPIs, which have no `log10` form:
a CPI that large needs well over 100 loci. `reporting`'s `sizeBiasLikelihoodRatio`,
`bayesianLikelihoodRatio` and `likelihoodRatioLowerBound` throw past `1e300`
(use `log10SizeBiasRandomMatchProbability`, `log10BayesianRandomMatchProbability`
or `log10RmpQuantiles`).

## Performance

### Alignment cost grows with the region, and with divergence

`align` computes a band around the diagonal and widens it (doubling) until
its score beats the best score any path leaving the band can reach, so the
result is the full-matrix optimum. For near-identical sequences (a haplotype
against its own reference) the first band of half-width 16 suffices and the
cost is roughly linear in the region length. Measured with eight
differences: 342 bp in 0.05 s, the ~1.1 kb control region in 0.15 s, 1600 bp
in 0.22 s, 5000 bp (`lineage.MAX_ALIGN_LENGTH`, the longest accepted) in
0.7 s. Divergent sequences widen the band towards the full quadratic
matrix, about 1.5 µs per cell: two unrelated 1100-base sequences take 3.6 s.
A band pass over `lineage.MAX_ALIGN_CELLS` (1,250,000 cells) throws instead,
which refuses unrelated sequences past about 1100 bases and, at 5000 bases,
anything needing a band wider than 124 (roughly 10% substitutions). Sliding
indels 3' is linear: 5000 bases against one take 0.14 s in all.

### Peeling slows as observed alleles accumulate

Pedigree size costs linear time: validation and peeling are iterative and
build their indexes once per call. Measured, one locus with four observed
alleles: a 2000-generation chain (4002 members) validates in 0.2 s and peels in
1.9 s; one couple with 1000 children peels in 0.6 s.

Cost grows steeply with the number of distinct alleles the typed members show
at a locus. Measured, one locus: 4 alleles 11 ms, 8 alleles 90 ms, 12 alleles
0.6 s, 16 alleles 1.9 s, 20 alleles 5.4 s. A large family at a highly
polymorphic locus can take seconds per locus. Lumping already limits this to
the observed alleles; `pedigree.MAX_PEDIGREE_ALLELES` (20) makes a locus with
more throw. At the limit, an untyped couple with ten typed children takes
2.8 s per locus. A typed parent contributes one genotype instead of `G`, so
families with typed parents are cheap: a 24-locus missing-person `pedigreeLR`
(10 members, 8 typed, 15-allele table) takes 0.16 s.

An eight-locus pedigree LR with four typed people takes about 0.1 s.
Peeling memory does not grow with pedigree size.

### Building a database or profile one entry at a time is quadratic

The copy-returning builders copy the whole container before writing, so a loop
of `withLocus` over `n` loci costs `1 + 2 + ... + n` entry copies. Measured,
10 alleles per locus: 25 loci 5.9 ms, 100 loci 78 ms, 400 loci 1077 ms. A real
STR panel of 15 to 45 loci takes tens of milliseconds, and a 44-marker
reference table converts in about 16-21 ms.

`profile(id, genotypes)` builds in one pass; prefer it to folding
`withGenotype` over a list. `FrequencyDb` has no bulk builder: constructing the
structs directly would skip the validation that rejects a misread table, which
is not worth the saving at this size.

### Iterating a list of large values copies each element

`for (def d in $dbs)` binds a copy of each element, and copying a `FrequencyDb`
copies its locus map; inside a loop over loci that is quadratic. Indexing the
list reads the element in place. The same holds for rotating buffers by
assignment and for `lists.push` / `lists.pop`, which return new lists. The
deck's own code indexes instead (`combineFrequencies`, `align`'s two score
rows, the pedigree graph walks).

### Building a large pedigree is quadratic

`withFounder` and `withChild` return copies, so an `n`-member pedigree copies
the member map `n` times: 1.5 s for 1000 members, 23 s for 4000.
`pedigree.pedigreeOf(name, people)` builds in one pass (4000 members in
0.05 s) with the same checks.

### The mixture LR enumerates genotypes

`mixture.j` sums over every genotype combination of the unknowns: about
`(k(k+1)/2)^u / u!` terms for `k` allele classes, 35-45 µs each. Measured
(`jennifer` 0.25.0-dev+43), one locus, Hd with `u` unknowns, θ = 0.01:
6 alleles, `u = 2`, 0.01 s; 10 alleles, `u = 2`, 0.1 s; 6 alleles, `u = 3`,
0.16 s; 8 alleles, `u = 3`, 0.7 s; 10 alleles, `u = 3`, 2.3 s. `locusLR`
evaluates both hypotheses: Hp (one known plus 3 unknowns) and Hd (3 unknowns)
at 10 alleles take 4.3 s per locus. A 24-locus three-person stain, six alleles
per locus, Hp suspect plus 2 unknowns against Hd 3 unknowns, takes 3.7 s
(107184 combinations). `MAX_UNKNOWNS` (3) and `MAX_LOCUS_ALLELES` (10) make
larger loci throw; `MAX_COMBINATIONS` (250000, about 11 s) caps a
whole-mixture call, checked up front. Score larger mixtures in groups of
loci (`restrictMixture`) and multiply the results.

### The permutation HWE test is linear in sample size times permutations

`hwePermutationTest` costs about 7 µs per person per permutation: 500 people
at 1000 permutations take 3.8 s, so `DEFAULT_PERMUTATIONS` (10000) takes about
40 s. Lower the count for screening and raise it for the loci near the
threshold. The LD tests cost about 5 µs per pair, person and permutation:
`ldPermutationTests` on 24 loci (276 pairs) of 2000 people takes about 1 min
at 20 permutations and about 8 h at the default. Memory is
linear in sample size whatever the number of alleles or genotype classes:
count tables up to 65536 cells (256 alleles, or 256 x 256 genotype classes)
are indexed directly, larger ones hashed into at least `2n` slots.

### The posterior bound is Monte Carlo

`reporting.log10RmpQuantiles` draws two or three gamma deviates per locus per
draw, about 16 µs each on the default binary: 10,000 draws over 24 loci take
about 12 s. Its answer carries Monte Carlo error, about
`sqrt(q (1 - q) / draws)` on the probability scale; fix the seed to reproduce
it.

### Mutation matrices are rebuilt per locus and call

Every `mutation.paternityIndexAt` / `duoPaternityIndexAt` builds its
matrices over the locus's whole allele table: `O(n²)` per model interpreted,
plus an `O(n³)` `linalg.solve` with PM. One PM model at a 120-allele locus
takes about 0.19 s, at 200 alleles 0.47 s; a 24-locus trio with 30 alleles
per locus and a PM paternal model 0.47 s. Nothing is cached across calls, so
a batch of alleged fathers rebuilds every matrix.
`mutation.MAX_STATIONARY_ALLELES` (200) caps the detectable alleles under PM.

### Database search cost

Measured on 5000 synthetic profiles x 24 loci, 20 alleles per locus, query
typed at all 24, `FULL_SIBLINGS` (`jennifer` 0.25.0-dev+43): `search.search`
1.55 s, `search.familialSearch` 1.6 s, both linear in the database size.
`familialSearch` memoises `kinship.kinshipLRAt` per locus and candidate
genotype, passing a one-locus table; calling `kinship.kinshipLR` per
candidate took 34 s. On `jennifer-tiny` the same run takes 6.2 s / 2.5 s, and
the cost grows faster than linearly: every allocation pays for garbage
collection over the live database. Both functions rank with `lists.sort` on
plain keys rather than `lists.sortBy`: a function that passes a callback
anywhere in its call tree deep-copies every compound parameter bound from a
top-level variable, which copied the whole database (64 MiB at 2000 x 24) per
call and cost `search` 25.5 s on tiny at 5000.

## Reporting

A number from this deck depends on assumptions it cannot state for you. Report
with it:

- the reference population and the table's source
- θ
- the loci used, especially if a profile was narrowed
- the alternative hypothesis (unrelated, a brother, an unknown man)
- for a mixture LR, both hypotheses' contributors, the drop-out and drop-in
  probabilities
- for a probability of paternity, the prior
- for a lineage marker, the database and its size
- for a sampling-uncertainty bound, the method (size-bias or posterior), the
  prior, the quantile, the number of draws and the seed
- for a verbal equivalent, the scale and its band edges (`VerbalEquivalent`
  holds them)
