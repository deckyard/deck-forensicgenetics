<!-- SPDX-License-Identifier: LGPL-3.0-only -->
<!-- SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev> -->

# forensicgenetics cheatsheet

Every export of `@deckyard/forensicgenetics`. The module pages
([forensicgenetics](forensicgenetics.md), [kinship](kinship.md),
[pedigree](pedigree.md), [lineage](lineage.md), [mixture](mixture.md),
[popstats](popstats.md), [panels](panels.md), [search](search.md),
[xstr](xstr.md), [reporting](reporting.md), [mutation](mutation.md)) list each
module's rules and errors.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/kinship.j" as kinship;
import "@deckyard/forensicgenetics/lineage.j" as lineage;
import "@deckyard/forensicgenetics/mixture.j" as mixture;
import "@deckyard/forensicgenetics/pedigree.j" as pedigree;
import "@deckyard/forensicgenetics/popstats.j" as popstats;
import "@deckyard/forensicgenetics/panels.j" as panels;
import "@deckyard/forensicgenetics/search.j" as search;
import "@deckyard/forensicgenetics/xstr.j" as xstr;
import "@deckyard/forensicgenetics/reporting.j" as reporting;
import "@deckyard/forensicgenetics/mutation.j" as mutation;
```

`tools/check-docs.sh` checks every signature, field and constant value on
this page against the code.

## `forensicgenetics`

| Function | Returns | Does |
|---|---|---|
| `forensicgenetics.alleleValue(allele)` | `float` | Repeat count of an allele label; throws if not numeric. |
| `forensicgenetics.isNumericAllele(allele)` | `bool` | Whether the label is digits with at most one `.`. |
| `forensicgenetics.db(name)` | `FrequencyDb` | Empty database. |
| `forensicgenetics.withLocus(fdb, locus, n, freqs)` | `FrequencyDb` | Copy with one locus added or replaced; validates the table. |
| `forensicgenetics.combineFrequencies(dbs, name)` | `FrequencyDb` | Tables pooled by allele counts over their shared loci. |
| `forensicgenetics.loci(fdb)` | `list of string` | Locus names, insertion order. |
| `forensicgenetics.hasLocus(fdb, locus)` | `bool` | Whether the locus is present. |
| `forensicgenetics.sampleSize(fdb, locus)` | `int` | N at the locus; 0 if unknown. |
| `forensicgenetics.minFrequency(fdb, locus)` | `float` | `5/(2N)`; 0.0 when N is 0. |
| `forensicgenetics.alleleFrequency(fdb, locus, allele)` | `float` | Frequency floored at `5/(2N)`, capped at 1. |
| `forensicgenetics.genotype(locus, a, b)` | `Genotype` | Genotype with the pair in canonical order. |
| `forensicgenetics.isHomozygous(g)` | `bool` | Whether both alleles are equal. |
| `forensicgenetics.alleles(g)` | `list of string` | Distinct alleles: one or two labels. |
| `forensicgenetics.profile(id, genotypes)` | `Profile` | Profile from one genotype per locus. |
| `forensicgenetics.withGenotype(p, g)` | `Profile` | Copy with the genotype set at its locus. |
| `forensicgenetics.profileLoci(p)` | `list of string` | Locus names, insertion order. |
| `forensicgenetics.hasGenotype(p, locus)` | `bool` | Whether the profile is typed at the locus. |
| `forensicgenetics.genotypeAt(p, locus)` | `Genotype` | Genotype at the locus; throws if untyped. |
| `forensicgenetics.commonLoci(fdb, p)` | `list of string` | Loci both cover, in profile order. |
| `forensicgenetics.restrict(p, keep)` | `Profile` | Copy keeping only the named loci. |
| `forensicgenetics.genotypeFrequency(fdb, g, theta)` | `float` | NRC II rec. 4.2 genotype frequency, capped at 1. |
| `forensicgenetics.randomMatchProbability(fdb, p, theta)` | `float` | RMP: product of `genotypeFrequency` over the profile. |
| `forensicgenetics.log10RandomMatchProbability(fdb, p, theta)` | `float` | `log10` RMP, summed per locus. |
| `forensicgenetics.matchLikelihoodRatio(fdb, p, theta)` | `float` | `1 / RMP`. |
| `forensicgenetics.randomMatchProbabilityAcross(dbs, weights, p, theta)` | `float` | `Σ w_i · RMP_i` over candidate populations. |
| `forensicgenetics.matchLikelihoodRatioAcross(dbs, weights, p, theta)` | `float` | `1 / randomMatchProbabilityAcross`. |
| `forensicgenetics.mixture(id, observed)` | `Mixture` | Mixture; repeated alleles collapse. |
| `forensicgenetics.mixtureLoci(m)` | `list of string` | Locus names, insertion order. |
| `forensicgenetics.inclusionProbability(fdb, locus, observed)` | `float` | PI = `(Σ p)²`, sum capped at 1. |
| `forensicgenetics.exclusionProbability(fdb, locus, observed)` | `float` | PE = `1 - PI`. |
| `forensicgenetics.cpi(fdb, m)` | `float` | CPI: product of PI over loci. |
| `forensicgenetics.cpe(fdb, m)` | `float` | CPE = `1 - CPI`. |
| `forensicgenetics.isIncluded(p, m)` | `bool` | Whether `p` is not excluded at any shared locus. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `forensicgenetics.THETA_GENERAL` | `float` | `0.01` | NRC II, broad well-mixed population. |
| `forensicgenetics.THETA_ISOLATED` | `float` | `0.03` | NRC II, small, isolated or endogamous population. |
| `forensicgenetics.THETA_NONE` | `float` | `0.0` | Plain product rule; not for casework. |
| `forensicgenetics.MIN_ALLELE_COUNT` | `float` | `5.0` | The 5 in `5/(2N)`. |
| `forensicgenetics.MAX_FREQUENCY_SUM` | `float` | `1.05` | Largest locus frequency sum `withLocus` accepts. |

### `forensicgenetics.LocusFrequencies`

One locus of a reference table.

| Field | Holds |
|---|---|
| `locus as string` | the locus name, e.g. "D3S1358" |
| `n as int` | individuals typed (2N alleles); 0 if unknown |
| `alleles as map of string to float` | allele label -> published frequency |

### `forensicgenetics.FrequencyDb`

A reference allele-frequency table, keyed by locus.

| Field | Holds |
|---|---|
| `name as string` | the population or dataset label |
| `loci as map of string to LocusFrequencies` | locus name -> its frequencies |

### `forensicgenetics.Genotype`

An unordered allele pair at one autosomal locus; build with `genotype`.

| Field | Holds |
|---|---|
| `locus as string` | the locus name |
| `a as string` | the lower allele in canonical order |
| `b as string` | the higher allele; equal to `a` for a homozygote |

### `forensicgenetics.Profile`

A single-source STR profile: one genotype per locus.

| Field | Holds |
|---|---|
| `id as string` | a caller-chosen label for the profile |
| `loci as map of string to Genotype` | locus name -> the genotype seen there |

### `forensicgenetics.Mixture`

The alleles seen at each locus of a mixed stain.

| Field | Holds |
|---|---|
| `id as string` | a caller-chosen label |
| `loci as map of string to list of string` | locus name -> the alleles observed |

## `kinship`

| Function | Returns | Does |
|---|---|---|
| `kinship.relationship(name, k0, k1, k2)` | `Relationship` | Relationship from `(k0, k1, k2)`; validated. |
| `kinship.transmissionProbability(g, allele)` | `float` | Mendelian transmission: 1.0, 0.5 or 0.0. |
| `kinship.mutationModel(rate, decay)` | `MutationModel` | Stepwise model; rate and decay in `[0, 1)`. |
| `kinship.withLocusRate(model, locus, rate)` | `MutationModel` | Copy with a per-locus rate. |
| `kinship.mutationProbability(model, locus, source, target)` | `float` | P(`source` transmitted as `target`). |
| `kinship.jointGenotypeProbability(db, ga, gb, rel)` | `float` | `k0·P(Ga)P(Gb) + k1·T1 + k2·[Ga==Gb]·P(Ga)`. |
| `kinship.kinshipLRAt(db, ga, gb, rel, alt)` | `float` | One-locus LR, `rel` over `alt`. |
| `kinship.pairLoci(db, a, b)` | `list of string` | Loci both profiles and the table share, in `a`'s order. |
| `kinship.kinshipLR(a, b, rel, db)` | `float` | LR of `rel` against `UNRELATED` over `pairLoci`. |
| `kinship.kinshipLRAgainst(a, b, rel, alt, db)` | `float` | LR of `rel` against `alt` over `pairLoci`. |
| `kinship.paternityIndexAt(db, mother, child, father, model)` | `float` | One-locus PI, `X / Y`. |
| `kinship.trioLoci(db, mother, child, father)` | `list of string` | Loci all three and the table share, in the child's order. |
| `kinship.paternityIndex(mother, child, allegedFather, db)` | `float` | CPI under `DEFAULT_MUTATION`. |
| `kinship.paternityIndexWith(mother, child, allegedFather, db, model)` | `float` | CPI under `model`. |
| `kinship.inconsistentLoci(mother, child, allegedFather, db)` | `list of string` | Loci with PI 0 under `NO_MUTATION`. |
| `kinship.probabilityOfPaternity(cpi)` | `float` | `W = CPI / (CPI + 1)`, prior 0.5. |
| `kinship.probabilityOfPaternityPrior(cpi, prior)` | `float` | `W` under a prior in `(0, 1)`. |
| `kinship.jointGenotypeProbabilityTheta(db, ga, gb, rel, theta)` | `float` | Joint probability, non-IBD alleles by the Balding-Nichols sampling formula. |
| `kinship.kinshipLRAtTheta(db, ga, gb, rel, alt, theta)` | `float` | One-locus LR with theta; equals `kinshipLRAt` at 0. |
| `kinship.kinshipLRTheta(a, b, rel, db, theta)` | `float` | `kinshipLR` with theta. |
| `kinship.kinshipLRAgainstTheta(a, b, rel, alt, db, theta)` | `float` | `kinshipLRAgainst` with theta. |
| `kinship.paternityIndexAtTheta(db, mother, child, father, model, theta)` | `float` | One-locus trio PI; `Y` draws given the mother's and father's alleles. |
| `kinship.paternityIndexTheta(mother, child, allegedFather, db, theta)` | `float` | Trio CPI with theta under `DEFAULT_MUTATION`. |
| `kinship.paternityIndexWithTheta(mother, child, allegedFather, db, model, theta)` | `float` | Trio CPI with theta under `model`. |
| `kinship.duoPaternityIndexAt(db, child, father, model, theta)` | `float` | One-locus duo PI (mother untyped), mutation and theta. |
| `kinship.duoPaternityIndex(child, allegedFather, db, model, theta)` | `float` | Duo CPI over the loci both carry. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `kinship.UNRELATED` | `Relationship` | `Relationship{name: "unrelated", k0: 1.0, k1: 0.0, k2: 0.0}` | (1, 0, 0); the default alternative. |
| `kinship.PARENT_CHILD` | `Relationship` | `Relationship{name: "parent-child", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0). |
| `kinship.FULL_SIBLINGS` | `Relationship` | `Relationship{name: "full siblings", k0: 0.25, k1: 0.5, k2: 0.25}` | (1/4, 1/2, 1/4). |
| `kinship.HALF_SIBLINGS` | `Relationship` | `Relationship{name: "half siblings", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `kinship.GRANDPARENT` | `Relationship` | `Relationship{name: "grandparent", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0), as half siblings. |
| `kinship.AVUNCULAR` | `Relationship` | `Relationship{name: "avuncular", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0), as half siblings. |
| `kinship.FIRST_COUSINS` | `Relationship` | `Relationship{name: "first cousins", k0: 0.75, k1: 0.25, k2: 0.0}` | (3/4, 1/4, 0). |
| `kinship.IDENTICAL` | `Relationship` | `Relationship{name: "identical", k0: 0.0, k1: 0.0, k2: 1.0}` | (0, 0, 1): same person or monozygotic twins. |
| `kinship.NO_MUTATION` | `MutationModel` | `MutationModel{rate: 0.0, decay: 0.0, perLocus: {}}` | Mendelian only. |
| `kinship.DEFAULT_MUTATION` | `MutationModel` | `MutationModel{rate: 0.002, decay: 0.1, perLocus: {}}` | Placeholder, not a published rate. |

### `kinship.Relationship`

A relationship, as its IBD coefficients.

| Field | Holds |
|---|---|
| `name as string` | a label for the relationship |
| `k0 as float` | the probability of sharing no allele identical by descent |
| `k1 as float` | the probability of sharing exactly one |
| `k2 as float` | the probability of sharing both |

### `kinship.MutationModel`

A stepwise STR mutation model. A `k`-step change in one direction has probability `rate · (1 - decay) · decay^(k-1) / 2`.

| Field | Holds |
|---|---|
| `rate as float` | the per-locus, per-meiosis mutation rate |
| `decay as float` | the geometric decay per additional repeat step, in [0, 1) |
| `perLocus as map of string to float` | locus -> rate, overriding `rate` |

## `lineage`

| Function | Returns | Does |
|---|---|---|
| `lineage.haplotypeFrequency(k, n)` | `float` | Point estimate `k / N`; not for reporting. |
| `lineage.haplotypeUpperBound(k, n, confidence)` | `float` | Clopper-Pearson one-sided upper bound. |
| `lineage.lineageMatchProbability(k, n)` | `float` | The bound at `DEFAULT_CONFIDENCE`. |
| `lineage.lineageLikelihoodRatio(k, n)` | `float` | `1 / lineageMatchProbability`. |
| `lineage.parseDifference(text)` | `Difference` | Parse `16189C`, `315.1C`, `249d` and variants. |
| `lineage.formatDifference(d)` | `string` | Designation string. |
| `lineage.mtHaplotype(id, ranges, text)` | `MtHaplotype` | Haplotype from whitespace-separated designations. |
| `lineage.mtDesignations(h)` | `list of string` | Designations in position order. |
| `lineage.mtMatches(a, b)` | `bool` | Same differences; ranges ignored. |
| `lineage.mtDistance(a, b)` | `int` | Size of the symmetric difference. |
| `lineage.mtRangesAgree(a, b)` | `bool` | Same ranges, any order. |
| `lineage.align(reference, sample, refStart)` | `list of Difference` | Needleman-Wunsch differences, indels 3'. |
| `lineage.mtHaplotypeFrom(id, reference, sample, refStart)` | `MtHaplotype` | `align` plus the covered range. |
| `lineage.yHaplotype(id, loci)` | `YHaplotype` | Y haplotype; alleles validated. |
| `lineage.yLoci(h)` | `list of string` | Locus names, insertion order. |
| `lineage.ySharedLoci(a, b)` | `list of string` | Shared loci, in `a`'s order. |
| `lineage.yMatches(a, b)` | `bool` | Exact at every shared locus. |
| `lineage.yStepDistance(a, b)` | `float` | Summed repeat differences over shared loci. |
| `lineage.yMatchesWithin(a, b, steps)` | `bool` | `yStepDistance <= steps`. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `lineage.DEFAULT_CONFIDENCE` | `float` | `0.95` | One-sided level of `lineageMatchProbability`. |
| `lineage.MAX_ALIGN_LENGTH` | `int` | `5000` | Longest sequence `align` accepts, in bases. |

### `lineage.Difference`

One difference from the reference sequence, in rCRS nomenclature.

| Field | Holds |
|---|---|
| `position as int` | the 1-based reference position |
| `insertion as int` | 0 for a substitution or deletion; 1, 2, ... for the n-th base inserted after `position` |
| `base as string` | the observed base; "" for a deletion |
| `kind as string` | "sub", "ins" or "del" |

### `lineage.MtHaplotype`

The sequenced ranges and the differences from the reference within them.

| Field | Holds |
|---|---|
| `id as string` | a caller-chosen label |
| `ranges as list of string` | the sequenced regions, e.g. ["16024-16365", "73-340"] |
| `differences as list of Difference` | the differences from the reference, in ascending position order |

### `lineage.YHaplotype`

A repeat count per locus; multi-copy loci as `"11,14"` or `"11/14"`.

| Field | Holds |
|---|---|
| `id as string` | a caller-chosen label |
| `loci as map of string to string` | locus name -> allele designation(s) |

## `mixture`

| Function | Returns | Does |
|---|---|---|
| `mixture.hypothesis(name, known, knownDropout, unknowns, unknownDropout)` | `Hypothesis` | Validated hypothesis: known contributors with their drop-out, `unknowns` sharing one drop-out. |
| `mixture.conditioningProfiles(hp, hd)` | `list of forensics.Profile` | Known contributors of both hypotheses, deduplicated; the θ conditioning set. |
| `mixture.minimumContributors(m)` | `int` | `ceil(max alleles at a locus / 2)`. |
| `mixture.restrictMixture(m, keep)` | `forensics.Mixture` | Copy keeping only the named loci. |
| `mixture.evidenceProbability(fdb, m, locus, h, typed, dropIn, theta)` | `float` | Probability of the evidence under `h` at one locus, θ conditioned on `typed`. |
| `mixture.locusLR(fdb, m, locus, hp, hd, dropIn, theta)` | `float` | LR of `hp` against `hd` at one locus. |
| `mixture.locusLRs(fdb, m, hp, hd, dropIn, theta)` | `map of string to float` | `locusLR` at every locus of the mixture. |
| `mixture.mixtureLR(fdb, m, hp, hd, dropIn, theta)` | `float` | Product of the locus LRs. |
| `mixture.log10MixtureLR(fdb, m, hp, hd, dropIn, theta)` | `float` | `log10` LR, summed per locus. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `mixture.MAX_UNKNOWNS` | `int` | `3` | Most unknown contributors per hypothesis. |
| `mixture.MAX_LOCUS_ALLELES` | `int` | `10` | Most distinct alleles (stain plus typed profiles) per locus. |

### `mixture.Hypothesis`

Who contributed to the stain, under one hypothesis. Build with `hypothesis`.

| Field | Holds |
|---|---|
| `name as string` | a label, used in error messages |
| `known as list of forensics.Profile` | the known contributors |
| `knownDropout as list of float` | drop-out probability per known contributor, in `known` order |
| `unknowns as int` | number of unknown contributors, 0 to `MAX_UNKNOWNS` |
| `unknownDropout as float` | drop-out probability shared by every unknown |

## `pedigree`

| Function | Returns | Does |
|---|---|---|
| `pedigree.pedigree(name)` | `Pedigree` | Empty pedigree. |
| `pedigree.withFounder(ped, id, sex)` | `Pedigree` | Copy with a founder added. |
| `pedigree.withChild(ped, id, sex, father, mother)` | `Pedigree` | Copy with a child of `father` and `mother`. |
| `pedigree.members(ped)` | `list of string` | Member ids, insertion order. |
| `pedigree.hasMember(ped, id)` | `bool` | Whether the id is a member. |
| `pedigree.personOf(ped, id)` | `Person` | The member; throws if absent. |
| `pedigree.isFounder(ped, id)` | `bool` | Whether the member has no parents in the pedigree. |
| `pedigree.childrenOf(ped, id)` | `list of string` | Children's ids, member order. |
| `pedigree.validate(ped)` | - | Throws on a missing parent, an ancestor cycle or a loop. |
| `pedigree.isSimple(ped)` | `bool` | Non-throwing `validate`. |
| `pedigree.locusLikelihood(ped, genotypes, db, locus)` | `float` | Likelihood at one locus, untyped members summed out. |
| `pedigree.pedigreeLoci(ped, profiles, db)` | `list of string` | Loci the table covers and a member is typed at. |
| `pedigree.pedigreeLikelihood(ped, profiles, db)` | `float` | Product of `locusLikelihood` over `pedigreeLoci`. |
| `pedigree.log10PedigreeLikelihood(ped, profiles, db)` | `float` | `log10` likelihood, summed per locus. |
| `pedigree.pedigreeLR(numerator, denominator, profiles, db)` | `float` | LR of two pedigrees, locus by locus. |
| `pedigree.typedMembers(ped, profiles)` | `list of string` | Members with a non-empty profile. |

### `pedigree.Person`

One member. A founder has `""` for both parents.

| Field | Holds |
|---|---|
| `id as string` | the member's identifier, unique within the pedigree |
| `sex as string` | "m", "f" or ""; not used by the calculation |
| `father as string` | the father's id, or "" for a founder |
| `mother as string` | the mother's id, or "" for a founder |

### `pedigree.Pedigree`

A pedigree: a set of members and the parentages joining them.

| Field | Holds |
|---|---|
| `name as string` | a label for the pedigree, used in error messages |
| `members as map of string to Person` | id -> member |

## `popstats`

| Function | Returns | Does |
|---|---|---|
| `popstats.expectedHeterozygosity(fdb, locus)` | `float` | `1 - Σ p²`. |
| `popstats.unbiasedHeterozygosity(fdb, locus)` | `float` | `2N / (2N - 1) · (1 - Σ p²)`; needs `n > 0`. |
| `popstats.matchProbability(fdb, locus)` | `float` | `Σ P(g)² = 2 (Σ p²)² - Σ p⁴` under HWE. |
| `popstats.powerOfDiscrimination(fdb, locus)` | `float` | `1 - matchProbability`. |
| `popstats.polymorphismInformationContent(fdb, locus)` | `float` | Botstein et al. (1980) PIC. |
| `popstats.typicalPaternityIndex(fdb, locus)` | `float` | `1 / (2 Σ p²)`. |
| `popstats.powerOfExclusionTrio(fdb, locus)` | `float` | Jamieson & Taylor (1997), mother typed. |
| `popstats.powerOfExclusionDuo(fdb, locus)` | `float` | Jamieson & Taylor (1997), mother untyped. |
| `popstats.locusStats(fdb, locus)` | `LocusStats` | Every frequency statistic of one locus. |
| `popstats.statsTable(fdb)` | `list of LocusStats` | `locusStats` per locus, table order. |
| `popstats.combinedStats(table)` | `CombinedStats` | Products over loci, with log10 forms. |
| `popstats.combinedPower(powers)` | `float` | `1 - Π (1 - x)`. |
| `popstats.log10CombinedComplement(powers)` | `float` | `Σ log10 (1 - x)`. |
| `popstats.observedHeterozygosity(genotypes)` | `float` | Heterozygotes / n. |
| `popstats.alleleCounts(genotypes)` | `map of string to int` | Allele copies; sums to 2n. |
| `popstats.estimateFrequencies(genotypes)` | `map of string to float` | Gene counting: copies / 2n. |
| `popstats.withEstimatedLocus(fdb, genotypes)` | `FrequencyDb` | Copy with the sample's counted locus, `n` = sample size. |
| `popstats.observedMatchProbability(genotypes)` | `float` | `Σ x_g²` over observed genotype shares. |
| `popstats.observedTypicalPaternityIndex(genotypes)` | `float` | `1 / (2h)`, `h` observed homozygosity. |
| `popstats.observedPowerOfExclusion(genotypes)` | `float` | Brenner & Morris (1990) `H² (1 - 2 H h²)`. |
| `popstats.hweExactTest(genotypes)` | `float` | Biallelic exact HWE p-value (Wigginton et al. 2005). |
| `popstats.hwePermutationTest(genotypes, permutations, seed)` | `float` | Seeded permutation HWE p-value, any allele count. |
| `popstats.fStatistics(samples)` | `FStatistics` | Weir & Cockerham (1984) `a`, `b`, `c` of one locus over population samples. |
| `popstats.fStatisticsTable(populations, loci)` | `list of FStatistics` | `fStatistics` per locus from profiles grouped by population. |
| `popstats.combineFStatistics(table)` | `FStatistics` | Components summed over loci (multilocus estimate). |
| `popstats.fst(stats)` | `float` | θ = `a / (a + b + c)`; needs two populations. |
| `popstats.fit(stats)` | `float` | F = `(a + b) / (a + b + c)`; needs two populations. |
| `popstats.fis(stats)` | `float` | f = `b / (b + c)`. |
| `popstats.inbreedingCoefficient(genotypes)` | `float` | f of one sample (r = 1). |
| `popstats.ldPermutationTest(first, second, permutations, seed)` | `float` | Seeded permutation test of genotypic LD between two loci. |
| `popstats.ldPermutationTests(profiles, permutations, seed)` | `list of LdTest` | `ldPermutationTest` for every locus pair of a profile set. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `popstats.MIN_FREQUENCY_SUM` | `float` | `0.95` | Smallest locus frequency sum the statistics accept. |
| `popstats.DEFAULT_PERMUTATIONS` | `int` | `10000` | Permutation count; SE about 0.002 at p = 0.05. |
| `popstats.MAX_PERMUTATIONS` | `int` | `1000000` | Largest count `hwePermutationTest` and the LD tests accept. |

### `popstats.LocusStats`

Forensic-efficiency parameters of one locus, from its frequencies under HWE.

| Field | Holds |
|---|---|
| `locus as string` | the locus name |
| `n as int` | individuals typed, as the table states; 0 if unknown |
| `alleleCount as int` | the number of alleles the table lists |
| `expectedHeterozygosity as float` | `1 - Σ p²` |
| `matchProbability as float` | `Σ P(g)²` over every genotype g |
| `powerOfDiscrimination as float` | `1 - matchProbability` |
| `pic as float` | the polymorphism information content |
| `typicalPaternityIndex as float` | `1 / (2 Σ p²)` |
| `powerOfExclusionTrio as float` | the share of random men a mother-child-man trio excludes |
| `powerOfExclusionDuo as float` | the share of random men a child-man duo excludes |

### `popstats.CombinedStats`

Per-locus statistics combined over independent loci.

| Field | Holds |
|---|---|
| `loci as int` | the number of loci combined |
| `matchProbability as float` | `Π MP` |
| `log10MatchProbability as float` | `Σ log10 MP` |
| `powerOfDiscrimination as float` | `1 - Π MP` |
| `powerOfExclusionTrio as float` | `1 - Π (1 - PE_trio)` |
| `log10NonExclusionTrio as float` | `Σ log10 (1 - PE_trio)` |
| `powerOfExclusionDuo as float` | `1 - Π (1 - PE_duo)` |
| `log10NonExclusionDuo as float` | `Σ log10 (1 - PE_duo)` |
| `typicalPaternityIndex as float` | `Π TPI` |
| `log10TypicalPaternityIndex as float` | `Σ log10 TPI` |

### `popstats.FStatistics`

Weir & Cockerham (1984) variance components, summed over alleles (and loci).

| Field | Holds |
|---|---|
| `locus as string` | the locus name; "" for a combination of loci |
| `loci as int` | the number of loci summed |
| `populations as int` | the number of population samples, r |
| `a as float` | the between-population component |
| `b as float` | the between-individual, within-population component |
| `c as float` | the within-individual component |

### `popstats.LdTest`

One locus pair from `ldPermutationTests`.

| Field | Holds |
|---|---|
| `first as string` | the first locus, in name order |
| `second as string` | the second locus |
| `n as int` | individuals typed at both |
| `pValue as float` | the permutation p-value; 1.0 when no table can be formed |

## `panels`

| Function | Returns | Does |
|---|---|---|
| `panels.panelNames()` | `list of string` | Names of the built-in panels. |
| `panels.hasPanel(name)` | `bool` | Whether a built-in panel has the name; case and separators ignored. |
| `panels.panel(name)` | `Panel` | Built-in panel by name; throws if unknown. |
| `panels.customPanel(name, loci)` | `Panel` | Caller-defined panel, loci normalised; duplicates throw. |
| `panels.sharedLoci(a, b)` | `list of string` | Loci both panels type, in `a`'s order. |
| `panels.intersect(a, b)` | `Panel` | `sharedLoci` as a panel `"<a> & <b>"`; throws if disjoint. |
| `panels.normalizeLocus(name)` | `string` | Canonical locus name (`THO1` -> `TH01`, `PENTA_E` -> `Penta E`). |
| `panels.normalizeProfile(p)` | `Profile` | Copy with every locus renamed; collisions throw. |
| `panels.normalizeDb(fdb)` | `FrequencyDb` | Copy with every locus renamed; collisions throw. |
| `panels.restrictProfile(p, pan)` | `Profile` | Panel loci the profile has, panel order, canonical names. |
| `panels.restrictDb(fdb, pan)` | `FrequencyDb` | Panel loci the table has, panel order, canonical names. |
| `panels.missingInProfile(p, pan)` | `list of string` | Panel loci the profile is not typed at. |
| `panels.missingInDb(fdb, pan)` | `list of string` | Panel loci the table does not cover. |

Built-in panels (autosomal loci only): `CODIS13` (13), `CODIS20` (20), `ESS`
(12), `GlobalFiler` (21), `PowerPlex Fusion` (22), `PowerPlex Fusion 6C` (23),
`PowerPlex ESX 17`, `NGM SElect`, `Investigator ESSplex SE QS` (16 each),
`Identifiler Plus`, `Investigator IDplex Plus` (15 each).

### `panels.Panel`

A named set of autosomal STR loci.

| Field | Holds |
|---|---|
| `name as string` | the panel's name, e.g. "GlobalFiler" |
| `source as string` | the document the locus list is taken from |
| `loci as list of string` | the autosomal STR loci, normalised names |
| `otherMarkers as list of string` | non-autosomal markers the kit also types; never used in a calculation |

## `search`

| Function | Returns | Does |
|---|---|---|
| `search.searchOptions(minLoci, maxMismatches, allowDropout)` | `SearchOptions` | Validated options; `minLoci >= 1`, `maxMismatches >= 0`. |
| `search.allelesShared(ga, gb)` | `int` | IBS 0, 1 or 2: multiset intersection of the allele pairs. |
| `search.compareProfiles(query, candidate, opts)` | `Comparison` | One candidate: loci compared, identical, drop-out, mismatches, IBS, verdict. |
| `search.search(query, database, opts)` | `list of Comparison` | Every candidate, ranked: verdict, mismatches, identical loci, shared alleles. |
| `search.hits(results)` | `list of Comparison` | The exact and partial candidates, order kept. |
| `search.familialSearch(query, database, rel, fdb, minLoci)` | `list of FamilialHit` | Candidates ranked by `kinship.kinshipLR` for `rel`; unscored last. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `search.VERDICT_EXACT` | `string` | `"exact"` | Every compared locus identical. |
| `search.VERDICT_PARTIAL` | `string` | `"partial"` | At most `maxMismatches` mismatching loci. |
| `search.VERDICT_MISMATCH` | `string` | `"mismatch"` | More mismatching loci than allowed. |
| `search.VERDICT_INSUFFICIENT` | `string` | `"insufficient"` | Fewer than `minLoci` loci compared. |

### `search.SearchOptions`

How a direct search classes a candidate.

| Field | Holds |
|---|---|
| `minLoci as int` | fewest loci compared for a verdict other than "insufficient" |
| `maxMismatches as int` | most mismatching loci a "partial" may have |
| `allowDropout as bool` | a homozygous query locus whose allele the candidate carries is drop-out, not a mismatch |

### `search.Comparison`

One candidate of a direct search.

| Field | Holds |
|---|---|
| `id as string` | the candidate's id |
| `index as int` | the candidate's position in the database list |
| `lociCompared as int` | loci both profiles are typed at |
| `matchingLoci as int` | compared loci with identical genotypes |
| `dropoutLoci as int` | compared loci counted as drop-out |
| `mismatchedLoci as list of string` | the mismatching loci, in query order |
| `ibs as map of string to int` | compared locus -> alleles shared (0, 1, 2) |
| `ibs0 as int` | compared loci sharing no allele |
| `ibs1 as int` | compared loci sharing one allele |
| `ibs2 as int` | compared loci sharing both |
| `sharedAlleles as int` | the sum of `ibs` |
| `verdict as string` | "exact", "partial", "mismatch" or "insufficient" |

### `search.FamilialHit`

One candidate of a familial search.

| Field | Holds |
|---|---|
| `id as string` | the candidate's id |
| `index as int` | the candidate's position in the database list |
| `lociCompared as int` | loci both profiles are typed at |
| `scored as bool` | false when fewer than `minLoci` loci were compared |
| `lr as float` | the kinship LR against unrelated; 0.0 when unscored |

## `xstr`

| Function | Returns | Does |
|---|---|---|
| `xstr.xRelationship(name, sexA, sexB, k0, k1, k2)` | `XRelationship` | Relationship from X-IBD coefficients and both sexes; validated. |
| `xstr.unrelated(sexA, sexB)` | `XRelationship` | (1, 0, 0) for the pair of sexes. |
| `xstr.coancestry(rel)` | `float` | X coancestry: `k2/2 + k1/4` (F-F), `k1` (M-M), `k1/2` (mixed). |
| `xstr.maleProfile(id, alleles)` | `XProfile` | Male profile from locus -> allele. |
| `xstr.femaleProfile(id, genotypes)` | `XProfile` | Female profile from one genotype per locus. |
| `xstr.haplotypeKey(alleles)` | `string` | Alleles joined by `HAPLOTYPE_SEPARATOR`. |
| `xstr.linkageGroup(name, loci, minFrequency, haplotypes)` | `LinkageGroup` | Linked loci with a haplotype frequency table; validated. |
| `xstr.haplotypeFrequency(group, alleles)` | `float` | Listed frequency floored at `minFrequency`. |
| `xstr.argusX12Groups()` | `map of string to list of string` | Argus X-12 linkage groups LG1-LG4 (QIAGEN handbook, Table 1). |
| `xstr.independentModel(db)` | `XFrequencies` | Every locus scored on its own from `db`. |
| `xstr.withLinkageGroup(xf, group)` | `XFrequencies` | Model with the group's loci scored as one haplotype. |
| `xstr.matchProbability(xf, p)` | `float` | Male `p`, female HW genotype frequency (phase summed), product over loci and groups. |
| `xstr.log10MatchProbability(xf, p)` | `float` | `log10` of `matchProbability`, summed. |
| `xstr.pairUnits(xf, a, b)` | `list of string` | Loci and group names a pair is scored over. |
| `xstr.kinshipLR(a, b, rel, xf)` | `float` | X kinship LR of `rel` against unrelated. |
| `xstr.kinshipLRAgainst(a, b, rel, alt, xf)` | `float` | X kinship LR of `rel` against `alt`. |
| `xstr.log10KinshipLRAgainst(a, b, rel, alt, xf)` | `float` | `log10` of `kinshipLRAgainst`, summed; throws on an exclusion. |
| `xstr.daughterTrioIndex(mother, daughter, allegedFather, xf, model)` | `float` | Paternity index for a daughter, mother typed. |
| `xstr.daughterDuoIndex(daughter, allegedFather, xf, model)` | `float` | Paternity index for a daughter, mother untyped. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `xstr.MALE` | `string` | `"male"` | One X allele per locus. |
| `xstr.FEMALE` | `string` | `"female"` | Two X alleles per locus. |
| `xstr.HAPLOTYPE_SEPARATOR` | `string` | `"-"` | Joins a haplotype key. |
| `xstr.FATHER_DAUGHTER` | `XRelationship` | `XRelationship{name: "father-daughter", sexA: "male", sexB: "female", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0). |
| `xstr.FATHER_SON` | `XRelationship` | `XRelationship{name: "father-son", sexA: "male", sexB: "male", k0: 1.0, k1: 0.0, k2: 0.0}` | (1, 0, 0): no X shared. |
| `xstr.MOTHER_SON` | `XRelationship` | `XRelationship{name: "mother-son", sexA: "female", sexB: "male", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0). |
| `xstr.MOTHER_DAUGHTER` | `XRelationship` | `XRelationship{name: "mother-daughter", sexA: "female", sexB: "female", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0). |
| `xstr.FULL_SISTERS` | `XRelationship` | `XRelationship{name: "full sisters", sexA: "female", sexB: "female", k0: 0.0, k1: 0.5, k2: 0.5}` | (0, 1/2, 1/2). |
| `xstr.FULL_BROTHERS` | `XRelationship` | `XRelationship{name: "full brothers", sexA: "male", sexB: "male", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.BROTHER_SISTER` | `XRelationship` | `XRelationship{name: "brother-sister", sexA: "male", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.PATERNAL_HALF_SISTERS` | `XRelationship` | `XRelationship{name: "paternal half sisters", sexA: "female", sexB: "female", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0), as mother-daughter. |
| `xstr.MATERNAL_HALF_SISTERS` | `XRelationship` | `XRelationship{name: "maternal half sisters", sexA: "female", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.MATERNAL_HALF_BROTHERS` | `XRelationship` | `XRelationship{name: "maternal half brothers", sexA: "male", sexB: "male", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.PATERNAL_GRANDMOTHER_GRANDDAUGHTER` | `XRelationship` | `XRelationship{name: "paternal grandmother-granddaughter", sexA: "female", sexB: "female", k0: 0.0, k1: 1.0, k2: 0.0}` | (0, 1, 0). |
| `xstr.MATERNAL_GRANDMOTHER_GRANDDAUGHTER` | `XRelationship` | `XRelationship{name: "maternal grandmother-granddaughter", sexA: "female", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.MATERNAL_GRANDFATHER_GRANDDAUGHTER` | `XRelationship` | `XRelationship{name: "maternal grandfather-granddaughter", sexA: "male", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.PATERNAL_AUNT_NIECE` | `XRelationship` | `XRelationship{name: "paternal aunt-niece", sexA: "female", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.PATERNAL_UNCLE_NIECE` | `XRelationship` | `XRelationship{name: "paternal uncle-niece", sexA: "male", sexB: "female", k0: 0.5, k1: 0.5, k2: 0.0}` | (1/2, 1/2, 0). |
| `xstr.MATERNAL_AUNT_NIECE` | `XRelationship` | `XRelationship{name: "maternal aunt-niece", sexA: "female", sexB: "female", k0: 0.25, k1: 0.75, k2: 0.0}` | (1/4, 3/4, 0). |

### `xstr.XProfile`

An X-STR profile; a male genotype has `a == b`.

| Field | Holds |
|---|---|
| `id as string` | a caller-chosen label |
| `sex as string` | `MALE` or `FEMALE` |
| `loci as map of string to Genotype` | locus name -> the genotype seen there |

### `xstr.XRelationship`

X-IBD coefficients for an ordered pair of sexes.

| Field | Holds |
|---|---|
| `name as string` | a label, read A to B ("father-daughter") |
| `sexA as string` | the first person's sex, `MALE` or `FEMALE` |
| `sexB as string` | the second person's sex |
| `k0 as float` | the probability that no X allele is shared IBD |
| `k1 as float` | the probability that exactly one is |
| `k2 as float` | the probability that two are (two females only) |

### `xstr.LinkageGroup`

Haplotype frequencies for tightly linked loci.

| Field | Holds |
|---|---|
| `name as string` | the group's label, e.g. "LG1" |
| `loci as list of string` | the linked loci, in key order |
| `minFrequency as float` | the floor for rare and unseen haplotypes, in [0, 1) |
| `haplotypes as map of string to float` | haplotype key -> frequency |

### `xstr.XFrequencies`

The X-STR frequency model.

| Field | Holds |
|---|---|
| `db as FrequencyDb` | single-locus X-STR allele frequencies |
| `groups as list of LinkageGroup` | the linkage groups, disjoint |
## `reporting`

| Function | Returns | Does |
|---|---|---|
| `reporting.sizeBiasAlleleFrequency(fdb, g, allele)` | `float` | `(x + c) / (2N + 4)`, `c` the copies two carriers of `g` add (2 or 4); no floor. |
| `reporting.sizeBiasGenotypeFrequency(fdb, g, theta)` | `float` | NRC II 4.2 at size-bias corrected frequencies. |
| `reporting.sizeBiasRandomMatchProbability(fdb, p, theta)` | `float` | Product of `sizeBiasGenotypeFrequency` over the profile. |
| `reporting.log10SizeBiasRandomMatchProbability(fdb, p, theta)` | `float` | Its log10, summed per locus. |
| `reporting.sizeBiasLikelihoodRatio(fdb, p, theta)` | `float` | `1 / sizeBiasRandomMatchProbability`. |
| `reporting.log10RmpQuantiles(fdb, p, theta, prior, quantiles, draws, seed)` | `list of float` | Posterior quantiles of log10 RMP under per-locus Dirichlet frequencies; seeded Monte Carlo. |
| `reporting.likelihoodRatioLowerBound(fdb, p, theta, prior, level, draws, seed)` | `float` | `1 / RMP` at the `level` posterior quantile of the RMP. |
| `reporting.verbalEquivalent(lr)` | `VerbalEquivalent` | ENFSI (2015) verbal band, wording and bounds; `lr > 0`. |
| `reporting.verbalEquivalentLog10(log10Lr)` | `VerbalEquivalent` | The same from log10 LR. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `reporting.PRIOR_UNIFORM` | `int` | `1` | Dirichlet prior `α = 1` per allele. |
| `reporting.PRIOR_ONE_OVER_K` | `int` | `2` | Dirichlet prior `α = 1/k` per allele. |
| `reporting.DEFAULT_DRAWS` | `int` | `10000` | Simulation size; quantile SE about `sqrt(q (1 - q) / draws)`. |
| `reporting.MAX_DRAWS` | `int` | `1000000` | Largest simulation size accepted. |
| `reporting.SUPPORTS_FIRST` | `string` | `"first"` | LR above 1. |
| `reporting.SUPPORTS_ALTERNATIVE` | `string` | `"alternative"` | LR below 1. |
| `reporting.SUPPORTS_NEITHER` | `string` | `"neither"` | LR exactly 1. |

### `reporting.VerbalEquivalent`

A likelihood ratio's band on the ENFSI (2015) scale; bounds on the LR scale.

| Field | Holds |
|---|---|
| `label as string` | "weak support" to "extremely strong support", or "no support" at LR = 1 |
| `supports as string` | `SUPPORTS_FIRST`, `SUPPORTS_ALTERNATIVE` or `SUPPORTS_NEITHER` |
| `statement as string` | the guideline's first phrasing, propositions exchanged below LR = 1 |
| `lower as float` | the band's lower LR bound |
| `lowerInclusive as bool` | whether `lower` belongs to the band |
| `upper as float` | the band's upper LR bound; 0.0 when `hasUpper` is false |
| `upperInclusive as bool` | whether `upper` belongs to the band |
| `hasUpper as bool` | false only for the open-ended top band |
## `mutation`

| Function | Returns | Does |
|---|---|---|
| `mutation.equalModel(rate)` | `MutationSpec` | Equal: `rate / (n-1)` to each other allele. |
| `mutation.proportionalModel(rate)` | `MutationSpec` | Proportional to the target's frequency; expected rate `rate`. |
| `mutation.stepwiseModel(rate, range)` | `MutationSpec` | Familias Stepwise: `range^steps` within a microvariant class, normalised to `rate`. |
| `mutation.extendedStepwiseModel(rate, range, rate2)` | `MutationSpec` | Stepwise plus `rate2` for changes between microvariant classes. |
| `mutation.stationaryModel(model)` | `MutationSpec` | Copy made stationary by the PM transformation. |
| `mutation.withLocusRate(model, locus, rate)` | `MutationSpec` | Copy with a per-locus rate. |
| `mutation.withLocusRates(model, rates)` | `MutationSpec` | Copy with several per-locus rates. |
| `mutation.locusAlleles(db, locus, typed, silent)` | `LocusAlleles` | A locus's alleles (table, typed, silent) and frequencies. |
| `mutation.mutationMatrix(model, alleles)` | `MutationMatrix` | The model's matrix over those alleles; rows sum to 1. |
| `mutation.transitionProbability(m, source, target)` | `float` | One matrix entry. |
| `mutation.overallRate(m, alleles)` | `float` | Expected rate `Σ π_i (1 - M_ii)`. |
| `mutation.isStationary(m, alleles)` | `bool` | Whether `π M = π` to 1e-9. |
| `mutation.parentageModel(paternal, maternal)` | `ParentageModel` | Sex-specific models, no silent allele. |
| `mutation.withSilentFrequency(model, locus, q)` | `ParentageModel` | Copy with a silent allele at one locus. |
| `mutation.paternityIndexAt(db, mother, child, father, model, theta)` | `float` | One-locus trio PI with both parents' mutation and silent alleles. |
| `mutation.paternityIndex(mother, child, allegedFather, db, model, theta)` | `float` | Trio CPI over the loci all three carry. |
| `mutation.duoPaternityIndexAt(db, child, father, model, theta)` | `float` | One-locus duo PI, paternal model and silent alleles. |
| `mutation.duoPaternityIndex(child, allegedFather, db, model, theta)` | `float` | Duo CPI over the loci both carry. |
| `mutation.kinshipLRAt(db, ga, gb, rel, alt, silent, theta)` | `float` | One-locus kinship LR with a silent allele. |
| `mutation.kinshipLR(a, b, rel, alt, db, silent, theta)` | `float` | Kinship LR with per-locus silent frequencies. |
| `mutation.aabb2003Counts()` | `list of MutationCount` | AABB 2003 mutation counts at 18 loci (STRBase). |
| `mutation.aabb2003Rates(which)` | `map of string to float` | Rates from those counts: `"total"`, `"paternal"` or `"maternal"`. |
| `mutation.kinshipModelWithRates(model, rates)` | `MutationModel` | A `kinship.MutationModel` with per-locus rates. |

| Constant | Type | Value | Means |
|---|---|---|---|
| `mutation.SILENT_ALLELE` | `string` | `"silent"` | Label of the silent allele in a matrix. |
| `mutation.MODEL_EQUAL` | `string` | `"equal"` | Equal model. |
| `mutation.MODEL_PROPORTIONAL` | `string` | `"proportional"` | Proportional model. |
| `mutation.MODEL_STEPWISE` | `string` | `"stepwise"` | Stepwise and extended stepwise. |
| `mutation.NO_MUTATION` | `MutationSpec` | `MutationSpec{model: "equal", rate: 0.0, range: 1.0, rate2: 0.0, stationary: false, perLocus: {}}` | The identity matrix. |
| `mutation.MENDELIAN` | `ParentageModel` | `ParentageModel{paternal: NO_MUTATION, maternal: NO_MUTATION, silent: {}}` | No mutation, no silent allele. |

### `mutation.MutationSpec`

A mutation model before it meets a locus's alleles.

| Field | Holds |
|---|---|
| `model as string` | `MODEL_EQUAL`, `MODEL_PROPORTIONAL` or `MODEL_STEPWISE` |
| `rate as float` | the per-meiosis mutation rate, in [0, 1) |
| `range as float` | Stepwise: the probability of `n+1` steps relative to `n`, above 0 |
| `rate2 as float` | Stepwise: the rate of non-integer-step (microvariant) changes |
| `stationary as bool` | whether to make the matrix stationary (PM) |
| `perLocus as map of string to float` | locus -> rate, overriding `rate` |

### `mutation.LocusAlleles`

A locus's alleles in matrix order and their frequencies.

| Field | Holds |
|---|---|
| `locus as string` | the locus name |
| `alleles as list of string` | the allele labels |
| `frequencies as list of float` | floored frequencies times `1 - q`; `q` for the silent allele |

### `mutation.MutationMatrix`

`rows[i][j]` = P(`alleles[i]` transmitted as `alleles[j]`).

| Field | Holds |
|---|---|
| `locus as string` | the locus name |
| `alleles as list of string` | the row and column labels |
| `rows as list of list of float` | the transition probabilities |

### `mutation.ParentageModel`

Sex-specific mutation models and per-locus silent frequencies.

| Field | Holds |
|---|---|
| `paternal as MutationSpec` | the model for a father's transmissions |
| `maternal as MutationSpec` | the model for a mother's transmissions |
| `silent as map of string to float` | locus -> silent allele frequency; absent means 0 |

### `mutation.MutationCount`

Published mutation counts at one locus.

| Field | Holds |
|---|---|
| `locus as string` | the locus, as `panels.normalizeLocus` names it |
| `maternalMutations as int` | mutations assigned to the mother |
| `maternalMeioses as int` | maternal meioses observed |
| `paternalMutations as int` | mutations assigned to the father |
| `paternalMeioses as int` | paternal meioses observed |
| `unassignedMutations as int` | mutations from either parent |
| `totalMutations as int` | all mutations |
| `totalMeioses as int` | all meioses |
