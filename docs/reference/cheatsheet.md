<!-- SPDX-License-Identifier: LGPL-3.0-only -->
<!-- SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev> -->

# forensicgenetics cheatsheet

Every export of `@deckyard/forensicgenetics`. The module pages
([forensicgenetics](forensicgenetics.md), [kinship](kinship.md),
[pedigree](pedigree.md), [lineage](lineage.md)) list each module's rules and
errors.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/kinship.j" as kinship;
import "@deckyard/forensicgenetics/lineage.j" as lineage;
import "@deckyard/forensicgenetics/pedigree.j" as pedigree;
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
