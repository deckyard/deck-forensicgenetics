# `mixture`

```jennifer
import "@deckyard/forensicgenetics/mixture.j" as mixture;
```

The semi-continuous (LRmix) mixture likelihood ratio, with drop-out and
drop-in. It takes the core's `Mixture`, `Profile` and `FrequencyDb`.
Signatures, constants and the `Hypothesis` fields are in the
[Cheatsheet](cheatsheet.md#mixture).

## Conventions

- Evidence is allele presence only (the core `Mixture`); peak heights are not
  used.
- Drop-out is per allele copy: `d` for a heterozygote's allele, `d²` for a
  homozygote's, and an allele carried by several contributors is lost only if
  every copy drops out. Each known contributor has its own `d`; all unknowns
  share one.
- Drop-in: `c · p_a` per observed allele no contributor carries; `1 - c` when
  there is none. `c` is in `[0, 1)`, every `d` in `[0, 1]`.
- θ is applied through the Balding-Nichols sampling formula, conditioned on
  `conditioningProfiles(hp, hd)`: every known contributor of either
  hypothesis. `θ = 0` is Hardy-Weinberg.
- Every locus of the mixture is used. Each must be in the table and typed in
  every known profile, or the call throws; narrow with `restrictMixture`.
  Profile loci absent from the mixture are ignored.
- An empty allele list at a locus (build the `Mixture` struct directly;
  `forensicgenetics.mixture` rejects it) is valid evidence of complete
  drop-out.
- Cost per locus is about `(k(k+1)/2)^u / u!` terms for `k` allele classes and
  `u` unknowns. `MAX_UNKNOWNS` (3) and `MAX_LOCUS_ALLELES` (10) bound it;
  timings are in
  [Limitations](../project/limitations.md#the-mixture-lr-enumerates-genotypes).

Guide: [Mixtures](../casework/mixtures.md#likelihood-ratio-with-drop-out-and-drop-in).
Formulae: [Formulae and sources](../project/formulae.md#mixtures-semi-continuous-likelihood-ratio).

## Errors

`Error{kind: "mixture"}`, thrown for:

- a hypothesis with a drop-out list of the wrong length, a probability outside
  `[0, 1]`, unknowns outside `0..MAX_UNKNOWNS`, no contributor, or one profile
  id twice
- one profile id with different genotypes in the two hypotheses
- `dropIn` outside `[0, 1)`, `theta` outside `[0, 1)`
- a mixture with no loci; a locus the mixture or the table lacks
- a known or typed profile not typed at a mixture locus; a known contributor
  missing from `typed` (`evidenceProbability`)
- more than `MAX_LOCUS_ALLELES` alleles at a locus, or floored frequencies of
  those alleles summing past 1
- `P(E | Hd) = 0` at a locus (`locusLR` and the functions built on it)
- a zero locus LR in `log10MixtureLR`

An allele the table does not list at a locus with `n = 0` surfaces as the
core's `"forensicgenetics"` error.
