# `xstr`

```jennifer
import "@deckyard/forensicgenetics/xstr.j" as xstr;
```

X-chromosomal STRs: male and female profiles, match probabilities, pairwise
kinship LRs from X-IBD coefficients, and paternity indices for a daughter.
Signatures, relationship constants and struct fields are in the
[Cheatsheet](cheatsheet.md#xstr).

## Conventions

- A profile carries its sex. A male's genotype at a locus has `a == b`: his
  one allele, not a homozygote. `maleProfile` builds it from
  `locus -> allele`; `femaleProfile` from core `Genotype`s.
- An `XRelationship` is ordered: person A has `sexA`, person B `sexB`, and
  the profiles passed must match. `FATHER_DAUGHTER` takes the man first.
  `k2` is 0 for any pair with a male.
- `XFrequencies` is the model. `independentModel(db)` scores every locus
  from a core `FrequencyDb` (floored at `5/(2N)` as in the core).
  `withLinkageGroup` makes a cluster's loci one marker, scored from the
  group's haplotype table: a male's haplotype frequency, a female's genotype
  summed over the haplotype pairs (phases) her genotypes allow.
- A haplotype key is the alleles in the group's locus order joined by `-`
  (`haplotypeKey(["17", "24", "12"])` is `"17-24-12"`). Unlisted haplotypes
  score at the group's `minFrequency`; with 0 they throw.
- A locus or group is scored when every person is typed at all of it; one
  person typed at it alone is skipped. A person typed at part of a group
  throws, as does a single locus every person carries that the database
  lacks.
- `kinshipLR` multiplies over loci and groups; `pairUnits` lists them, a
  group by its name.
- `daughterTrioIndex` and `daughterDuoIndex` take a `kinship.MutationModel`
  for the alleged father: within a group, mutation is per locus and
  independent. The mother does not mutate.
- Theta is 0 throughout.

Guide: [X-chromosomal STRs](../casework/x-chromosome.md).

## Errors

`Error{kind: "xstr"}`, thrown for:

- an unknown sex; a male profile with two alleles at a locus; an empty allele
  or a repeated locus
- coefficients outside [0, 1], not summing to 1, or `k2 > 0` with a male
- a relationship whose sexes differ from the profiles'
- a malformed linkage group (no loci, repeated loci, a floor outside [0, 1),
  a key with the wrong number of alleles, a frequency outside [0, 1], a sum
  past `forensicgenetics.MAX_FREQUENCY_SUM`), a group name or locus already in
  the model, an unlisted haplotype with no floor
- no shared locus, a shared locus nothing covers, a partly typed group
- an alternative hypothesis giving the observations probability zero, the
  `log10` of an exclusion, a maternal exclusion in a trio

An unlisted allele at a locus with `n == 0` raises the core's
`"forensicgenetics"` error; a non-numeric allele under a mutation model,
`kinship`'s `"kinship"` error.
