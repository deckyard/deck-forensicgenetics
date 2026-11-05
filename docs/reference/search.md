# `search`

```jennifer
import "@deckyard/forensicgenetics/search.j" as search;
```

Direct and familial database searching: one query `Profile` against a
`list of Profile`. Signatures, verdict constants and struct fields are in the
[Cheatsheet](cheatsheet.md#search).

## Conventions

- Loci are compared by name as written. Normalise both sides with
  `panels.normalizeProfile` when they come from different sources.
- Only loci both profiles carry are compared; `lociCompared` counts them.
- `allelesShared` is IBS: the size of the multiset intersection of the two
  allele pairs (`10/10` vs `10/11` is 1).
- With `allowDropout`, a query homozygote `a/a` against a candidate
  heterozygote carrying `a` is a drop-out locus, not a mismatch. The test is
  directional: the query is the evidence.
- Verdicts: `"insufficient"` if `lociCompared < minLoci`; else `"exact"` if
  every compared locus is identical; else `"partial"` if mismatching loci
  `<= maxMismatches`; else `"mismatch"`.
- `search` ranks by verdict, then fewer mismatches, more identical loci, more
  shared alleles, then database order.
- `familialSearch` scores each candidate sharing at least `minLoci` loci with
  the value of `kinship.kinshipLR` (theta 0), ranks by LR descending, and puts
  unscored candidates last with `lr` 0.0. An LR of 0.0 on a scored candidate
  is an exclusion (for example parent-child with no shared allele).

Guide: [Database searching](../casework/database-search.md).

## Errors

`Error{kind: "search"}`, thrown for:

- `minLoci` below 1 or a negative `maxMismatches`, including an options
  literal that bypassed `searchOptions`
- an empty query profile
- `allelesShared` on genotypes at different loci

`familialSearch` also raises `kinship`'s `"kinship"` error when the query and
a candidate share a locus the table does not cover.
