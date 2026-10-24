# `lineage`

```jennifer
import "@deckyard/forensicgenetics/lineage.j" as lineage;
```

mtDNA and Y-STR haplotypes. Independent of the rest of the deck. Signatures,
constants and struct fields are in the [Cheatsheet](cheatsheet.md#lineage).

## Conventions

- Counts: `N > 0`, `k` in `[0, N]`, confidence in `(0.5, 1)`. The bound is
  exactly 1.0 at `k == N`, so the LR there is 1.
- Positions are 1-based rCRS coordinates in `[1, 999999]`; insertion indices
  in `[1, 9999]`. Parsing ignores case and surrounding whitespace.
- Differences sort by position, then insertion index.
- `mtMatches` and `mtDistance` ignore sequenced ranges; `mtRangesAgree`
  compares them as sets.
- `align` upper-cases its input, ignores whitespace, accepts IUPAC codes, and
  slides indels to the 3' end of homopolymer runs. No reference is bundled.
  EMPOP's phylogenetic alignment rules are not applied.
- Y-STR multi-copy loci are written `"11,14"` or `"11/14"` and compared
  sorted. Only shared loci are compared.

Guide: [Lineage markers](../casework/lineage.md). Formulae:
[Lineage markers](../project/formulae.md#lineage-markers).

## Errors

`Error{kind: "lineage"}`, thrown for:

- a non-positive `N`, a `k` outside `[0, N]`, or a confidence outside
  `(0.5, 1)`
- a malformed difference designation, a position outside `[1, 999999]`, an
  insertion index of 0 or over four digits, or a deleted insertion
- an unknown `kind` in `formatDifference`
- in `align` and `mtHaplotypeFrom`: an empty sequence, a non-positive
  `refStart`, a non-IUPAC character, or a sequence longer than
  `MAX_ALIGN_LENGTH`
- an empty locus name, a locus with no allele, or a Y-STR allele that is not a
  repeat count
- two Y haplotypes sharing no locus, or with different copy counts at a shared
  locus
- a negative `steps`
