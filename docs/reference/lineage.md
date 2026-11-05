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
  in `[1, 9999]`. Position 0 is accepted for an insertion only (`0.1G`, a base
  before position 1, as `align` emits it). Parsing ignores case and
  surrounding whitespace. A trailing `D` is a deletion, never the IUPAC code
  D.
- A haplotype names each position and insertion index at most once.
- Differences sort by position, then insertion index.
- `mtMatches` and `mtDistance` ignore sequenced ranges; `mtRangesAgree`
  compares them as sets (duplicates ignored).
- `align` upper-cases its input, ignores whitespace, accepts IUPAC codes, and
  returns the full-matrix optimum (the band widens until it is proven). Each
  indel run slides as a block to the 3' end of a homopolymer or tandem-repeat
  run. No reference is bundled. EMPOP's phylogenetic alignment rules are not
  applied.
- Y-STR multi-copy loci are written `"11,14"` or `"11/14"` and compared
  sorted. Only shared loci are compared. Step distances are exact decimal
  differences; alleles with different fractional parts (`17.2`, `18`) are not
  stepwise neighbours. `mtDistance` merges the two sorted difference lists; a
  struct literal that is unsorted or repeats a key is compared as designation
  sets, with the same result. `yMatchesWithin` is false where `yStepDistance` would
  throw for copy counts or fractional parts.

Guide: [Lineage markers](../casework/lineage.md). Formulae:
[Lineage markers](../project/formulae.md#lineage-markers).

## Errors

`Error{kind: "lineage"}`, thrown for:

- a non-positive `N`, a `k` outside `[0, N]`, or a confidence outside
  `(0.5, 1)`
- a malformed difference designation, a position outside `[1, 999999]`
  (other than an insertion at 0), an insertion index of 0 or over four
  digits, or a deleted insertion
- two designations of one haplotype at the same position and insertion index
- an unknown `kind` in `formatDifference`
- in `align` and `mtHaplotypeFrom`: an empty sequence, a non-positive
  `refStart`, a `refStart` that puts the last reference base past 999999, a
  non-IUPAC character, a sequence longer than
  `MAX_ALIGN_LENGTH`, or a band pass over `MAX_ALIGN_CELLS` (sequences too
  divergent)
- an empty locus name, a locus with no allele, or a Y-STR allele that is not a
  repeat count (1-6 digits, optionally `.` and 1-6 digits: `11.`, `.5` and
  `1234567` are not)
- two Y haplotypes sharing no locus; in `yStepDistance`, different copy counts
  at a shared locus or paired alleles with different fractional parts
- a negative `steps`
