# `emp`

```jennifer
import "@deckyard/forensicgenetics/emp.j" as emp;
```

EMPOP's EMP files: mtDNA haplotypes read from and written to text. Imports
`lineage` only. Signatures, constants and struct fields are in the
[Cheatsheet](cheatsheet.md#emp).

## Conventions

- Text in, text out: the module does no I/O. Read a file with
  `fs.readString` and pass its content.
- A haplotype is a `lineage.MtHaplotype`: the sample ID is its `id`, the
  ranges of the `#!` line above it its `ranges`.
- Lines end in LF or CRLF; a leading byte-order mark, blank lines and empty
  trailing cells are accepted. `writeEmp` writes LF.
- `comments` holds the comment lines above the first `#!` line. Comments
  further down, a haplotype commented out among them, are skipped and not
  written back.
- Sample IDs hold letters, digits, `-`, `_` and `/` only, and each appears
  once. A haplogroup is one word, `?` if unknown. A count is 0 to
  `MAX_COUNT`; 0 keeps the entry but excludes it from `empSampleSize` and
  `empCount`.
- Ranges are `start-end` within 1-16569; a start above the end runs across
  the origin (`16024-576`). Several ranges on one `#!` line are separated by
  spaces.
- Differences use `lineage`'s notation, one per cell. Beyond it, as EMPOP's
  EMPCheck: positions within 1-16569 and within the ranges, and insertions
  numbered `.1`, `.2`, ... without a gap. A difference that names the rCRS
  base itself is not caught: no reference is bundled.
- `writeEmp` writes differences in range order (16xxx before 73 for
  `16024-576`), deletions as `DEL`, and an empty last cell for a haplotype
  without differences, as EMPOP's template does.
- `empCount` counts the entries with the same differences and the same ranges
  (as sets), as `lineage.mtMatches` and `lineage.mtRangesAgree` decide.

Guide: [Lineage markers](../casework/lineage.md#empops-emp-files).

## Errors

`Error{kind: "limit"}`, thrown by `readEmp` past `MAX_ENTRIES` haplotypes.

`Error{kind: "emp"}`, thrown for:

- in `readEmp`, naming the line: a haplotype before the first `#!` line, a
  `#!` line without a range, a malformed or out-of-genome range, a sample ID
  with other characters or seen before, an empty or multi-word haplogroup, a
  count that is not a whole number in 0-9999, two variants in one cell, a
  difference `lineage.parseDifference` refuses, a position named twice, a
  position outside 1-16569 or outside the ranges, an insertion `.n` without
  `.n-1`
- in `empEntry`, naming the sample, and in `writeEmp`, naming the entry
  (1-based): the same, and in `writeEmp` a repeated sample ID or a comment
  holding a line break
- an unknown difference `kind` in `empDesignation`

A malformed difference throws `"lineage"` in `lineage`; here it is `"emp"`,
with the line in front of `lineage`'s message.
