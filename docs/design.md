# Design

## Scope

Forensic and statistical genetics as pure computation over profiles and
allele-frequency tables the caller supplies:

- STR random match probabilities with the NRC II theta correction, and
  across an uncertain population;
- mixture CPI / CPE and the semi-continuous (LRmix-style) mixture LR;
- pairwise kinship and trio / duo paternity, with mutation models and silent
  alleles following Familias;
- whole-pedigree likelihoods by Elston-Stewart peeling;
- mtDNA and Y-STR lineage statistics, rCRS difference notation and a banded
  alignment that produces it, and EMPOP's EMP files of mtDNA haplotypes;
- X-STR kinship and daughter paternity;
- population statistics (forensic efficiency, HWE tests, F-statistics, LD
  tests), kit panels, direct and familial database search, and reporting aids
  (size-bias correction, a posterior bound on the RMP, the ENFSI verbal
  scale).

The deck does no I/O, so it needs no capability and runs on `jennifer-tiny`.
The one file format it knows, EMP, it takes and returns as text.

Out of scope:

- **Reading published tables and run files.** Reader decks build the deck's
  values: `@deckyard/strider` turns STRidER frequency files into a
  `FrequencyDb`. The core never imports a reader deck. EMP is the exception
  (`emp.j`): a line of it is a `lineage.MtHaplotype` in the notation
  `lineage` already parses, so a separate deck would hold little but a call
  into this one. Fragment analysis
  (capillary-electrophoresis traces to called alleles, over the `abif` deck)
  is planned as separate work.
- **Fully continuous mixture models** (peak heights, EuroForMix-style
  deconvolution); the mixture LR stops at the semi-continuous model.
- **Casework use.** The deck is for research, teaching and cross-checks; it
  is not a validated forensic tool. [Limitations](project/limitations.md)
  states this and every statistical boundary in detail.

## Decisions

- **Twelve modules, each importable alone.** The core
  (`forensicgenetics.j`) holds the shared types (`FrequencyDb`,
  `LocusFrequencies`, `Genotype`, `Profile`, `Mixture`), the frequency floor
  and NRC II genotype frequencies. `kinship`, `pedigree`, `mixture`,
  `lineage`, `popstats`, `panels`, `search`, `xstr`, `reporting` and
  `mutation` sit on it; `lineage` imports nothing from the deck, and `emp`
  only `lineage`. A downstream deck that needs only mtDNA notation imports
  `lineage.j` alone.
- **EMP files are read strictly.** EMPOP's own checker, EMPCheck, refuses a
  malformed file rather than guessing, and so does `readEmp`, naming the line.
  It runs the checks that need no reference: notation, positions within
  1-16569 and the ranges, insertions numbered without a gap, unique sample
  IDs. EMPOP's template marks the ranges with a placeholder (`#! SEQUENCE
  RANGE 1`), which is refused too. Differences are written in range order,
  as EMPOP lists them, so a file read and written again keeps its lines.
- **Pure `.j`, a deck and not core.** Probability arithmetic over maps fits
  the language. The two places where the interpreter's speed binds use
  built-ins: the Clopper-Pearson bound comes from `stats.proportionCi`,
  and the stationary mutation matrix is solved with `linalg.solve`
  (interpreted, the cubic loop took 2.3 s at 120 alleles).
- **Frequencies are stored as published.** `alleleFrequency` applies the
  `MIN_ALLELE_COUNT / (2N)` floor (5/2N) and a cap of 1 on read. `n == 0`
  means no floor, and an unlisted allele at such a locus throws. `withLocus`
  refuses a sum above `MAX_FREQUENCY_SUM` (1.05) as a misread table; partial
  tables summing below 1 are allowed.
- **Alleles are strings** (`"9.3"`, `"X"`); `"15"` and `"15.0"` are
  different alleles. `genotype` puts the pair in a canonical order.
- **Never drop a locus silently.** Match probabilities, kinship, paternity and
  pedigree calculations throw when a typed locus is missing from the table,
  because dropping it can hide an exclusion. Callers narrow deliberately with
  `commonLoci` and `restrict`.
- **An uncertain population is not a pooled one.** `combineFrequencies`
  pools by counts, but the docs steer to `randomMatchProbabilityAcross`:
  pooling error is not conservative in either direction.
- **Theta where the formulas carry it.** Genotype frequencies are
  Balding-Nichols (NRC II 4.2). Kinship and paternity have `...Theta`
  variants using the sampling formula; at theta 0 they delegate to the
  originals, so results are identical. Pedigree peeling and X-STRs use theta
  0, since peeling does not carry the joint founder sample the formula needs.
- **Stable functions are frozen.** The original kinship functions and the
  `MutationModel` struct are a downstream contract: a struct literal missing
  a field is a runtime error, so new models are a separate `MutationSpec` and
  `ParentageModel` rather than new fields.
- **Overflow is an error, with a way out.** Products across loci are checked
  against the float maximum before multiplying and throw naming the `log10`
  form (`log10KinshipLR`, `log10PedigreeLR`, `log10MixtureLR`, ...), since
  the interpreter raises a runtime error on overflow.
- **Placeholders are labelled.** `DEFAULT_MUTATION` is explicitly not a
  published rate. The AABB 2003 mutation counts were transcribed from an
  archived copy of STRBase and are checked row by row against their totals.
- **Sources over notation.** The Familias models follow the Familias R
  package source (including R's half-to-even rounding of microvariant
  classes), not a paper's notation. X-STR relationship constants are derived
  from X transmission and checked against published coancestry tables.
  Formulas whose source could not be read (X-STR power of exclusion) are not
  implemented. [Formulae and sources](project/formulae.md) lists every
  formula with its reference.
- **Own random streams.** The permutation tests and the posterior bound use
  their own seeded Park-Miller generator, so they neither read nor disturb
  `math.rand`, and a seed reproduces a result.
- **Exact where float would mislead.** Y-STR step distances, mutation step
  differences and the posterior's nearest-rank quantile use `decimal` where a
  float would round to the wrong integer.
- **Performance shapes the code.** Jennifer copies values: iterating a list
  of databases, folding `withLocus` per entry, passing nested lists, or
  calling into another module with that module's structs all copy whole
  values. The modules index instead of iterating, write fields directly,
  wrap message vectors in structs, avoid callbacks in functions that take the
  database (a function value in the call tree makes compound arguments copy),
  and memoise per locus and genotype in familial search. These patterns are
  commented where they occur and should be kept.

## Errors

Every module throws `Error{kind: <module>}` through its own `fail`: kinds
`"forensicgenetics"`, `"kinship"`, `"pedigree"`, `"mixture"`, `"lineage"`, `"emp"`,
`"popstats"`, `"panels"`, `"search"`, `"xstr"`, `"reporting"` and
`"mutation"`. A module calling another rethrows that module's error unchanged
unless documented otherwise (the X-STR daughter indices rethrow
`"kinship"` as `"xstr"`; `emp` rethrows `"lineage"` as `"emp"` with the line
number in front). `@deckyard/strider` relies on the core kind
`"forensicgenetics"` for `withLocus`.

The caps below throw the generic `Error{kind: "limit"}` (since 0.5.0; they
were the module's kind before). Ranges on a caller's own parameters
(`MAX_PERMUTATIONS`, `MAX_DRAWS`, theta, rates) and malformed tables
(`MAX_FREQUENCY_SUM`) are validation and keep the module's kind.

## Security and limits

Input is profiles and tables, often from other tools; the caps keep a
pathological input from running for hours or exhausting memory:

- Mixtures: `MAX_UNKNOWNS` (3), `MAX_LOCUS_ALLELES` (10) and
  `MAX_COMBINATIONS` (250,000 for a whole call, counted before any locus is
  evaluated; about 11 s).
- Pedigrees: `MAX_PEDIGREE_ALLELES` (20) observed alleles per locus; simple
  pedigrees only (loops are refused). Peeling is iterative, so depth needs no
  cap.
- Alignment: `MAX_ALIGN_LENGTH` (5000) and `MAX_ALIGN_CELLS` (1.25 million
  per pass), so oversize or too divergent input fails instead of appearing to
  hang. The band doubles until a bound proves the banded result equals the
  full matrix.
- X-STR: `MAX_GROUP_HETEROZYGOUS_LOCI` (12) per linkage group.
- EMP files: `MAX_ENTRIES` (200,000) haplotypes per file, several times
  EMPOP's whole database; reading that many takes about 100 s.
- Mutation: `MAX_STATIONARY_ALLELES` (200) for the stationary model.
- Monte Carlo: `MAX_PERMUTATIONS` and `MAX_DRAWS` (1,000,000 each), checked
  as a range on the caller's parameter (module kind).
- Population statistics refuse a locus whose listed frequencies sum below
  `MIN_FREQUENCY_SUM` (0.95).

## Testing and validation

- Twelve overlays, about 520 tests, no fixtures and no network. The EMP
  tests use files written by hand; the reader was also run over two files
  EMPOP publishes (273 and 202 haplotypes), which read and write back to the
  same lines, but those files are not in the repository.
- Expectations are derived by hand or from published sources, never captured
  from the code: the Genepop documentation's F-statistics example, exact
  enumeration of LD re-pairings, closed-form kinship values that peeling must
  reproduce to 1e-9, a full-width aligner the banded one must agree with, and
  published coancestry tables for X-STRs.
- Cross-module checks: functions that delegate (theta 0, silent frequency 0,
  familial search against `kinshipLR`) are asserted equal, not just close.
- [Testing and contributing](project/testing.md) describes how expectations
  are chosen.
- Not validated against Familias, LRmix or EuroForMix output on casework
  data.

## Known limitations and open questions

[Limitations](project/limitations.md) is the full list. The main ones: no
theta in pedigrees or X-STRs, simple pedigrees only, mixtures stop at the
semi-continuous LR, EMPOP's phylogenetic alignment rules are not applied,
linked X-STRs are scored without recombination, and the default mutation
rates are placeholders.

## Dependencies on the interpreter

Floor `>=0.27.0` for `jennifer` and `jennifer-tiny`, set by the UTF-8 codec
of `encoding.encode`, `strings.replicate` and `testing.assertThrowsArgs`;
`decimal` needs 0.26.0. It uses
`stats` (the Clopper-Pearson interval) and `linalg` from the built-in libraries. Two of
`jennifer`'s planned milestones were found here and are worked around in the
code: M30.3 (borrow nested collections; peeling copied 2.3 GiB on a
300-generation chain) and M30.8 (borrow through named callbacks; `search`
copied its whole database per call).

Every call that takes an `emp.EmpFile` copies it across the module
boundary (the M30 copy-elision work removes this): at 5000 haplotypes about
0.06 s on `jennifer`, 1.7 s on `jennifer-tiny`.

Two decks depend on this one by tag: `@deckyard/strider` (`FrequencyDb`
construction, `MAX_FREQUENCY_SUM`, the core error kind) and `@deckyard/mitoleaf`
(`lineage`'s difference notation and distance). Their behaviour changes only
after a release here.
