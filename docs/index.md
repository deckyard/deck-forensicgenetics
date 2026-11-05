# @deckyard/forensicgenetics

Forensic and statistical genetics for [Jennifer](https://jennifer-lang.dev/).
The deck works on profiles (an unordered allele pair per autosomal STR locus,
plus mtDNA and Y-STR haplotypes) and on reference frequency or count data. It
does not handle sequences, apart from aligning an mtDNA region to name its
differences.

This deck is part of [Deckyard](https://github.com/deckyard), community decks for [Jennifer](https://jennifer-lang.dev/).

## What it computes

| Question | Statistic | Guide |
| --- | --- | --- |
| Could the stain have come from someone else? | RMP with NRC II theta, single-source LR | [Match probability](casework/match-probability.md) |
| How much of the population is not excluded from a mixture? | PI / PE, CPI / CPE | [Mixtures](casework/mixtures.md) |
| Did the suspect contribute to the mixture? | semi-continuous LR with drop-out and drop-in | [Mixtures](casework/mixtures.md#likelihood-ratio-with-drop-out-and-drop-in) |
| Are these two people related, and how? | kinship LR from IBD coefficients | [Kinship](casework/kinship.md) |
| Is this man the child's father? | trio paternity index with a stepwise mutation model | [Paternity](casework/paternity.md) |
| What if a parent mutated, or an allele is silent? | Familias mutation models (Equal, Proportional, Stepwise, Extended stepwise, PM-stationary), sex-specific, published rates; silent alleles | [Mutation models and silent alleles](casework/mutation-models.md) |
| What does the whole family say? | pedigree likelihood by Elston-Stewart peeling | [General pedigrees](casework/pedigrees.md) |
| How common is this haplotype? | Clopper-Pearson upper bound on a database count | [Lineage markers](casework/lineage.md) |
| How informative is this table, and does its sample fit HWE? | He, PD, PIC, TPI, PE; HWE exact and permutation tests | [Population statistics](concepts/population-statistics.md) |
| How different are the sub-populations, and are the loci independent? | Weir & Cockerham F_ST / F_IT / F_IS; genotypic LD permutation test | [Population structure and linkage](concepts/population-structure.md) |
| Who in a database matches, or is related to, this profile? | direct comparison with IBS counts; familial ranking by kinship LR | [Database searching](casework/database-search.md) |
| Are these women paternal half sisters? Is he this girl's father? | X-STR kinship LR by X-IBD and sex, daughter trio / duo index, linkage-group haplotypes | [X-chromosomal STRs](casework/x-chromosome.md) |
| How much does the database's size matter, and how is the LR worded? | size-bias correction, Dirichlet posterior bound on the RMP, ENFSI verbal scale | [Reporting a match](casework/reporting.md) |

## Modules

| Module | Import | Covers |
| --- | --- | --- |
| `forensicgenetics.j` | `@deckyard/forensicgenetics/` | profiles, `FrequencyDb`, match probability, mixtures |
| `kinship.j` | `@deckyard/forensicgenetics/kinship.j` | relationship LRs, paternity |
| `pedigree.j` | `@deckyard/forensicgenetics/pedigree.j` | whole-pedigree likelihoods |
| `lineage.j` | `@deckyard/forensicgenetics/lineage.j` | mtDNA and Y-STR |
| `mixture.j` | `@deckyard/forensicgenetics/mixture.j` | mixture LR with drop-out and drop-in |
| `popstats.j` | `@deckyard/forensicgenetics/popstats.j` | forensic-efficiency table, HWE tests, F-statistics, LD tests |
| `panels.j` | `@deckyard/forensicgenetics/panels.j` | kit locus panels, locus-name normalisation |
| `search.j` | `@deckyard/forensicgenetics/search.j` | direct and familial database search |
| `xstr.j` | `@deckyard/forensicgenetics/xstr.j` | X-STR match probability, kinship, daughter paternity |
| `reporting.j` | `@deckyard/forensicgenetics/reporting.j` | sampling uncertainty of an RMP, verbal scale for an LR |
| `mutation.j` | `@deckyard/forensicgenetics/mutation.j` | Familias mutation models, silent alleles, AABB rates |

Every module except `lineage` imports the core (`search`, `xstr` and
`mutation` also `kinship`); `lineage` imports nothing from the
deck. The core reads no published frequency format: a reader deck such as
[`@deckyard/strider`](https://github.com/deckyard/deck-strider-frequency-converter)
(STRidER XML) converts one into a `FrequencyDb`.

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;

def db as forensicgenetics.FrequencyDb init forensicgenetics.withLocus(
    forensicgenetics.db("example"), "D3S1358", 1034,
    {"15": 0.2456, "16": 0.2314, "17": 0.2007, "18": 0.3223});
def suspect as forensicgenetics.Profile init forensicgenetics.profile("S1", [
    forensicgenetics.genotype("D3S1358", "15", "16")
]);
def lr as float init forensicgenetics.matchLikelihoodRatio($db, $suspect,
    forensicgenetics.THETA_GENERAL);
```

Every module throws `Error` with its own `kind`: `"forensicgenetics"`,
`"kinship"`, `"pedigree"`, `"lineage"`, `"mixture"`, `"popstats"`, `"panels"`,
`"search"`, `"xstr"`, `"reporting"` or `"mutation"`.

Start with [Installing the deck](getting-started/installation.md) and the
[Quick start](getting-started/quick-start.md). Every formula and its source is
in [Formulae and sources](project/formulae.md); every export is in the
[Cheatsheet](reference/cheatsheet.md). Read [Limitations](project/limitations.md)
before relying on a number.

## Not a validated forensic tool

This deck implements published estimators and is tested against them, but it
has not been through the validation any jurisdiction requires of software used
in casework. Treat it as a library for research, teaching and cross-checking,
and validate independently before it informs a real conclusion.

The frequency table the examples build (`examples/lib/exampledb.j`) is
synthetic: made-up numbers that exist only to exercise the examples. Real
reference data comes from [STRidER](https://strider.online/) and elsewhere; see
[Reference frequency data](concepts/frequency-data.md).

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
