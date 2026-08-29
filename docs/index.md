# @jennifer/forensicgenetics

Forensic and statistical genetics for [Jennifer](https://jennifer-lang.dev/).

This deck works on **profiles** - an unordered allele pair per autosomal STR
locus, plus uniparental haplotypes for mtDNA and Y-STR - and on reference
**frequency or count** data. It does not work on sequences; that is the
sequence deck's job, and the two share no code.

## What it answers

Four questions, each with its own module:

**Could this stain have come from someone else?**
A random match probability from Hardy-Weinberg genotype frequencies with the
NRC II theta correction, multiplied across a panel - and its reciprocal, the
single-source likelihood ratio. For a mixed stain, the combined probability of
inclusion and exclusion instead. See [Match probability](casework/match-probability.md)
and [Mixtures](casework/mixtures.md).

**Are these two people related, and how?**
A likelihood ratio built from the pair's IBD coefficients, testing one
relationship against another. See [Kinship](casework/kinship.md).

**Is this man the child's father?**
A trio calculation over mother, child and alleged father, with a stepwise
mutation model so that one inconsistent locus is not a false exclusion. See
[Paternity](casework/paternity.md).

**What about everyone else in the family?**
The likelihood of a whole pedigree by Elston-Stewart peeling, which sums out
every untyped member - the question that arises the moment the person you
actually want to test is unavailable. See [General pedigrees](casework/pedigrees.md).

Plus [lineage markers](casework/lineage.md): mtDNA and Y-STR haplotypes, where
the estimate is a database count rather than a product over loci.

## The shape of it

```jennifer
import "@jennifer/forensicgenetics/" as forensics;

def db as forensics.FrequencyDb init forensics.withLocus(
    forensics.db("example"), "D3S1358", 1034,
    {"15": 0.2456, "16": 0.2314, "17": 0.2007, "18": 0.3223});

def suspect as forensics.Profile init forensics.profile("S1", [
    forensics.genotype("D3S1358", "15", "16"),
    forensics.genotype("TH01", "6", "9.3")
]);

def lr as float init forensics.matchLikelihoodRatio($db, $suspect,
    forensics.THETA_GENERAL);
```

Five modules, three of which you can ignore until you need them:

| Module | Import | Topic |
| --- | --- | --- |
| `forensicgenetics.j` | `@jennifer/forensicgenetics/` | Profiles, the generic frequency database, match probability, mixtures |
| `kinship.j` | `@jennifer/forensicgenetics/kinship.j` | Relationship LRs, paternity |
| `pedigree.j` | `@jennifer/forensicgenetics/pedigree.j` | Whole-pedigree likelihoods |
| `lineage.j` | `@jennifer/forensicgenetics/lineage.j` | mtDNA and Y-STR |

`kinship` and `pedigree` build on the core; `lineage` is independent and can be
used on its own. The core imports nothing - it defines the generic
`FrequencyDb` and reads no published format, so a reader deck for any source
converts into the same type. [`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter)
does that for STRidER-style XML.

Every module raises `Error` with its own `kind` - `"forensicgenetics"`,
`"kinship"`, `"pedigree"`, `"lineage"` - so a `catch` can tell an unusable
frequency table from an impossible pedigree.

## Where to go next

New here: [Installing the deck](getting-started/installation.md), then the
[Quick start](getting-started/quick-start.md).

Looking for a formula: [Formulae and sources](project/formulae.md) collects
every estimator in one place with its citation.

About to rely on a number: read [Limitations](project/limitations.md) first.
It is short, and some of it will change what you do.

## Not a validated forensic tool

This deck implements published estimators and is tested against them, but it
has not been through the validation any jurisdiction requires of software used
in casework. Treat it as a library for research, teaching and cross-checking,
and validate independently before it informs a real conclusion.

The bundled frequency table is synthetic - made-up numbers that exist only to
exercise the reader and the examples. Real reference data comes from
[STRidER](https://strider.online/) and elsewhere; see
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
