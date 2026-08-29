# @jennifer/forensicgenetics

[![test](https://github.com/jennifer-language/deck-forensicgenetics/actions/workflows/test.yml/badge.svg)](https://github.com/jennifer-language/deck-forensicgenetics/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://jennifer-language.github.io/deck-forensicgenetics/)
[![release](https://img.shields.io/github/v/release/jennifer-language/deck-forensicgenetics)](https://github.com/jennifer-language/deck-forensicgenetics/releases/latest)

Forensic and statistical genetics for [Jennifer](https://jennifer-lang.dev/) -
the statistical-genetics sibling of the sequence deck. It works on **profiles**
(an unordered allele pair per autosomal STR locus, plus uniparental haplotypes
for mtDNA and Y-STR) and on reference **frequency / count** data, not on
sequences.

Requires Jennifer 0.25.0 or newer. Pure Jennifer over the system libraries: no
deck dependencies, and nothing that stops it running on `jennifer-tiny`.

## Documentation

**<https://jennifer-language.github.io/deck-forensicgenetics/>**

The manual is the whole story - installing, the statistical model behind each
estimator, worked examples, the complete API reference, and the limitations
that matter before you quote a number.

| | |
| --- | --- |
| Start here | [Quick start](https://jennifer-language.github.io/deck-forensicgenetics/getting-started/quick-start.html), [Installing the deck](https://jennifer-language.github.io/deck-forensicgenetics/getting-started/installation.html) |
| Casework | [match probability](https://jennifer-language.github.io/deck-forensicgenetics/casework/match-probability.html), [mixtures](https://jennifer-language.github.io/deck-forensicgenetics/casework/mixtures.html), [kinship](https://jennifer-language.github.io/deck-forensicgenetics/casework/kinship.html), [paternity](https://jennifer-language.github.io/deck-forensicgenetics/casework/paternity.html), [pedigrees](https://jennifer-language.github.io/deck-forensicgenetics/casework/pedigrees.html), [lineage markers](https://jennifer-language.github.io/deck-forensicgenetics/casework/lineage.html) |
| API reference | [forensicgenetics](https://jennifer-language.github.io/deck-forensicgenetics/reference/forensicgenetics.html), [kinship](https://jennifer-language.github.io/deck-forensicgenetics/reference/kinship.html), [pedigree](https://jennifer-language.github.io/deck-forensicgenetics/reference/pedigree.html), [lineage](https://jennifer-language.github.io/deck-forensicgenetics/reference/lineage.html) |
| The project | [Limitations](https://jennifer-language.github.io/deck-forensicgenetics/project/limitations.html), [Formulae and sources](https://jennifer-language.github.io/deck-forensicgenetics/project/formulae.html), [Testing and contributing](https://jennifer-language.github.io/deck-forensicgenetics/project/testing.html) |

## At a glance

```jennifer
import "@jennifer/forensicgenetics/" as forensics;

def db as forensics.FrequencyDb init forensics.withLocus(
    forensics.db("example"), "D3S1358", 1034,
    {"15": 0.2456, "16": 0.2314, "17": 0.2007, "18": 0.3223});
def suspect as forensics.Profile init forensics.profile("S1", [
    forensics.genotype("D3S1358", "15", "16"),
    forensics.genotype("TH01", "6", "9.3")
]);

def rmp as float init forensics.randomMatchProbability($db, $suspect,
    forensics.THETA_GENERAL);
def lr as float init forensics.matchLikelihoodRatio($db, $suspect,
    forensics.THETA_GENERAL);
```

| Module | Import | What it does |
| --- | --- | --- |
| `forensicgenetics.j` | `@jennifer/forensicgenetics/` | Vendor-independent core types, Hardy-Weinberg with the NRC II theta correction, RMP, single-source LR, mixture CPI / CPE |
| `kinship.j` | `@jennifer/forensicgenetics/kinship.j` | IBD coefficients, relationship likelihood ratios, trio paternity with a stepwise mutation model |
| `pedigree.j` | `@jennifer/forensicgenetics/pedigree.j` | General-pedigree likelihoods by Elston-Stewart peeling |
| `lineage.j` | `@jennifer/forensicgenetics/lineage.j` | mtDNA and Y-STR haplotypes, Clopper-Pearson counting estimates, rCRS difference nomenclature |

The core reads no published format at all - it defines the generic
`FrequencyDb`, and a reader deck such as
[`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter)
converts a published table into it. See
[Reference frequency data](https://jennifer-language.github.io/deck-forensicgenetics/concepts/frequency-data.html).

## Not a validated forensic tool

This deck implements published estimators and is tested against them, but it
has not been through the validation any jurisdiction requires of software used
in casework. Treat it as a library for research, teaching and cross-checking,
and validate independently before it informs a real conclusion.

The bundled frequency table is synthetic - made-up numbers that exist only to
exercise the examples. The full list of caveats is in
[Limitations](https://jennifer-language.github.io/deck-forensicgenetics/project/limitations.html);
it is short, and some of it will change what you do.

**Do not use this deck for a legal proceeding** - a prosecution, a defence, a
paternity ruling, an immigration decision, or any other legal determination.
Those require validated, accredited software operated under an accredited
process, and this is not that.

**No warranty, no liability.** This software is provided as is, without
warranty of any kind, express or implied. You use it entirely at your own risk.
The author accepts no responsibility and no liability for any claim, damage,
loss, or other consequence arising from this software or from any number it
produces. Checking that a result is correct, appropriate and lawful for your
purpose is yours alone. This restates, and does not narrow or replace, the
disclaimers in the [licence](LICENSE); where they differ, the licence governs.

Full text: [Not a validated forensic tool](https://jennifer-language.github.io/deck-forensicgenetics/#not-a-validated-forensic-tool).

## Licence

LGPL-3.0-only - see [LICENSE](LICENSE).
