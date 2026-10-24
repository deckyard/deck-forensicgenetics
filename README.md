# @deckyard/forensicgenetics

[![test](https://github.com/deckyard/deck-forensicgenetics/actions/workflows/test.yml/badge.svg)](https://github.com/deckyard/deck-forensicgenetics/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-github%20pages-blue)](https://deckyard.github.io/deck-forensicgenetics/)
[![release](https://img.shields.io/github/v/release/deckyard/deck-forensicgenetics)](https://github.com/deckyard/deck-forensicgenetics/releases/latest)

Forensic and statistical genetics for [Jennifer](https://jennifer-lang.dev/):
STR random match probabilities with the NRC II theta correction, mixture
CPI / CPE, kinship and trio paternity likelihood ratios, whole-pedigree
likelihoods by Elston-Stewart peeling, and mtDNA / Y-STR haplotype statistics.
It works on profiles and reference allele frequencies, not on sequences.
Requires Jennifer 0.25.0 or newer; no deck dependencies, and it runs on
`jennifer-tiny`.

## Install

Copy the deck's `src/*.j`, without the `*_test.j` overlays, into
`vendor/deckyard/forensicgenetics/`. A release archive is
`deckyard-forensicgenetics-<version>.tar.gz`, with the sources under
`deckyard-forensicgenetics-<version>/src/`. Details:
[Installing the deck](https://deckyard.github.io/deck-forensicgenetics/getting-started/installation.html).

## Example

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;

def db as forensicgenetics.FrequencyDb init forensicgenetics.withLocus(
    forensicgenetics.db("example"), "D3S1358", 1034,
    {"15": 0.2456, "16": 0.2314, "17": 0.2007, "18": 0.3223});
def suspect as forensicgenetics.Profile init forensicgenetics.profile("S1", [
    forensicgenetics.genotype("D3S1358", "15", "16")
]);

def rmp as float init forensicgenetics.randomMatchProbability($db, $suspect,
    forensicgenetics.THETA_GENERAL);
def lr as float init forensicgenetics.matchLikelihoodRatio($db, $suspect,
    forensicgenetics.THETA_GENERAL);
```

## Documentation

<https://deckyard.github.io/deck-forensicgenetics/> - concepts, casework
guides, the formulae with their sources, the API reference, and the
[limitations](https://deckyard.github.io/deck-forensicgenetics/project/limitations.html)
to read before quoting a number.

## Not a validated forensic tool

The deck implements published estimators and is tested against them, but it
has not been through the validation any jurisdiction requires of casework
software. Use it for research, teaching and cross-checking, and validate
independently before a number informs a real conclusion. The frequency table
the examples build is synthetic.

**Not for legal or court use**: not for a prosecution, a defence, a paternity
ruling, an immigration decision, or any other legal determination. Those
require validated, accredited software operated under an accredited process.

**No warranty, no liability.** Provided as is, without warranty of any kind,
express or implied; you use it at your own risk. The author accepts no
liability for any claim, damage, loss or other consequence arising from this
software or any number it produces. Checking that a result is correct,
appropriate and lawful for your purpose is yours alone. This restates, and
does not narrow or replace, the disclaimers in the [licence](LICENSE); where
they differ, the licence governs.

Full text: [Not a validated forensic tool](https://deckyard.github.io/deck-forensicgenetics/#not-a-validated-forensic-tool).

## Licence

LGPL-3.0-only, see [LICENSE](LICENSE). Copyright (C) 2026 mplx &lt;mplx@jennifer-lang.dev&gt;.
