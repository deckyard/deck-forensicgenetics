# Installing the deck

## Requirements

Jennifer 0.25.0 or newer. Every source file declares the floor, and an older
interpreter refuses to load the deck:

```jennifer
# pragma-jennifer-version: >=0.25.0
```

The deck has no deck dependencies and runs on both `jennifer` and
`jennifer-tiny`.

## Into a vendor tree

The manifest follows the [jvc deck
specification](https://github.com/jennifer-language/jvc/blob/main/docs/deck-spec.md):
a consumer vendors the deck's `src/` directory, without the `*_test.j`
overlays, to `vendor/deckyard/forensicgenetics/`. Until
[`jvc`](https://github.com/jennifer-language/jvc) does that for you, copy the
files by hand; the `@scope/package` resolver reads a hand-built tree the same
way.

```
your-project/
  app.j
  vendor/
    deckyard/
      forensicgenetics/
        forensicgenetics.j
        kinship.j
        pedigree.j
        lineage.j
```

The entry file must sit at the root of `vendor/deckyard/forensicgenetics/`.

From a release archive, `deckyard-forensicgenetics-<version>.tar.gz`, which
holds one directory `deckyard-forensicgenetics-<version>/` with `src/`, the
manifest, `README.md` and `LICENSE`:

```sh
v=0.2.0
mkdir -p vendor/deckyard/forensicgenetics
tar -xzf deckyard-forensicgenetics-$v.tar.gz
cp deckyard-forensicgenetics-$v/src/*.j vendor/deckyard/forensicgenetics/
rm vendor/deckyard/forensicgenetics/*_test.j
```

From a checkout:

```sh
mkdir -p vendor/deckyard/forensicgenetics
git clone https://github.com/deckyard/deck-forensicgenetics /tmp/deck
cp /tmp/deck/src/*.j vendor/deckyard/forensicgenetics/
rm vendor/deckyard/forensicgenetics/*_test.j
```

Then:

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/kinship.j" as kinship;
```

The trailing `/` resolves to the entry file,
`forensicgenetics/forensicgenetics.j`; sibling modules are named outright.

The resolver takes the vendor root from `--vendor`, then `JENNIFER_VENDOR`,
then the nearest `vendor/` directory above the program being run.

## By relative path

Inside a checkout, relative imports need no vendor tree, and the examples use
them:

```jennifer
import "./src/forensicgenetics.j" as forensicgenetics;
import "./src/kinship.j" as kinship;
```

Struct identity is keyed by the file's canonical path, so a `Profile` built
through either import is the same type.

## Checking the install

From the deck root:

```sh
for f in src/*_test.j; do jennifer test "$f"; done
```

All four overlays report `0 failed`. One example per topic runs against the
synthetic table in `examples/lib/exampledb.j`:

```sh
jennifer run examples/matchprobability.j   # RMP, LR, and a mixture's CPI / CPE
jennifer run examples/paternity.j          # trio paternity across a mutation
jennifer run examples/pedigree.j           # deficient paternity through two grandparents
jennifer run examples/lineage.j            # mtDNA nomenclature and Y-STR comparison
```

## Repository layout

| Path | What it is |
| --- | --- |
| `src/*.j` | the deck: `forensicgenetics.j`, `kinship.j`, `pedigree.j`, `lineage.j` |
| `src/*_test.j` | white-box test overlays, run with `jennifer test` |
| `examples/` | four runnable programs, one per topic |
| `examples/lib/exampledb.j` | the synthetic frequency table the examples build |
| `docs/` | this documentation |
| `LICENSE` | the LGPL-3.0 text, with the GPL-3.0 it incorporates |
| `deck.json` | the deck manifest |
