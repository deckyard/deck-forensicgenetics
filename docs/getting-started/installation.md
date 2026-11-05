# Installing the deck

## Requirements

Jennifer 0.26.0 or newer. Every source file declares the floor, and an older
interpreter refuses to load the deck:

```jennifer
# pragma-jennifer-version: >=0.26.0
```

The deck has no deck dependencies and runs on both `jennifer` and
`jennifer-tiny`.

## Into a vendor tree

The manifest follows the [jvc deck
specification](https://github.com/jennifer-language/jvc/blob/main/docs/deck-spec.md):
a consumer vendors the deck's `src/` directory, without the `*_test.j`
overlays, to `vendor/deckyard/forensicgenetics/`.
[`jvc`](https://github.com/jennifer-language/jvc) does that:

```sh
jvc add @deckyard/forensicgenetics ">=0.4.0"
jvc install
```

A tree copied by hand resolves the same way:

```
your-project/
  app.j
  vendor/
    deckyard/
      forensicgenetics/
        forensicgenetics.j
        kinship.j
        pedigree.j
        ...                  # the other eight modules
```

The entry file must sit at the root of `vendor/deckyard/forensicgenetics/`.

From a release archive, `deckyard-forensicgenetics-<version>.tar.gz`, which
holds one directory `deckyard-forensicgenetics-<version>/` with `src/`, the
manifest, `README.md` and `LICENSE`:

```sh
v=0.4.0
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

Every overlay reports `0 failed`. One example per topic runs against the
synthetic table in `examples/lib/exampledb.j`:

```sh
jennifer run examples/matchprobability.j   # RMP, LR, and a mixture's CPI / CPE
jennifer run examples/paternity.j          # trio paternity across a mutation
jennifer run examples/pedigree.j           # deficient paternity through two grandparents
jennifer run examples/popstats.j           # forensic-efficiency table, HWE permutation test
jennifer run examples/lineage.j            # mtDNA nomenclature and Y-STR comparison
jennifer run examples/mixture.j            # mixture LR with drop-out and drop-in
jennifer run examples/search.j             # direct and familial database search
jennifer run examples/reporting.j          # size-bias LR, posterior bound, ENFSI verbal scale
jennifer run examples/mutation.j           # Familias mutation models, silent alleles
jennifer run examples/xstr.j               # X-STR kinship and daughter paternity
jennifer run examples/popstructure.j       # F-statistics and LD tests
```

## Repository layout

| Path | What it is |
| --- | --- |
| `src/*.j` | the deck: the core `forensicgenetics.j` and ten sibling modules |
| `src/*_test.j` | white-box test overlays, run with `jennifer test` |
| `examples/` | eleven runnable programs, one per topic |
| `examples/lib/exampledb.j` | the synthetic frequency table the examples build |
| `docs/` | this documentation |
| `LICENSE` | the LGPL-3.0 text, with the GPL-3.0 it incorporates |
| `deck.json` | the deck manifest |
