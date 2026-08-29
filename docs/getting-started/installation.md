# Installing the deck

## Requirements

Jennifer 0.25.0 or newer. Every source file declares it, and the interpreter
enforces the floor at read time - an older build refuses to load the deck
rather than failing somewhere inside it:

```jennifer
# pragma-jennifer-version: >=0.25.0
```

Nothing else. The deck has no deck dependencies and no Go, so it runs on both
the default `jennifer` binary and the embeddable `jennifer-tiny`.

## Into a vendor tree

[`jvc`](https://github.com/jennifer-language/jvc), the Jennifer package
manager, has not cut a release yet. Until it does, the deck installs by being
copied into a project-local `vendor/` tree - the `@scope/package` resolver has
shipped, so a hand-populated tree resolves exactly as a managed one will. This
deck's manifest already follows the [jvc deck
specification](https://github.com/jennifer-language/jvc/blob/main/docs/deck-spec.md),
so the layout below is the one jvc will produce:

```
your-project/
  app.j
  vendor/
    jennifer/
      forensicgenetics/          <- the contents of this deck's src/, flattened
        forensicgenetics.j
        kinship.j
        pedigree.j
        lineage.j
```

Only `src/` is vendored: the manifest stays in the deck repository, and the
`*_test.j` overlays are not installed.

The installed tree is **flat**: the entry file must sit at its root, because
`@jennifer/forensicgenetics/` appends the package-named entry to the package
directory. In the repository, and in the release tarball, those files live in
`src/`; an install copies their contents up one level.

The release tarball is shaped the way the spec requires - a single wrapping
directory holding `src/`, the manifest, the README and the licence - so
unpacking it is two steps rather than one:

```sh
mkdir -p vendor/jennifer/forensicgenetics
tar -xzf forensicgenetics-0.1.0.tar.gz
cp forensicgenetics-0.1.0/src/*.j vendor/jennifer/forensicgenetics/
rm vendor/jennifer/forensicgenetics/*_test.j
```

The overlays ship in the tarball so an install-time check can run them; they
are removed here because a consumer's tree carries library code only.

From a checkout it is the same copy, without the download:

```sh
mkdir -p vendor/jennifer/forensicgenetics
git clone https://github.com/jennifer-language/deck-forensicgenetics /tmp/deck
cp /tmp/deck/src/*.j vendor/jennifer/forensicgenetics/
rm vendor/jennifer/forensicgenetics/*_test.j
```

Now the imports resolve with no path:

```jennifer
import "@jennifer/forensicgenetics/" as forensics;
import "@jennifer/forensicgenetics/kinship.j" as kinship;
```

The trailing `/` on the first one is not a typo. `@jennifer/forensicgenetics/`
expands to the package-named entry file, `forensicgenetics/forensicgenetics.j`. A sibling
module is named outright.

### Where the vendor root comes from

The resolver looks for `--vendor`, then `JENNIFER_VENDOR`, then the nearest
`vendor/` directory above the program being run. Keeping `vendor/` beside
`app.j` is enough; nothing is installed globally and each project owns its own
copy.

## By relative path

For a checkout you are working *on* rather than *with*, plain relative imports
avoid the vendor tree entirely:

```jennifer
import "./src/forensicgenetics.j" as forensics;
import "./src/kinship.j" as kinship;
```

This is what the bundled examples use, which is why they run straight out of a
clone. Both forms produce the same types: module struct identity is keyed by
the file's canonical path, so a `Profile` built through one import is the same
type as a `Profile` built through the other.

## Checking the install

From the deck root:

```sh
for f in src/*_test.j; do jennifer test "$f"; done
```

All four overlays should report `0 failed`. Then run an example, which
exercises the whole calculation chain against the bundled synthetic table.
There is one per topic:

```sh
jennifer run examples/matchprobability.j   # RMP, LR, and a mixture's CPI / CPE
jennifer run examples/paternity.j          # trio paternity across a mutation
jennifer run examples/pedigree.j           # deficient paternity through two grandparents
jennifer run examples/lineage.j            # mtDNA nomenclature and Y-STR comparison
```

If the tests pass but an example cannot find `testdata/frequencies.xml`, you
are not in the deck root - the examples read that path relative to the working
directory.

## What is in the repository

| Path | What it is |
| --- | --- |
| `src/*.j` | the deck: `forensicgenetics.j`, `kinship.j`, `pedigree.j`, `lineage.j` |
| `src/*_test.j` | white-box test overlays, run with `jennifer test` |
| `examples/` | four runnable programs, one per topic |
| `examples/lib/exampledb.j` | a **synthetic** frequency table the examples build |
| `docs/` | this documentation |
| `LICENSE` | the LGPL-3.0 licence text, with the GPL-3.0 it incorporates |
| `deck.json` | the deck manifest |

The `.j` files in `src/` are the whole deliverable. Everything else supports
them.
