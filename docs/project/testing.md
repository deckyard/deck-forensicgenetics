# Testing and contributing

## Running the tests

Each module has a co-located white-box overlay; `jennifer test` loads the
module beside it, so an overlay calls private helpers as well as exports.
`jennifer test` takes one file:

```sh
for f in src/*_test.j; do jennifer test "$f"; done

jennifer test src/forensicgenetics_test.j --coverage
jennifer test src/pedigree_test.j --filter='.*Peeling.*'
jennifer test src/kinship_test.j --format=tap
```

181 tests across four overlays. Statement coverage of `forensicgenetics.j` is
99.2%; the two uncovered lines are the zero-RMP guards in the likelihood-ratio
functions, which only a product underflowing to zero reaches.

## How expectations are chosen

No expected value is captured from this implementation. Every number is
derived independently:

- NRC II genotype frequencies by hand from recommendation 4.2, written as
  literals with the arithmetic in a comment.
- Kinship LRs from the textbook closed forms: `1/(2p)` for parent/child,
  `(1 + 1/p)²/4` for two identical homozygotes,
  `1/4 + (1 + pa + pb)/(8 pa pb)` for two identical heterozygotes.
- Clopper-Pearson bounds by bisecting the exact binomial CDF outside Jennifer,
  and at `k = 0` against the closed form `1 - α^(1/N)`.

A formula regression therefore fails the suite instead of being written into
it.

`src/kinship_test.j` and `src/pedigree_test.j` use one locus `"L"` with alleles
`10`/`11`/`12`/`13` at `.1`/`.2`/`.3`/`.4` and `N = 10000`, which puts the
`5/(2N)` floor at 0.00025, below every allele, and keeps each expectation to a
line of arithmetic.

## Cross-validation between modules

`pedigree.j` computes likelihoods by peeling; `kinship.j` computes the same
quantities from closed-form IBD expressions. They share nothing but the allele
frequencies, so agreement is evidence. `src/pedigree_test.j` asserts that a
pedigree LR reproduces the closed form, to a relative tolerance of 1e-9, for:

| Hypothesis pair | Checked against |
| --- | --- |
| full siblings vs unrelated (two genotype pairs) | `kinship.FULL_SIBLINGS` |
| half siblings vs unrelated | `kinship.HALF_SIBLINGS` |
| parent/child vs unrelated | `kinship.PARENT_CHILD` |
| grandparent through an untyped parent | `kinship.GRANDPARENT` |
| aunt through an untyped father | `kinship.AVUNCULAR` |
| trio against an unrelated-man alternative | `kinship.paternityIndexWith(mother, child, allegedFather, db, kinship.NO_MUTATION)` |

Other asserted invariants:

- Symmetry: `kinshipLRAt(a, b, ...) == kinshipLRAt(b, a, ...)`.
- Normalisation: summing `jointGenotypeProbability` over every genotype of the
  second person gives `P(first person's genotype)`.
- Mutation totals: `mutationProbability` over every non-zero step in both
  directions sums to `rate`.
- `10^log10RandomMatchProbability` equals `randomMatchProbability`.
- Banded alignment equals a full-width pass, byte for byte, on randomised
  sequences with substitutions and indels, including cases that force the band
  to widen.

## CI gates

```sh
tools/check-all.sh
```

runs what CI runs and ends with the list of failures: `jennifer fmt -l`,
`jennifer lint` (fails on a warning; 100-column limit, nesting limit of 4, no
unused imports, no method calls in interpolation), the header, overlay,
coverage and docs guards, the four suites, and every example.
`jennifer fmt -w <files>` fixes formatting.

`tools/check-docs.sh` compares the [Cheatsheet](../reference/cheatsheet.md)
with the code: every signature, field and constant value, and every docblock's
`@param` / `@field` names, written `@param name {type} description`.
`tools/check-coverage.sh` fails when an export is not called from an overlay.

## Adding to the deck

| It works on | Module |
| --- | --- |
| a profile and allele frequencies | `forensicgenetics.j` |
| a pair of people and a relationship | `kinship.j` |
| three or more people joined by parentage | `pedigree.j` |
| a haploid, non-recombining marker | `lineage.j` |

Dependencies run `forensicgenetics.j ← {kinship.j, pedigree.j}`; `lineage.j`
and the core import nothing from the deck. A reader deck imports the core to
produce a `FrequencyDb`; the core must not learn about a reader.

Conventions:

- Builders return copies and never mutate an argument.
- Throw rather than guess, with a message naming what was wrong. Never drop a
  locus silently or return zero where no answer exists.
- One error `kind` per module.
- A docblock on every export, with `@param`, `@return` and `@throws`.
- A test with a hand-derived expectation.

## The documentation

The book is built with [Grimoire](https://grimoire.jennifer-lang.dev/) from
`docs/`, with `docs/SUMMARY.md` as the outline and `grimoire.toml` as the
configuration:

```sh
grimoire build            # writes site/
grimoire serve --watch    # preview on :8080
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" \
    ghcr.io/jennifer-language/grimoire build
```

A `SUMMARY.md` entry naming a missing file fails the build.
