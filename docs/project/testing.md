# Testing and contributing

## Running the tests

Each module has a co-located white-box overlay. `jennifer test` loads the
module beside it, so an overlay can reach private helpers as well as the
exported surface.

```sh
for f in src/*_test.j; do jennifer test "$f"; done
```

One at a time, with coverage:

```sh
jennifer test src/forensicgenetics_test.j --coverage
jennifer test src/pedigree_test.j --filter='.*Peeling.*'
jennifer test src/kinship_test.j --format=tap
```

Current state: 151 tests across four overlays, all passing, with 92% statement
coverage on the core module - the uncovered lines are almost entirely
`fail(...)` branches for input that the type system already makes hard to
produce.

## How the expectations are chosen

The rule the suite follows is that **no expected value is captured from this
implementation**. Every number is derived independently:

- **NRC II genotype frequencies** are computed by hand from recommendation 4.2
  and written into the test as literals, with the arithmetic in a comment.
- **Kinship likelihood ratios** are checked against the textbook closed forms -
  `1/(2p)` for parent/child, `(1 + 1/p)²/4` for two identical homozygotes,
  `1/4 + (1 + pa + pb)/(8 pa pb)` for two identical heterozygotes.
- **Clopper-Pearson bounds** are computed by bisecting the exact binomial CDF
  outside Jennifer, and the `k = 0` case is additionally checked against its
  closed form `1 - α^(1/N)`.

A regression in a formula therefore fails the suite rather than being blessed
by it.

### The fixtures are deliberately tiny

`src/kinship_test.j` and `src/pedigree_test.j` use one locus, `"L"`, with alleles
`10`/`11`/`12`/`13` at frequencies `.1`/`.2`/`.3`/`.4` and `N = 10000`. The
large `N` puts the `5/(2N)` floor at 0.00025 so it never interferes, and the
round frequencies make every expectation hand-derivable in a line of
arithmetic.

## Cross-validation between modules

The strongest checks in the suite are the ones where two modules that share no
code must agree.

`pedigree.j` computes likelihoods by Elston-Stewart peeling - summing over
untyped members' genotypes by passing vectors along the pedigree tree.
`kinship.j` computes the same quantities from closed-form IBD expressions. The
two derivations have nothing in common beyond the allele frequencies, so
agreement is real evidence rather than a tautology.

`src/pedigree_test.j` asserts that a pedigree LR reproduces the closed form for:

| Hypothesis pair | Checked against |
| --- | --- |
| siblings vs unrelated | `kinship.FULL_SIBLINGS` |
| half-siblings vs unrelated | `kinship.HALF_SIBLINGS` |
| parent/child vs unrelated | `kinship.PARENT_CHILD` |
| grandparent through an untyped parent | `kinship.GRANDPARENT` |
| aunt through an untyped father | `kinship.AVUNCULAR` |
| trio through an unrelated-man alternative | `kinship.paternityIndexWith(..., NO_MUTATION)` |

all to within a relative tolerance of 1e-9.

### Other invariants worth knowing

- **Symmetry.** `kinshipLRAt(a, b, ...) == kinshipLRAt(b, a, ...)`. Easy to get
  wrong in the IBD-1 term, so it is asserted rather than assumed.
- **Normalisation.** Summing `jointGenotypeProbability` over every genotype the
  second person could have must give back `P(first person's genotype)` - the
  check that the IBD-1 term is a proper conditional distribution.
- **Mutation totals.** Summing `mutationProbability` over every non-zero step
  in both directions must give back exactly `rate`.
- **Log/product agreement.** `10^log10RandomMatchProbability` must equal
  `randomMatchProbability`.
- **Banded vs full-width alignment.** `align` computes only a band around the
  diagonal for speed. That shortcut is checked the only way worth checking in a
  module whose job is producing a haplotype designation: against a full-width
  pass, on randomised sequences with substitutions and indels, demanding
  byte-identical output - including cases that force the band to widen.

## Style and CI gates

The pipeline runs, and a change must pass:

```sh
jennifer fmt -w src/*.j examples/*.j    # format in place
git diff --exit-code                    # ...and nothing should have changed
jennifer lint src/*.j examples/*.j  # exits 1 on any finding, info level included
```

(Newer interpreters also have `jennifer fmt -l`, which lists unformatted files
directly. The pipeline uses the format-and-diff form because it works on every
version and shows exactly which lines were wrong.)

`jennifer lint` is stricter than it looks - a 100-column limit, a block-nesting
limit of 4, no unused imports, and no method calls inside string interpolation
slots. Run `jennifer fmt -w src/*.j examples/*.j` before committing.

The four examples are also executed in CI, so a change that breaks the
end-to-end path fails even if every unit test still passes.

## Adding to the deck

Where a new estimator goes:

| It works on | Module |
| --- | --- |
| a profile and allele frequencies | `forensicgenetics.j` |
| a pair of people and a relationship | `kinship.j` |
| three or more people joined by parentage | `pedigree.j` |
| a haploid, non-recombining marker | `lineage.j` |

The dependency direction is `forensicgenetics.j ← {kinship.j, pedigree.j}`,
with `lineage.j` an independent leaf and the core importing nothing. Keep it
that way: a reader deck imports the core to produce a `FrequencyDb`, and the
core must never learn about a reader.

Conventions the existing code follows:

- **Copy-returning builders.** Nothing mutates its argument; value semantics do
  the rest.
- **Throw rather than guess.** A missing locus, a mismatched pair, an
  impossible pedigree - each gets a message naming the specific thing that was
  wrong. Silently dropping a locus or returning zero where no answer exists is
  the failure mode this deck is most careful about.
- **One error `kind` per module**, so a caller can tell a bad frequency table
  from a bad pedigree.
- **A docblock on every exported name**, with `@param`, `@return` and `@throws`.
- **A test with a hand-derived expectation**, not a captured one.

## The documentation

This site is built with [Grimoire](https://grimoire.jennifer-lang.dev/) from
the Markdown in `docs/`, with `docs/SUMMARY.md` as the outline and
`grimoire.toml` as the configuration.

```sh
grimoire build            # writes site/
grimoire serve --watch    # preview on :8080
```

or, with nothing installed but a container runtime:

```sh
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" \
    ghcr.io/jennifer-language/grimoire build
```

A `SUMMARY.md` entry naming a file that does not exist makes the build exit 1,
so a broken link in the outline fails CI.
