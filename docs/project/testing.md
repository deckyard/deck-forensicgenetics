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

495 tests across eleven overlays. Statement coverage of `forensicgenetics.j` is
99.2%; the two uncovered lines are the zero-RMP guard in `reciprocalLR`, which
only a product underflowing to zero reaches, and the final `return 0` of
`compareAlleles`, unreachable because identical labels return early.

## How expectations are chosen

No expected value is captured from this implementation. Every number is
derived independently:

- NRC II genotype frequencies by hand from recommendation 4.2, written as
  literals with the arithmetic in a comment.
- Kinship LRs from the textbook closed forms: `1/(2p)` for parent/child,
  `(1 + 1/p)²/4` for two identical homozygotes,
  `1/4 + (1 + pa + pb)/(8 pa pb)` for two identical heterozygotes.
- Kinship and paternity with theta from closed forms derived by hand from the
  Balding-Nichols sampling formula at `θ = 0.1`, each factor written out in
  the comment; the `IDENTICAL` LR against `1 / genotypeFrequency` (NRC II 4.2,
  a separate code path); and every `...Theta` function at `θ = 0` against its
  original with exact equality.
- Mixture LRs from hand sums over the unknowns' genotypes (drop-out, drop-in,
  the lumped allele class), the `d = 0`, `c = 0` inclusion-exclusion forms
  `(pa+pb)^(2u) - pa^(2u) - pb^(2u)`, and the identity with
  `matchLikelihoodRatio` at θ 0, 0.01 and 0.03.
- Population statistics from the power sums of `.1/.2/.3/.4` and `.5/.5`,
  written out in the overlay; HWE p-values from the Levene-Haldane
  distribution enumerated by hand for samples of three and four people.
- F-statistics from Weir & Cockerham's `a`, `b`, `c` worked by hand per
  allele on samples of two to five people (fractions written out), the
  identical-populations closed form `θ = -(1 - Ho/(2 He)) / (n - 1)`, and the
  Genepop documentation's worked example `Fmulti.txt` (per-locus and
  multilocus FIS / FST / FIT to the 4 places it prints).
- LD p-values from the hypergeometric distribution of 2 x 2 genotype tables
  by hand, and from exact enumeration of all `n!` re-pairings (n = 7) in the
  overlay, which the Monte Carlo estimate must match within 4 SE.
- Clopper-Pearson bounds by bisecting the exact binomial CDF outside Jennifer,
  and at `k = 0` against the closed form `1 - α^(1/N)`.
- Kit panels from the set relations and locus counts the standards and
  manuals state (GlobalFiler = CODIS 20 + SE33, ESX 17 = ESS + D2S1338,
  D16S539, D19S433, SE33), written out again in the overlay.
- Direct-search counts by hand from genotypes written beside each candidate.
- X-STR LRs and daughter paternity indices by hand from the X-IBD forms
  (`5`, `10`, `5.5`, `14.375` and so on, each with its arithmetic); every
  relationship constant's X coancestry against Gomes et al. (2020), Table 5;
  linkage-group probabilities summed by hand over the female's phases;
  wide groups (up to `MAX_GROUP_HETEROZYGOUS_LOCI` heterozygous loci, every
  haplotype at the floor `f`) against closed forms (full sisters
  `1 / (2 · 4^h · f)`, daughter trio `1 / f`); the
  Argus X-12 groups against the QIAGEN handbook's Table 1.
- Size-bias frequencies by hand from `(x + 2)/(2N + 4)` and `(x + 4)/(2N + 4)`
  on a table with `N = 100`, and with θ against the core's
  `genotypeFrequency` on a table holding the corrected frequencies.
- Posterior quantiles against laws with known quantiles: Beta(2, 2) (CDF
  `3x² - 2x³`), the arcsine law Beta(1/2, 1/2) (`sin²(qπ/2)`), and
  Beta(1/3, 8/3) inverted through `math.regBetaI`. Each tolerance is about
  4-5 Monte Carlo standard deviations of a nearest-rank quantile,
  `sqrt(q (1 - q) / M) / f(x_q)`, worked out in the comment.
- Verbal bands at every edge of the ENFSI (2015) table, on both sides of 1.
- Mutation matrices entry by entry from the FamiliasLocus definitions on the
  `.1/.2/.3/.4` locus (Equal `.001`, Proportional `alpha = .01`, Stepwise
  `.016/7`, ...); PM on two alleles by hand (`d = (2, 2/3)`); trio and duo
  indices with mutation, a maternal mutation and silent alleles from the
  joint genotype sums written out in the overlay, at `θ = 0` and `0.1`;
  AABB counts against their own printed totals and STRBase's rounded
  percentages.

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
| full siblings on an `N = 20` table with an allele at its `5/(2N)` floor | `kinship.FULL_SIBLINGS`, and the closed form `1/4 + 1/(8 p)` |

The agreement holds wherever `pedigree.j` returns a value: a locus whose
observed, floored frequencies sum past 1 throws instead, since rescaling them
would move a floored allele below its floor and break the agreement, and
keeping them would make an untyped founder's prior sum to more than 1.

Other asserted invariants:

- Symmetry: `kinshipLRAt(a, b, ...) == kinshipLRAt(b, a, ...)`.
- Normalisation: summing `jointGenotypeProbability` over every genotype of the
  second person gives `P(first person's genotype)`.
- Mutation totals: `mutationProbability` over every non-zero step in both
  directions sums to `rate`.
- `popstats.matchProbability` equals `Σ P(g)²` over the core's
  `genotypeFrequency` at `θ = 0`.
- `hweExactTest` (Wigginton recurrence) equals a direct Levene-formula sum
  through `lgamma` to 1e-9; `hwePermutationTest` is within 4 standard errors
  of the exact value; the generator reproduces the check value ISO C++ requires of
  `std::minstd_rand` (state 399268537 after 10000 steps from seed 1). The
  hashed count table (forced with a dense limit of 0) gives the dense run's
  p-values and LD hits exactly; with every allele or genotype distinct each
  re-pairing is as probable as the observed table, so p = 1 exactly at 1000
  people. `ldPermutationTests` gives, bit for bit, the p-values of each pair
  encoded on its own and run on one stream.
- `10^log10RandomMatchProbability` equals `randomMatchProbability`.
- `reporting`: the size-bias correction shrinks to the plug-in value as `N`
  grows (relative 1e-6 at `N = 10⁸`); the posterior median of log10 RMP is
  within 3e-4 of the plug-in value at `N = 10⁶`; the posterior interval
  narrows with `N`; gamma deviates have the right mean and variance, and
  splitting a draw across calls changes no deviate; the
  generator reproduces the `std::minstd_rand` check value; a table summing
  above 1 gets posterior mass `2N + 1`.
- `search.familialSearch` reproduces `kinship.kinshipLR` exactly for every
  candidate and relationship tested; a panel-restricted profile scores as
  `forensicgenetics.restrict` over the same loci.
- X-STR: two females score as `kinship.kinshipLRAt` with the same
  coefficients; the daughter duo index without mutation equals the
  father-daughter LR for every genotype pair; the joint probability summed
  over the second person's genotypes is the first's marginal; a one-locus
  linkage group equals the independent model; a group in linkage equilibrium
  gives the independent match probability, and the same mother-daughter LR
  (one meiosis).
- `mutation.j` with no mutation and no silent allele equals
  `kinship.paternityIndexAtTheta`, `kinship.duoPaternityIndexAt` and
  `kinship.kinshipLRAtTheta` exactly (`==`); its trio index equals a
  brute-force enumeration of every ordered parental and child allele pair,
  Familias-style, with silent alleles, sex-specific models and PM. Models
  reduce to each other at boundary values (Stepwise at range 1 and
  Proportional on uniform frequencies are Equal; PM leaves Proportional
  unchanged); every matrix row sums to 1; PM output is stationary with the
  expected rate kept; phenotype probabilities with a silent allele sum to the
  marginal.
- Banded alignment equals a full-width pass, byte for byte, on randomised
  sequences with substitutions and indels, including cases that force the band
  to widen, and equals an independent full-matrix Needleman-Wunsch in the
  overlay (`fullMatrix`) on those and on a 40-base shifted window, where the
  optimum leaves a band whose traceback never touches its edge.

## CI gates

```sh
tools/check-all.sh
```

runs what CI runs and ends with the list of failures: `jennifer fmt -l`,
`jennifer lint` (fails on a warning; 100-column limit, nesting limit of 4, no
unused imports, no method calls in interpolation), the header, overlay,
coverage and docs guards, every suite, and every example.
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
| a mixed stain and contributor hypotheses | `mixture.j` |
| a frequency table or genotype sample as a whole | `popstats.j` |
| locus names and kit panels | `panels.j` |
| one profile against a list of profiles | `search.j` |
| X-chromosomal profiles, sex-specific kinship | `xstr.j` |
| sampling uncertainty or wording of a reported number | `reporting.j` |
| mutation matrices, sex-specific models, silent alleles | `mutation.j` |

Dependencies run `forensicgenetics.j ← {kinship.j, pedigree.j, mixture.j, popstats.j, panels.j}`
`forensicgenetics.j ← reporting.j` and
`{forensicgenetics.j, kinship.j} ← {search.j, xstr.j, mutation.j}`; `lineage.j`
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
