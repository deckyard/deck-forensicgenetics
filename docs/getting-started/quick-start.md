# Quick start

The snippets assume the deck is installed under `vendor/` (see
[Installing the deck](installation.md)).

```jennifer
import "@deckyard/forensicgenetics/" as forensicgenetics;
import "@deckyard/forensicgenetics/kinship.j" as kinship;
import "@deckyard/forensicgenetics/pedigree.j" as pedigree;
```

## 1. A frequency table

```jennifer
def db as forensicgenetics.FrequencyDb init forensicgenetics.withLocus(
    forensicgenetics.db("example"), "TH01", 1000,
    {"6": 0.25, "7": 0.15, "8": 0.30, "9.3": 0.30});
```

`1000` is N, the number of individuals typed at the locus. It sets the
`5/(2N)` floor for rare and unlisted alleles. A reader deck converts a
published table instead, for example
[`@deckyard/strider`](https://github.com/deckyard/deck-strider-frequency-converter)
for STRidER XML. See [Reference frequency data](../concepts/frequency-data.md).

## 2. A match probability

```jennifer
def suspect as forensicgenetics.Profile init forensicgenetics.profile("S1", [
    forensicgenetics.genotype("D3S1358", "15", "16"),
    forensicgenetics.genotype("vWA", "17", "17"),
    forensicgenetics.genotype("TH01", "6", "9.3")
]);

def theta as float init forensicgenetics.THETA_GENERAL;
def rmp as float init forensicgenetics.randomMatchProbability($db, $suspect, $theta);
def lr as float init forensicgenetics.matchLikelihoodRatio($db, $suspect, $theta);
```

Allele labels are strings, and the pair is unordered. The table must cover
every locus of the profile; this one does not, so the call throws (see
[Errors](#errors)). For a long panel use `log10RandomMatchProbability`.

`jennifer run examples/matchprobability.j` runs the full version.

## 3. A paternity trio

```jennifer
def cpi as float init kinship.paternityIndex($mother, $child, $allegedFather, $db);
def w as float init kinship.probabilityOfPaternity($cpi);
```

`paternityIndex` uses `kinship.DEFAULT_MUTATION`, so one locus where the
alleged father cannot have transmitted the child's paternal allele lowers the
index instead of zeroing it. List those loci, or compute without mutation:

```jennifer
def bad as list of string init kinship.inconsistentLoci($mother, $child,
    $allegedFather, $db);
def strict as float init kinship.paternityIndexWith($mother, $child,
    $allegedFather, $db, kinship.NO_MUTATION);
```

`jennifer run examples/paternity.j` runs a trio with one inconsistent locus.

## 4. A relationship

```jennifer
def lr as float init kinship.kinshipLR($a, $b, kinship.FULL_SIBLINGS, $db);
def lr2 as float init kinship.kinshipLRAgainst($a, $b,
    kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, $db);
```

The first tests full siblings against unrelated, the second against half
siblings.

## 5. A pedigree

State both hypotheses as pedigrees and take their ratio. Here the alleged
father `FA` is unavailable and the case rests on his parents:

```jennifer
def h1 as pedigree.Pedigree init pedigree.pedigree("son of GF x GM");
$h1 = pedigree.withFounder($h1, "GF", "m");
$h1 = pedigree.withFounder($h1, "GM", "f");
$h1 = pedigree.withFounder($h1, "MO", "f");
$h1 = pedigree.withChild($h1, "FA", "m", "GF", "GM");   # untyped
$h1 = pedigree.withChild($h1, "CH", "m", "FA", "MO");

def h2 as pedigree.Pedigree init pedigree.pedigree("unrelated man");
$h2 = pedigree.withFounder($h2, "GF", "m");
$h2 = pedigree.withFounder($h2, "GM", "f");
$h2 = pedigree.withFounder($h2, "MO", "f");
$h2 = pedigree.withFounder($h2, "UN", "m");
$h2 = pedigree.withChild($h2, "CH", "m", "UN", "MO");

def lr as float init pedigree.pedigreeLR($h1, $h2, $profiles, $db);
```

`$profiles` is a `map of string to Profile` keyed by member id. `FA` has no
profile, so peeling sums over his genotypes. `jennifer run examples/pedigree.j`
runs this case.

## Errors

Invalid input throws; nothing returns a quiet wrong answer. Each module has its
own `kind`: `"forensicgenetics"`, `"kinship"`, `"pedigree"`, `"lineage"`.

```jennifer
try {
    def rmp as float init forensicgenetics.randomMatchProbability($db, $p, 0.01);
} catch (e) {
    io.printf("%s: %s\n", $e.kind, $e.message);
}
```

`randomMatchProbability` throws on a profile locus the table does not cover,
rather than skipping it and inflating the RMP. To score the shared loci only:

```jennifer
def narrowed as forensicgenetics.Profile init forensicgenetics.restrict($p,
    forensicgenetics.commonLoci($db, $p));
```
