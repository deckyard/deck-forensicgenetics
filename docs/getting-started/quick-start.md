# Quick start

Four calculations, each self-contained. Run them from the deck root so the
relative paths resolve.

## 1. A frequency table

Every calculation needs reference allele frequencies. The deck reads no
published format itself - a reader deck converts one
([`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter) does
STRidER-style XML):

```jennifer
import "@jennifer/forensicgenetics/" as forensics;
import "@mplx/strider/" as strider;

def db as forensics.FrequencyDb init strider.load($xml, "STRidER AT");
```

or build one directly, which is what the tests and examples do:

```jennifer
def db as forensics.FrequencyDb init forensics.withLocus(
    forensics.db("example"), "TH01", 1000,
    {"6": 0.25, "7": 0.15, "8": 0.30, "9.3": 0.30});
```

The `1000` is N, the number of individuals typed at that locus. It is not
decoration: it sets the `5/(2N)` floor that keeps a rare allele from driving a
match probability to zero. See [Reference frequency data](../concepts/frequency-data.md).

## 2. A match probability

```jennifer
def suspect as forensics.Profile init forensics.profile("S1", [
    forensics.genotype("D3S1358", "15", "16"),
    forensics.genotype("vWA", "17", "17"),
    forensics.genotype("TH01", "6", "9.3")
]);

def theta as float init forensics.THETA_GENERAL;
def rmp as float init forensics.randomMatchProbability($db, $suspect, $theta);
def lr as float init forensics.matchLikelihoodRatio($db, $suspect, $theta);
```

`rmp` is the chance an unrelated person from the reference population carries
this profile; `lr` is its reciprocal, the weight of evidence. For a long panel
use `log10RandomMatchProbability`, which cannot underflow.

Allele designations are **strings**, not numbers - `"9.3"` at TH01 is a real
allele and must not collapse into `"9"`. The pair is unordered:
`genotype("TH01", "9.3", "6")` and `genotype("TH01", "6", "9.3")` are the same
value.

Run the whole thing: `jennifer run examples/matchprobability.j`

## 3. A paternity trio

```jennifer
import "../kinship.j" as kinship;

def cpi as float init kinship.paternityIndex($mother, $child, $allegedFather, $db);
def w as float init kinship.probabilityOfPaternity($cpi);
```

`paternityIndex` uses `DEFAULT_MUTATION`, so a single locus where the child
carries an allele the alleged father cannot have transmitted costs the index
rather than zeroing it. To see which loci those are:

```jennifer
def bad as list of string init kinship.inconsistentLoci($mother, $child,
    $allegedFather, $db);
```

and to compute strictly, with no allowance for mutation:

```jennifer
def strict as float init kinship.paternityIndexWith($mother, $child,
    $allegedFather, $db, kinship.NO_MUTATION);
```

Run it: `jennifer run examples/paternity.j`

## 4. A relationship

```jennifer
def lr as float init kinship.kinshipLR($a, $b, kinship.FULL_SIBLINGS, $db);
```

tests full siblings against unrelated. To test one relationship against
another rather than against unrelatedness:

```jennifer
def lr as float init kinship.kinshipLRAgainst($a, $b,
    kinship.FULL_SIBLINGS, kinship.HALF_SIBLINGS, $db);
```

## 5. A pedigree

When the person you want to test is unavailable, state both hypotheses as
pedigrees and take their ratio. Here the alleged father is dead, and the case
rests on his parents:

```jennifer
import "../pedigree.j" as pedigree;

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

`FA` is never typed; peeling sums over every genotype he could have had.
`profiles` is a `map of string to Profile` keyed by pedigree member id.

Run it: `jennifer run examples/pedigree.j`

## Errors

Everything validates its input and throws rather than returning a quiet wrong
answer. Each module uses its own `kind`:

```jennifer
try {
    def rmp as float init forensics.randomMatchProbability($db, $p, 0.01);
} catch (e) {
    io.printf("%s: %s\n", $e.kind, $e.message);
}
```

`"forensicgenetics"`, `"kinship"`, `"pedigree"`, `"lineage"`.

One rule worth knowing early: `randomMatchProbability` **refuses** a profile
locus the database does not cover, rather than skipping it. A silently dropped
locus inflates the RMP without saying so. Narrow the profile deliberately if
that is what you want:

```jennifer
def narrowed as forensics.Profile init forensics.restrict($p,
    forensics.commonLoci($db, $p));
```
