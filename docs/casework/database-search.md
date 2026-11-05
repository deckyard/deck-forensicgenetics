# Database searching

Comparing one profile with many: a direct search for the same person, a
familial search for a relative. Both need the two sides to name loci the same
way and to be compared over a known set of loci; `panels.j` handles both.

## Locus names and kits

Profiles and tables from different labs and software spell loci differently.
`panels.normalizeLocus` maps the common variants to the deck's convention:

```jennifer
import "@deckyard/forensicgenetics/panels.j" as panels;

def a as string init panels.normalizeLocus("THO1");     # "TH01"
def b as string init panels.normalizeLocus("PENTA_E");  # "Penta E"
def c as string init panels.normalizeLocus("VWA");      # "vWA"
```

`normalizeProfile` and `normalizeDb` rename every locus of a profile or table.

A stain typed with one kit and a reference typed with another overlap only
partly. `sharedLoci` gives the overlap; `restrictProfile` and `restrictDb`
narrow to a panel:

```jennifer
def overlap as panels.Panel init panels.intersect(
    panels.panel("Identifiler Plus"), panels.panel("ESS"));
# seven loci: the original (pre-2009) ESS
def stain as forensicgenetics.Profile init panels.restrictProfile($raw, $overlap);
def gaps as list of string init panels.missingInDb($db, $overlap);
```

Check `missingInDb` before a calculation: `randomMatchProbability` and the
kinship LRs throw on an uncovered locus rather than drop it.

## Direct search

```jennifer
import "@deckyard/forensicgenetics/search.j" as search;

def opts as search.SearchOptions init search.searchOptions(10, 1, true);
def ranked as list of search.Comparison init search.search($stain, $database, $opts);
def leads as list of search.Comparison init search.hits($ranked);
```

`searchOptions(minLoci, maxMismatches, allowDropout)` has no default; the
values are a matter of database policy and query quality. Each `Comparison`
reports the loci compared, identical loci, drop-out loci, the mismatching loci
by name, the IBS state per locus and the verdict.

- `minLoci` guards against a "match" on two or three loci.
- `maxMismatches` tolerates typing errors and mutation (for a relative's
  reference) at a stated number of loci.
- `allowDropout` treats a homozygous stain locus `a/a` against a reference
  `a/b` as possible drop-out of `b`. A stain allele the reference lacks is
  always a mismatch.

A partial candidate is a lead for re-examination, not a match; report it with
the RMP over the loci that do match, computed in the usual way.

## Familial search

```jennifer
import "@deckyard/forensicgenetics/kinship.j" as kinship;

def kin as list of search.FamilialHit init search.familialSearch(
    $stain, $database, kinship.FULL_SIBLINGS, $db, 15);
```

Candidates are ranked by `kinship.kinshipLR` for the chosen relationship
against unrelated (see [Kinship](kinship.md)). In a database of thousands,
unrelated people reach high sibling LRs by chance, so the ranking yields a
shortlist for Y-STR or further testing, not an identification. Half siblings,
grandparents and avuncular relatives give identical LRs.

`examples/search.j` runs the full sequence on 200 synthetic references typed
at eight loci. The donor ranks first in both searches; the planted full
sibling ranks 19th, behind unrelated references that share alleles by chance.
Eight loci are too few for a familial search.
