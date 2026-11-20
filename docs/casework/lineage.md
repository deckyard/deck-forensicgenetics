# Lineage markers

mtDNA and the Y chromosome are haploid, non-recombining and inherited from one
parent. The markers travel as one haplotype, so there is no Hardy-Weinberg
genotype frequency and no product across loci. Everyone in the maternal (mtDNA)
or paternal (Y) line shares the haplotype: a lineage marker can exclude a
person but cannot single one out. The match probability is a count: `k`
sightings among `N` database profiles.

## Counting estimates

```jennifer
lineage.haplotypeFrequency(2, 5000);           # 0.0004, the point estimate k/N
lineage.lineageMatchProbability(2, 5000);      # 95% upper bound
lineage.haplotypeUpperBound(2, 5000, 0.99);    # another confidence level
lineage.lineageLikelihoodRatio(2, 5000);       # 1 / the 95% bound
```

Do not report `k/N`: at `k = 0` it is zero, which no finite database supports.
Report the Clopper-Pearson exact-binomial one-sided upper bound: the largest
frequency at which `k` sightings in `N` are still plausible at the stated
confidence. At `k = 0` it is `1 - α^(1/N)` (the "rule of three" at 95%); at
`k = N` it is exactly 1.0, and the LR is 1. `N` must be positive, `k` in
`[0, N]`, and the confidence in `(0.5, 1)`.

`examples/lineage.j`:

```
  k      N          k/N   95% upper bound       LR
  0   5000   0.00000000        0.00059897   1669.5
  1   5000   0.00020000        0.00094842   1054.4
  2   5000   0.00040000        0.00125862    794.5
  5   5000   0.00100000        0.00210145    475.9
```

A haplotype never seen in 5000 profiles is worth an LR of about 1670.

The deck bundles and queries no database. Take `k` and `N` from the same search
of whichever database you used: [EMPOP](https://empop.online/) for mtDNA,
[YHRD](https://yhrd.org/) for Y-STR, or an in-house set.

## mtDNA difference nomenclature

An mtDNA haplotype is the list of differences from the rCRS over the sequenced
ranges:

```jennifer
def h as lineage.MtHaplotype init lineage.mtHaplotype("evidence",
    ["16024-16365", "73-340"], "16189C 16223T 263G 315.1C");
```

| Form | Meaning |
| --- | --- |
| `16189C` | substitution at 16189 to C |
| `A263G` | the same form with the reference base, which is dropped |
| `315.1C` | first base inserted after 315 |
| `249d` | deletion of 249 (`249D`, `249DEL` and `249-` also parse) |
| `16093Y` | an IUPAC ambiguity code, for a heteroplasmic position |

Parsing is case-insensitive. Positions are 1-based rCRS coordinates in
`[1, 999999]`; the bound keeps the fixed-width sort key correct and is far
above the 16569-base genome. Insertion indices run from 1 to 9999. `0.1G`
(position 0) is an insertion before position 1, the form `align` emits for
one before the first base of a region starting at 1. `D` after a position is
always a deletion; the IUPAC code D (A/G/T) cannot be written. `mtHaplotype`
throws when two designations name one position and insertion index
(`263G 263G`, `263G 263A`).

```jennifer
lineage.parseDifference("315.1C");    # Difference{position: 315, insertion: 1, ...}
lineage.formatDifference($d);         # "315.1C"
lineage.mtDesignations($h);           # ["263G", "315.1C", "16189C", "16223T"]
```

Differences sort by position, then insertion index: `315` < `315.1` < `315.2`
< `316`.

### Comparing haplotypes

```jennifer
lineage.mtMatches($a, $b);        # same differences?
lineage.mtDistance($a, $b);       # size of the symmetric difference
lineage.mtRangesAgree($a, $b);    # same sequenced ranges, in any order?
```

`mtMatches` ignores the ranges. An HVR-I-only haplotype and a whole control
region one with the same listed differences are not a like-for-like match, so
check `mtRangesAgree` as well. In `mtDistance` a designation only one side
carries counts once; a position both carry with different bases counts twice.

### From a sequence

```jennifer
def diffs as list of lineage.Difference init lineage.align($reference, $sample, 16024);
def h as lineage.MtHaplotype init lineage.mtHaplotypeFrom("evidence",
    $reference, $sample, 16024);
```

`align` is a global Needleman-Wunsch alignment (match +1, mismatch -1,
gap -2) with each indel run slid as a block to the 3' end of any homopolymer
or tandem-repeat run, so an extra C in the 311-315 stretch reads as `315.1C`,
not `311.1C`, and an extra AC in `GACACT` as `5.1A 5.2C`. Pass the rCRS
(NC_012920.1) region you sequenced, with its first position as `refStart`; no
reference is bundled. `mtHaplotypeFrom` records the range
`refStart`-`(refStart + length - 1)`; both throw when that end passes
999999, the largest position a designation may name.

- EMPOP's phylogenetic alignment rules are not applied. A complex indel may
  be placed differently from EMPOP; check it before quoting.
- Only a band around the diagonal is computed, widened until its score beats
  every path that leaves it, so the result is the full-matrix optimum.
  Near-identical sequences cost roughly linear time; divergent ones approach
  the full O(nm). Regions over `MAX_ALIGN_LENGTH` (5000 bases), and band passes
  over `MAX_ALIGN_CELLS` (1,250,000 cells), throw. Timings are in
  [Limitations](../project/limitations.md#alignment-cost-grows-with-the-region-and-with-divergence).

## EMPOP's EMP files

EMPOP takes and exports mtDNA population data as EMP files: tab-delimited
text, `#` comments, a `#!` line naming the sequenced range(s) for the
haplotypes below it, and one haplotype per line (sample ID, haplogroup,
count, then one difference per cell):

```
# Example population, 3 samples
#! 16024-576
S1	H1c	1	16519C	263G	315.1C	477C
S2	H5	1	16304C	263G	315.1C	456T
S3	H2a2a1	1	
```

The `emp` module reads such a file into `lineage.MtHaplotype` values and
writes them back. It does no I/O; pass the file's text:

```jennifer
import "@deckyard/forensicgenetics/emp.j" as emp;
import "@deckyard/forensicgenetics/lineage.j" as lineage;

def db as emp.EmpFile init emp.readEmp(fs.readString("population.emp"));
def evidence as lineage.MtHaplotype init lineage.mtHaplotype("evidence",
    ["16024-576"], "16304C 263G 315.1C 456T");
def k as int init emp.empCount($db, $evidence);
def p as float init lineage.lineageMatchProbability($k, emp.empSampleSize($db));
```

`empSampleSize` sums the counts (an entry with count 0 is excluded), and
`empCount` sums those of the entries with the same differences over the same
ranges. A haplotype typed over other ranges than the file's is not counted:
compare like with like.

`readEmp` is strict and names the line of what it refuses: a malformed range,
sample ID or count, a repeated sample, a difference `lineage` cannot parse, a
position outside the ranges or past 16569, an insertion `315.2C` without
`315.1C`. Those are the checks of EMPOP's EMPCheck that need no reference
sequence. To submit your own data, build entries and write the file:

```jennifer
def h as lineage.MtHaplotype init lineage.mtHaplotype("S4", ["16024-576"],
    "16189C 73G 263G 315.1C");
def file as emp.EmpFile init emp.EmpFile{
    comments: ["Example population, 1 sample"],
    entries: [emp.empEntry($h, "?", 1)]
};
def text as string init emp.writeEmp($file);
```

`writeEmp` lists the differences in range order, so `16189C` comes before
`73G` in a `16024-576` control region, as EMPOP writes them.

## Y-STR haplotypes

A repeat count per locus. Multi-copy loci (DYS385a/b, DYS389I/II reported
together) are written `"11,14"` or `"11/14"`:

```jennifer
def h as lineage.YHaplotype init lineage.yHaplotype("suspect", {
    "DYS19": "14",
    "DYS385": "11,14",
    "DYS389I": "13",
    "DYS390": "24"
});

lineage.yMatches($a, $b);            # exact at every shared locus
lineage.yStepDistance($a, $b);       # summed repeat-count differences
lineage.yMatchesWithin($a, $b, 1.0); # within one step
lineage.ySharedLoci($a, $b);         # the loci compared
```

Only shared loci are compared, so different panels compare over their overlap;
no shared locus throws. Multi-copy loci are sorted first, so `"11,14"` equals
`"14,11"`. An allele that is not a repeat count (1-6 digits, optionally `.`
and 1-6 more) throws, at construction and on comparison.
Step distances are exact (decimal): `17.2` to `18.2` is one step, while `17.2`
and `18` differ in their fractional part and no number of whole-repeat steps
joins them. `yStepDistance` throws there and on different copy counts at a
shared locus; `yMatchesWithin` returns false, as `yMatches` does.

`yMatchesWithin(a, b, 1.0)` keeps a paternal line together across one
single-step mutation, which `yMatches` would split:

```
  vs paternal relative    exact = false  steps = 1.0  within one step = true
  vs unrelated man        exact = false  steps = 8.0  within one step = false
```

## What a lineage match is worth

A match says the sample came from someone in this maternal or paternal line,
which includes every brother, son and maternal cousin. The LR reflects how
common the haplotype is, not how unusual the person is. Combining it with an
autosomal LR is valid only where the markers are independent, which the deck
does not check.
