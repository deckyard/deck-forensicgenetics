# Lineage markers

Mitochondrial DNA and the Y chromosome behave differently from autosomal STRs
in three ways that change the whole statistical treatment. They are **haploid**
- one copy, no allele pair. They are **non-recombining**, so the markers travel
together as a single unit. And they are **uniparentally inherited** - mtDNA
from the mother, the Y from the father.

Consequences:

- There is no Hardy-Weinberg equilibrium to invoke. No `p²`, no `2pq`.
- There is no product rule across loci. The markers are not independent; they
  are one thing.
- Everyone in a maternal (or paternal) line shares the haplotype. A lineage
  marker can exclude an individual, but it can never identify one.

So the match probability is a **direct count**: how often does this haplotype
appear in a reference database?

## The counting estimate

```jennifer
lineage.haplotypeFrequency(2, 5000);   # 0.0004 - the point estimate k/N
```

This is the number *not* to quote. It is unstable for rare haplotypes: a
haplotype seen once in five thousand, or never, gets a point estimate that no
finite database can actually support. At `k = 0` it is zero, which would mean
the haplotype cannot exist.

## The conservative estimate

The **Clopper-Pearson exact-binomial upper bound**: the largest population
frequency that would still make `k` sightings in `N` profiles unsurprising at
the stated confidence.

```jennifer
lineage.lineageMatchProbability(2, 5000);              # 95% upper bound
lineage.haplotypeUpperBound(2, 5000, 0.99);            # another confidence level
lineage.lineageLikelihoodRatio(2, 5000);               # 1 / the bound
```

It is well defined at `k = 0`, where it reduces to `1 - α^(1/N)` - the "rule of
three" at 95% - and it is exactly 1.0 at `k = N`.

```
k       N       k/N          95% upper bound   LR
0       5000    0.00000000   0.00059897        1669.5
1       5000    0.00020000   0.00094842        1054.4
2       5000    0.00040000   0.00125862        794.5
5       5000    0.00100000   0.00210145        475.9
```

The `k = 0` row is the whole argument for the bound: a haplotype never seen in
five thousand people is worth about 1670, not infinity.

Internally the one-sided bound at confidence `c` is the upper end of a
two-sided Clopper-Pearson interval at level `2c - 1`, the standard identity.

## Database-independent by design

The deck bundles no reference database and does not query one. It supplies the
estimator; you supply `k` and `N` from whichever database you actually
searched - [EMPOP](https://empop.online/) for mtDNA,
[YHRD](https://yhrd.org/) for Y-STR, or an in-house set.

That keeps the statistics honest about their source. The `N` you pass must be
the size of the database you searched, and the `k` must come from the same
search.

## mtDNA: difference nomenclature

An mtDNA haplotype is not reported as a sequence. It is reported as the list of
**differences from the reference** - the rCRS - which is far shorter and is
what databases index.

```jennifer
def h as lineage.MtHaplotype init lineage.mtHaplotype("evidence",
    ["16024-16365", "73-340"], "16189C 16223T 263G 315.1C");
```

| Form | Meaning |
| --- | --- |
| `16189C` | substitution at 16189 to C |
| `A263G` | the same, with the redundant reference base; it is dropped |
| `315.1C` | first base inserted after 315 |
| `249d` | deletion of 249 (`249D`, `249DEL`, `249-` also parse) |
| `16093Y` | an IUPAC ambiguity code, for a heteroplasmic position |

```jennifer
lineage.parseDifference("315.1C");    # Difference{position: 315, insertion: 1, ...}
lineage.formatDifference($d);         # back to "315.1C"
lineage.mtDesignations($h);           # ["263G", "315.1C", "16189C"]
```

Differences sort by position, then by insertion index, so `315` < `315.1` <
`315.2` < `316`.

### The ranges matter

```jennifer
lineage.mtMatches($a, $b);        # same difference set?
lineage.mtDistance($a, $b);       # how many designations differ
lineage.mtRangesAgree($a, $b);    # were they sequenced over the same regions?
```

`mtMatches` compares difference lists only. Two haplotypes with identical
listed differences are **not** a like-for-like match if one covers only HVR-I
and the other the whole control region - the second was simply looked at more
closely. Check `mtRangesAgree` too; the deck keeps the two questions separate
rather than guessing.

`mtDistance` counts the symmetric difference: a position one carries and the
other lacks counts once, and a position both carry with different bases counts
twice.

### Rendering a sequence

`align` turns a raw sequence into that nomenclature - a global
Needleman-Wunsch alignment (match +1, mismatch -1, gap -2) with indels slid to
the 3' end of any homopolymer run, which is what makes an extra C in the
311-315 stretch read as `315.1C` rather than `311.1C`.

```jennifer
def diffs as list of lineage.Difference init lineage.align($reference, $sample, 16024);
def h as lineage.MtHaplotype init lineage.mtHaplotypeFrom("evidence",
    $reference, $sample, 16024);
```

**No reference is bundled.** Pass the rCRS (NC_012920.1) region you sequenced,
with `refStart` its first position, so the differences come out in rCRS
coordinates.

Two caveats:

- **EMPOP's phylogenetic alignment rules are not applied.** This renders a
  plain optimal alignment. EMPOP shifts indels to the position the mtDNA
  phylogeny prefers, which can differ for a complex indel. Check such a
  designation before quoting it.
- **Cost grows with the region.** Only a band around the diagonal is computed,
  which is where the optimal path lies when a haplotype is compared with its own
  reference, so the usual case is roughly linear: 0.34 s for a 342 bp HVR-I
  region, 1.5 s for the ~1.1 kb control region. Sequences that genuinely diverge
  widen the band, up to the full quadratic matrix. `align` refuses a region
  longer than `MAX_ALIGN_LENGTH` (5000 bases) rather than appearing to hang.

Reference positions are bounded to six digits, far above the 16569-base
mitochondrial genome. That bound is what keeps difference designations sorting
correctly, since they are ordered through a fixed-width key.

## Y-STR haplotypes

A Y haplotype is a repeat count per locus. Multi-copy loci - DYS385a/b, and
DYS389I/II when reported together - hold their alleles comma-separated:

```jennifer
def h as lineage.YHaplotype init lineage.yHaplotype("suspect", {
    "DYS19": "14",
    "DYS385": "11,14",
    "DYS389I": "13",
    "DYS390": "24"
});
```

`"11,14"` and `"14,11"` are the same haplotype; both sides are sorted before
comparison. Differing copy counts at a shared locus are an error rather than a
guess.

```jennifer
lineage.yMatches($a, $b);            # exact at every shared locus
lineage.yStepDistance($a, $b);       # total repeat-unit difference
lineage.yMatchesWithin($a, $b, 1.0); # the single-step comparison
lineage.ySharedLoci($a, $b);         # which loci were compared
```

Only shared loci are compared, so two haplotypes typed on different panels
still compare over their overlap. The **single-step** form is what keeps a
paternal lineage together across one mutation: a relative differing by one
repeat at one locus is still the same lineage, and `yMatches` alone would
wrongly separate them.

```
  vs paternal relative    exact = false  steps = 1.0  within one step = true
  vs unrelated man        exact = false  steps = 8.0  within one step = false
```

## What lineage evidence is worth

A lineage match says *this sample came from someone in this maternal (or
paternal) line*. That line includes every brother, every son, every maternal
cousin. The likelihood ratio reflects how common the haplotype is, not how
unusual the person is - and it does not narrow the line to one member.

Combining a lineage LR with an autosomal one is legitimate only where the
markers are genuinely independent, which needs care beyond what this deck does
for you.

## Related

- [Kinship](kinship.md) - where autosomal markers cannot separate half-sibs,
  uncles and grandparents, lineage markers sometimes can
- [`lineage` reference](../reference/lineage.md)
