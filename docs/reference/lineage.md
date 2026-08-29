# `lineage`

```jennifer
import "@jennifer/forensicgenetics/lineage.j" as lineage;
```

mtDNA and Y-STR haplotypes. An independent leaf module - it does not depend on
`forensicgenetics` and can be used on its own. Errors are `Error{kind: "lineage"}`.

## Constants

| Name | Type | Value | Meaning |
| --- | --- | --- | --- |
| `DEFAULT_CONFIDENCE` | `float` | `0.95` | the one-sided confidence level for a conservative estimate |
| `MAX_ALIGN_LENGTH` | `int` | `5000` | the longest region `align` accepts, in bases |

## Types

```jennifer
def struct Difference {
    position as int,      # the reference position
    insertion as int,     # 0 for a substitution or deletion; 1, 2, ... for an insertion
    base as string,       # the observed base; "" for a deletion
    kind as string        # "sub", "ins" or "del"
};

def struct MtHaplotype {
    id as string,
    ranges as list of string,            # e.g. ["16024-16365", "73-340"]
    differences as list of Difference    # ascending position order
};

def struct YHaplotype {
    id as string,
    loci as map of string to string      # locus -> allele(s), multi-copy comma-separated
};
```

## Counting estimators

| Signature | Returns |
| --- | --- |
| `haplotypeFrequency(k as int, n as int)` | `float` - the point estimate `k/N` |
| `haplotypeUpperBound(k as int, n as int, confidence as float)` | `float` - one-sided Clopper-Pearson upper bound |
| `lineageMatchProbability(k as int, n as int)` | `float` - the bound at `DEFAULT_CONFIDENCE` |
| `lineageLikelihoodRatio(k as int, n as int)` | `float` - `1 / lineageMatchProbability` |

`N` must be positive and `k` in `[0, N]`; `confidence` must be in `(0.5, 1)`.
The bound is exactly `1.0` at `k == N`. Quote the bound, not `k/N` - see
[Lineage markers](../casework/lineage.md).

## rCRS difference nomenclature

| Signature | Returns |
| --- | --- |
| `parseDifference(text as string)` | `Difference` |
| `formatDifference(d)` | `string` |

Accepted forms: `16189C`, `A263G` (the reference base is dropped), `315.1C`,
and the deletions `249d`, `249D`, `249DEL`, `249-`. IUPAC ambiguity codes parse,
for heteroplasmic positions.

Positions must be in `[1, 999999]` - far above the 16569-base mitochondrial
genome. The bound exists because designations sort through a fixed-width key,
which a wider position would order wrongly.

## mtDNA haplotypes

| Signature | Returns |
| --- | --- |
| `mtHaplotype(id as string, ranges as list of string, text as string)` | `MtHaplotype` - `text` is whitespace-separated designations |
| `mtDesignations(h)` | `list of string` - in position order |
| `mtMatches(a, b)` | `bool` - identical difference sets |
| `mtDistance(a, b)` | `int` - size of the symmetric difference |
| `mtRangesAgree(a, b)` | `bool` - same sequenced regions, in any order |

`mtMatches` compares difference lists only. Check `mtRangesAgree` as well: two
haplotypes with the same listed differences are not a like-for-like match if one
was sequenced over less of the molecule.

## Alignment

| Signature | Returns |
| --- | --- |
| `align(reference as string, sample as string, refStart as int)` | `list of Difference` |
| `mtHaplotypeFrom(id as string, reference as string, sample as string, refStart as int)` | `MtHaplotype` - `align` plus the range |

Global Needleman-Wunsch (match +1, mismatch -1, gap -2), with indels slid to
the 3' end of any homopolymer run. Whitespace in either sequence is ignored; a
non-nucleotide character throws, and so does a region longer than
`MAX_ALIGN_LENGTH`.

Only a band around the diagonal is computed, widening automatically if the
traceback touches its edge, so the result is always the true global optimum and
never an approximation - the test overlay checks that against a full-width pass
on randomised input. For near-identical sequences the cost is roughly linear:
0.34 s for 342 bases, 1.5 s for the ~1.1 kb control region.

**No reference is bundled** - pass the rCRS region you sequenced.
**EMPOP's phylogenetic alignment rules are not applied.**

## Y-STR haplotypes

| Signature | Returns |
| --- | --- |
| `yHaplotype(id as string, loci as map of string to string)` | `YHaplotype` |
| `yLoci(h)` | `list of string` |
| `ySharedLoci(a, b)` | `list of string` - in the first haplotype's order |
| `yMatches(a, b)` | `bool` - exact at every shared locus |
| `yStepDistance(a, b)` | `float` - summed absolute repeat differences |
| `yMatchesWithin(a, b, steps as float)` | `bool` - distance at most `steps` |

Multi-copy loci are written `"11,14"` (or `"11/14"`) and compared as sorted
sets, so writing order does not matter. Differing copy counts at a shared locus
throw rather than being guessed at. Comparing haplotypes that share no locus
throws.

Y-STR alleles must be repeat counts. A non-numeric designation is rejected at
construction and again on comparison, with a `lineage` error naming the allele.

## Errors

`Error{kind: "lineage"}` for: a non-positive `N`, a `k` outside `[0, N]`, a
confidence outside `(0.5, 1)`; a malformed difference designation or a position
outside `[1, 999999]`; an empty sequence, a non-positive `refStart`, a
non-nucleotide character, or a region longer than `MAX_ALIGN_LENGTH` in
`align`; an empty locus name or allele list; a Y-STR allele that is not a
repeat count; no shared locus between two Y haplotypes; differing copy counts
at a shared locus; and a negative step tolerance.
