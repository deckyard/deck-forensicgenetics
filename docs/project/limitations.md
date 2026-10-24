# Limitations

## Not a validated forensic tool

The deck implements published estimators and is tested against them, but it has
not been through the validation any jurisdiction requires of software used in
casework. Treat it as a library for research, teaching and cross-checking, and
validate independently before it informs a real conclusion.

### Not for legal or court use

**Do not use this deck for a legal proceeding.** It is not intended for, and is
not recommended for, any calculation that informs a prosecution, a defence, a
paternity ruling, an immigration decision, a disaster-victim identification, or
any other legal or official determination. Jurisdictions require validated,
accredited software for that purpose, operated under an accredited process;
this deck is neither. If a number from here would ever reach a court, reproduce
it with validated software and let that result be the one you rely on.

### No warranty, no liability

This software is provided **as is**, without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose, and non-infringement. **You use it entirely
at your own risk.**

The author accepts no responsibility and no liability for any claim, damage,
loss, or other consequence arising from this software or from any number it
produces, however that number was reached and whoever relied on it. Verifying
that a result is correct, appropriate, and lawful for your purpose is your
responsibility alone.

This restates - and does not narrow, extend, or replace - the warranty and
liability disclaimers of the
[LGPL-3.0-only](https://www.gnu.org/licenses/lgpl-3.0.html) licence this deck
is released under. Where the two differ, the licence governs.

Everything below is a specific, technical limitation. None of it narrows the
two statements above.

## Data

### The example frequency table is synthetic

`examples/lib/exampledb.j` builds a table of made-up numbers. Every figure an
example prints illustrates the arithmetic, not a real population. Use a real
table, for example from [STRidER](https://strider.online/), for anything you
intend to rely on.

### Pooling populations is not an uncertain population

`combineFrequencies` merges tables into one. When the inputs' frequencies
differ, the merged table violates Hardy-Weinberg (the Wahlund effect), and the
match probability computed from it is biased downwards, overstating the
evidence. For "the donor was German or Austrian" use
`randomMatchProbabilityAcross`. In the test suite's example the two differ by a
third at one locus, and the gap compounds across a profile. See
[Reference frequency data](../concepts/frequency-data.md#when-the-population-is-uncertain).

### Reading a published table is a separate deck

The core reads no file format.
[`@deckyard/strider`](https://github.com/deckyard/deck-strider-frequency-converter)
converts STRidER-style XML and has its own limitations, chief among them that
it has not been validated against a real STRidER export.

### No reference sequence is bundled

`lineage.align` takes the rCRS (NC_012920.1) region you sequenced as an
argument.

## Statistical scope

### No theta in kinship or pedigree calculations

[`kinship.j`](../casework/kinship.md) and
[`pedigree.j`](../casework/pedigrees.md) use `θ = 0`. The sub-population
correction does not extend to the IBD-1 term by substituting an allele
frequency; it needs the full recursive sampling formula. Applying θ to some
terms and not others would answer neither question. θ is applied where it is
exact, in the single-source match probability. See
[Theta](../concepts/theta.md).

### Kinship is for outbred pairs

The three IBD coefficients describe non-inbred individuals. An inbred pedigree
needs the nine condensed Jacquard coefficients.

### Peeling handles simple pedigrees only

A consanguineous mating, or any other loop in the individual/family graph, is
rejected by `pedigree.validate`. Loops need cutset conditioning, which is not
implemented. A half-sibling pedigree is not a loop.

### Mixtures stop at CPI / CPE

CPI / CPE needs no contributor count but is weak: it does not use the
suspect's profile, ignores peak heights, and fails under allele dropout. Full
deconvolution (multi-contributor models with drop-in, dropout and peak
heights, as in EuroForMix) is not implemented. See
[Mixtures](../casework/mixtures.md).

### EMPOP's phylogenetic alignment rules are not applied

`lineage.align` returns a plain optimal alignment with indels slid 3'. EMPOP
places indels where the mtDNA phylogeny prefers, which can differ for a
complex indel. Check such a designation against EMPOP before quoting it.

### The default mutation rates are placeholders

`DEFAULT_MUTATION` is `rate = 0.002`, `decay = 0.1`: the right order of
magnitude for paternal autosomal STR mutation, not a published figure. When a
mutation model rescues a trio from a single inconsistency, the CPI is
proportional to the rate at that locus, and observed per-locus rates vary
about tenfold. Set published rates with `kinship.withLocusRate`.

## Performance

### Alignment cost grows with the region, and with divergence

`align` computes a band around the diagonal, so for near-identical sequences
(a haplotype against its own reference) the cost is roughly linear in the
region length. Measured: 342 bp in 0.08 s, the ~1.1 kb control region in
0.29 s, 1600 bp in 0.39 s, all within 50 MB; 14.6 times the length took 15.7
times as long, up to the cap. Divergent sequences widen the band towards the
full quadratic matrix. `align` refuses a region longer than
`lineage.MAX_ALIGN_LENGTH` (5000 bases); at that cap a near-identical pair
takes about 1.2 s.

An eight-locus pedigree LR with four typed people takes about 0.1 s, and
peeling memory does not grow with pedigree size.

### Peeling slows as observed alleles accumulate

Cost grows steeply with the number of distinct alleles the typed members show
at a locus. Measured, one locus: 4 alleles 11 ms, 8 alleles 90 ms, 12 alleles
0.6 s, 16 alleles 1.9 s, 20 alleles 5.4 s. A large family at a highly
polymorphic locus can take seconds per locus. Lumping already limits this to
the observed alleles; beyond that, trim the panel.

### Building a database or profile one entry at a time is quadratic

The copy-returning builders copy the whole container before writing, so a loop
of `withLocus` over `n` loci costs `1 + 2 + ... + n` entry copies. Measured,
10 alleles per locus: 25 loci 5.9 ms, 100 loci 78 ms, 400 loci 1077 ms. A real
STR panel of 15 to 45 loci takes tens of milliseconds, and a 44-marker
reference table converts in about 16-21 ms.

`profile(id, genotypes)` builds in one pass; prefer it to folding
`withGenotype` over a list. `FrequencyDb` has no bulk builder: constructing the
structs directly would skip the validation that rejects a misread table, which
is not worth the saving at this size.

### Iterating a list of large values copies each element

`for (def d in $dbs)` binds a copy of each element, and copying a `FrequencyDb`
copies its locus map; inside a loop over loci that is quadratic. Indexing the
list reads the element in place. The same holds for rotating buffers by
assignment and for `lists.push` / `lists.pop`, which return new lists. The
deck's own code indexes instead (`combineFrequencies`, `align`'s two score
rows, the pedigree graph walks).

### Building a large pedigree is quadratic

`withFounder` and `withChild` return copies, so an `n`-member pedigree copies
the member map `n` times: unnoticeable at tens of members, noticeable at
hundreds. The likelihood functions are unaffected.

### Long panels underflow

`randomMatchProbability` and `pedigreeLikelihood` return plain products. Use
`log10RandomMatchProbability`, `log10PedigreeLikelihood`, or `pedigreeLR`,
which takes the ratio locus by locus.

## Reporting

A number from this deck depends on assumptions it cannot state for you. Report
with it:

- the reference population and the table's source
- θ
- the loci used, especially if a profile was narrowed
- the alternative hypothesis (unrelated, a brother, an unknown man)
- for a probability of paternity, the prior
- for a lineage marker, the database and its size
