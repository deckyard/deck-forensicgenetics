# Limitations

Read this before quoting a number. Some of it will change what you do.

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

### The bundled frequency table is synthetic

`testdata/frequencies.xml` contains **made-up numbers**, generated only to
exercise the reader and the examples. The file says so in its own header
comment. Every figure printed by an example is therefore an illustration of the
arithmetic, not a statement about any real population.

The examples build their table in `examples/lib/exampledb.j`; every number in
it is made up. Download a real table from [STRidER](https://strider.online/)
before computing anything you intend to rely on.

### Pooling populations is not the same as an uncertain population

`combineFrequencies` merges tables into one. That merged table violates
Hardy-Weinberg by construction whenever the inputs' frequencies differ - the
Wahlund effect - so the match probability computed from it is biased
*downwards*, which overstates the evidence.

For "the donor was German or Austrian", use `randomMatchProbabilityAcross`,
which averages the match probability computed within each population rather
than averaging the frequencies. On the worked example in the test suite the two
differ by a third at a single locus, and the gap compounds across a profile.

Reserve `combineFrequencies` for the case where the merged group genuinely is
the reference population.

### Reading a published table is a separate deck

The core reads no published format. [`@mplx/strider`](https://github.com/mplx/deck-strider-frequency-converter)
converts STRidER-style XML, and carries its own limitations - chief among them
that it has not been validated against a real STRidER export.

### No reference sequence is bundled

`lineage.align` takes the rCRS (NC_012920.1) region you sequenced as an
argument rather than shipping a copy. This also keeps the module
database-independent, which is deliberate.

## Statistical scope

### Theta is not applied in kinship or pedigree calculations

[`kinship.j`](../casework/kinship.md) and
[`pedigree.j`](../casework/pedigrees.md) use the plain product rule, `θ = 0`.

The sub-population correction does not extend to the IBD-1 term - the joint
probability of two related people's genotypes under co-ancestry needs the full
recursive sampling formula, not a swapped-in allele frequency. Applying θ to
some terms and not others would produce a number that answers neither question,
so those modules are uniformly theta-free.

θ **is** applied where it is exact: the single-source match probability. See
[Theta](../concepts/theta.md).

### Kinship is for outbred pairs

The three IBD coefficients describe **non-inbred** individuals. An inbred
pedigree needs the nine condensed Jacquard coefficients, which the pairwise
form does not carry.

### Peeling handles simple pedigrees only

A consanguineous mating, or any other loop in the individual/family graph,
makes it cyclic. `pedigree.validate` **rejects** such a pedigree rather than
mis-summing it. Loops need cutset conditioning, which is not implemented.

A half-sibling pedigree - one father, two mothers - is not a loop and is
handled normally.

### Mixture support stops at CPI / CPE

[Mixtures](../casework/mixtures.md) gives the combined probability of inclusion
and exclusion, which needs no assumption about contributor count but is
correspondingly weak: it does not use the suspect's profile, ignores peak
heights, and degrades badly under allele dropout.

Full deconvolution - multi-contributor models with drop-in, dropout and peak
heights, the EuroForMix space - is a larger, later effort layered on this base.
It is not here.

### EMPOP's phylogenetic alignment rules are not applied

`lineage.align` renders a plain optimal alignment with indels slid 3'. EMPOP
shifts indels to the position the mtDNA phylogeny prefers, which can differ for
a complex indel. Check such a designation against EMPOP before quoting it.

### The default mutation rates are placeholders

`DEFAULT_MUTATION` is `rate = 0.002`, `decay = 0.1` - the right order of
magnitude for paternal autosomal STR mutation, not a published figure.

This matters more than it looks. When a mutation model rescues a trio from a
single inconsistency, the resulting CPI is directly proportional to the rate
assumed at that locus, and observed per-locus rates vary over roughly an order
of magnitude. Substitute published rates with
[`withLocusRate`](../reference/kinship.md).

## Performance

### Alignment cost grows with the region, and with divergence

`align` computes only a band around the diagonal, so for near-identical
sequences - a haplotype against its own reference, which is the case it exists
for - the cost is roughly linear in the region length rather than quadratic.
Measured: 342 bp in 0.34 s, the ~1.1 kb control region in 1.5 s, 1600 bp in
2.7 s, all inside 50 MB.

Two sequences that genuinely diverge force the band wider, up to the full
matrix, and then the cost is the quadratic one. `align` refuses a region longer
than `lineage.MAX_ALIGN_LENGTH` (5000 bases) rather than appearing to hang; at
that cap a near-identical pair takes about 15 s.

Everything else in the deck is fast. An eight-locus, four-typed-person pedigree
LR runs in about a tenth of a second, and peeling memory is flat regardless of
pedigree size.

### Peeling slows as observed alleles accumulate

Cost is driven by the number of **distinct alleles the typed members show at a
locus**, growing steeply in the resulting genotype count. Measured, for one
locus: 4 alleles 11 ms, 8 alleles 90 ms, 12 alleles 0.6 s, 16 alleles 1.9 s,
20 alleles 5.4 s. A large reference family at a highly polymorphic locus can
therefore take seconds per locus. Allele lumping already caps this at what is
actually observed; beyond that, the panel is the thing to trim.

### Building a large pedigree is quadratic

`withFounder` and `withChild` return copies, as every builder in this deck does,
so assembling an `n`-member pedigree copies the member map `n` times. That is
invisible for the tens of members a real case has and noticeable in the
hundreds. The likelihood functions themselves are unaffected.

### Long panels underflow

`randomMatchProbability` and `pedigreeLikelihood` return ordinary products and
will underflow on a very long panel. Use `log10RandomMatchProbability`,
`log10PedigreeLikelihood`, or `pedigreeLR`, which takes the ratio locus by
locus and stays in range.

## Reporting

Numbers from this deck are conditional on assumptions the deck cannot state for
you. Whatever you write down, write down alongside it:

- which reference population, and where the table came from
- which θ
- which loci contributed - especially if a profile was narrowed
- what the alternative hypothesis actually is (unrelated? a brother? an unknown
  man?)
- for a probability of paternity, which prior
- for a lineage marker, which database and how large it was

## Related

- [Formulae and sources](formulae.md) - what is implemented, and from where
- [Testing and contributing](testing.md) - how the estimators are checked
