# Profiles and genotypes

## Alleles are strings

An STR allele designation counts repeats of a short motif, and most of them are
whole numbers - but not all. `"9.3"` at TH01 is nine full repeats plus a
three-base partial one; it is common, it is not `9`, and it is not `9.5`.
`"15.2"` at D19S433 is the same idea. So the deck stores allele labels as
**strings**, verbatim:

```jennifer
def g as forensics.Genotype init forensics.genotype("TH01", "6", "9.3");
```

If they were numbers, `"15"` and `"15.0"` would silently become one key,
`"9.3"` would invite floating-point comparison, and off-ladder designations
like `"OL"` or amelogenin's `"X"` and `"Y"` would have nowhere to live.

Where repeat-count arithmetic is genuinely wanted - the stepwise mutation model
is the only place - ask for it explicitly:

```jennifer
forensics.alleleValue("9.3");        # 9.3
forensics.isNumericAllele("9.3");    # true
forensics.isNumericAllele("X");      # false
forensics.alleleValue("X");          # throws
```

## Genotype

A genotype is an unordered pair at one locus. `genotype` puts the pair in
canonical order, so the two ways of writing it are the same value and compare
equal:

```jennifer
forensics.genotype("TH01", "9.3", "6") == forensics.genotype("TH01", "6", "9.3")
```

Ordering is by repeat count when both alleles have one, and lexicographic
otherwise - so a panel that mixes `"15"` with `"X"` still has a total order,
and `"9.3"` sorts after `"6"` rather than before it as a string comparison
would have it.

```jennifer
forensics.isHomozygous(forensics.genotype("TH01", "6", "6"));   # true
forensics.alleles(forensics.genotype("TH01", "6", "6"));        # ["6"]
forensics.alleles(forensics.genotype("TH01", "6", "7"));        # ["6", "7"]
```

`alleles` returns the *distinct* alleles: one label for a homozygote, two for a
heterozygote. That is the set that matters when asking whether a mixture covers
someone.

## Profile

A profile is one genotype per locus, keyed by locus name:

```jennifer
def p as forensics.Profile init forensics.profile("suspect", [
    forensics.genotype("D3S1358", "15", "16"),
    forensics.genotype("vWA", "17", "17"),
    forensics.genotype("TH01", "6", "9.3")
]);
```

Listing a locus twice is an error, not a last-wins overwrite. The map is
insertion-ordered, so `profileLoci` returns the loci in the order they were
given.

```jennifer
forensics.profileLoci($p);              # ["D3S1358", "vWA", "TH01"]
forensics.hasGenotype($p, "FGA");       # false
forensics.genotypeAt($p, "TH01");       # the Genotype; throws if untyped
```

To add or replace a genotype, `withGenotype` returns a copy:

```jennifer
def q as forensics.Profile init forensics.withGenotype($p,
    forensics.genotype("FGA", "22", "24"));
```

### Value semantics

Jennifer copies on assignment and on argument passing, and this deck leans on
that: every builder is copy-returning and nothing mutates its argument.

```jennifer
def q as forensics.Profile init forensics.withGenotype($p, $g);
# $p is unchanged
```

That means a profile can be handed to several calculations without any of them
being able to disturb it, and a `FrequencyDb` can be shared freely.

### Narrowing a profile

A profile and a frequency table need not cover the same panel. The deck will
not guess:

```jennifer
def shared as list of string init forensics.commonLoci($db, $p);
def narrowed as forensics.Profile init forensics.restrict($p, $shared);
```

`commonLoci` is the intersection in profile order; `restrict` keeps exactly the
named loci and throws if the profile is untyped at one of them. Doing this
deliberately is the point - see
[Match probability](../casework/match-probability.md) for why
`randomMatchProbability` refuses to do it for you.

## Mixture

A mixed stain is not a profile. There is no genotype to record, only the set of
alleles seen at each locus, with no assumption about how many people
contributed them or how they pair up:

```jennifer
def stain as forensics.Mixture init forensics.mixture("stain", {
    "D3S1358": ["15", "16", "17"],
    "vWA": ["16", "17"],
    "TH01": ["6", "7", "9.3"]
});
```

Repeated alleles in a locus's list are collapsed; an empty list is an error.
`isIncluded` asks the non-exclusion question directly:

```jennifer
forensics.isIncluded($suspect, $stain);   # bool
```

It skips loci the mixture does not cover and loci where the profile is untyped,
and returns false as soon as one locus shows the profile carries an allele the
stain does not. See [Mixtures](../casework/mixtures.md) for what that is worth
numerically.
