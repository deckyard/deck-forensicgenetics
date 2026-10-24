# Profiles and genotypes

## Alleles are strings

Allele labels are stored verbatim as strings. `"9.3"` at TH01 and `"15.2"` at
D19S433 are ordinary designations; `"15"` and `"15.0"` are different keys; and
`"OL"`, `"X"` and `"Y"` need no special case.

```jennifer
def g as forensicgenetics.Genotype init forensicgenetics.genotype("TH01", "6", "9.3");
```

Where repeat-count arithmetic is meant (only the stepwise mutation model needs
it), convert explicitly:

```jennifer
forensicgenetics.alleleValue("9.3");        # 9.3
forensicgenetics.isNumericAllele("9.3");    # true
forensicgenetics.isNumericAllele("X");      # false
forensicgenetics.alleleValue("X");          # throws
```

## Genotype

`genotype` puts the pair in canonical order, so both writing orders give equal
values:

```jennifer
forensicgenetics.genotype("TH01", "9.3", "6") == forensicgenetics.genotype("TH01", "6", "9.3")
```

Two numeric alleles order by repeat count (`"6"` before `"9.3"`); otherwise
the order is lexicographic.

```jennifer
forensicgenetics.isHomozygous(forensicgenetics.genotype("TH01", "6", "6"));   # true
forensicgenetics.alleles(forensicgenetics.genotype("TH01", "6", "6"));        # ["6"]
forensicgenetics.alleles(forensicgenetics.genotype("TH01", "6", "7"));        # ["6", "7"]
```

`alleles` returns the distinct alleles.

## Profile

One genotype per locus:

```jennifer
def p as forensicgenetics.Profile init forensicgenetics.profile("suspect", [
    forensicgenetics.genotype("D3S1358", "15", "16"),
    forensicgenetics.genotype("vWA", "17", "17"),
    forensicgenetics.genotype("TH01", "6", "9.3")
]);

forensicgenetics.profileLoci($p);              # ["D3S1358", "vWA", "TH01"]
forensicgenetics.hasGenotype($p, "FGA");       # false
forensicgenetics.genotypeAt($p, "TH01");       # the Genotype; throws if untyped
```

A locus listed twice throws. Loci keep insertion order.

`withGenotype` returns a copy with one genotype added or replaced:

```jennifer
def q as forensicgenetics.Profile init forensicgenetics.withGenotype($p,
    forensicgenetics.genotype("FGA", "22", "24"));
# $p is unchanged
```

Every builder in the deck returns a copy and leaves its argument alone, so a
profile or a `FrequencyDb` can be passed to several calculations safely.

### Narrowing a profile

The deck does not drop loci a table lacks; narrow the profile explicitly:

```jennifer
def shared as list of string init forensicgenetics.commonLoci($db, $p);
def narrowed as forensicgenetics.Profile init forensicgenetics.restrict($p, $shared);
```

`commonLoci` returns the intersection in profile order. `restrict` keeps the
named loci and throws if the profile is untyped at one. See
[Match probability](../casework/match-probability.md#missing-loci-are-refused-not-skipped).

## Mixture

A mixed stain records the alleles seen per locus, with no assumption about the
number of contributors:

```jennifer
def stain as forensicgenetics.Mixture init forensicgenetics.mixture("stain", {
    "D3S1358": ["15", "16", "17"],
    "vWA": ["16", "17"],
    "TH01": ["6", "7", "9.3"]
});

forensicgenetics.isIncluded($suspect, $stain);   # bool
```

Repeated alleles collapse; an empty list throws. For `isIncluded` and the
statistics see [Mixtures](../casework/mixtures.md).
