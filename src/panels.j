# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

/**
 * Built-in autosomal STR locus panels - the CODIS and European core sets and
 * the autosomal loci of common commercial kits - with a locus-name normaliser
 * and helpers that narrow a `Profile` or `FrequencyDb` to a panel.
 *
 * Locus names follow the convention the rest of the deck uses: `D3S1358`,
 * `D8S1179` (D-numbered loci in upper case), `TH01`, `vWA`, `FGA`, `CSF1PO`,
 * `TPOX`, `Penta D`, `Penta E`, `SE33`. `normalizeLocus` maps common variants
 * (`THO1`, `VWA`, `PentaE`, `PENTA_E`, `d8s1179`, `HUMTH01`, `ACTBP2`) onto
 * it; `normalizeProfile` and `normalizeDb` rename a whole profile or table so
 * that loci typed under different spellings line up.
 *
 * Each panel lists only the autosomal STR loci. Amelogenin, Y-STRs (`DYS391`,
 * `DYS570`, `DYS576`) and the Y-indel are not autosomal and take no part in
 * an RMP or kinship LR; a kit's such markers are kept in
 * `Panel.otherMarkers` for reference. Loci are listed in one fixed order (the
 * D-numbered loci by chromosome, then the named loci), not in the kit's dye or
 * size order.
 *
 * Sources, one per panel (the `source` field repeats them):
 *
 * - `CODIS13`: Budowle B, Moretti TR, Niezgoda SJ, Brown BL (1998), "CODIS
 *   and PCR-based short tandem repeat loci: law enforcement tools", Proc. 2nd
 *   European Symposium on Human Identification, Promega, pp. 73-88.
 * - `CODIS20`: Hares DR (2015), "Selection and implementation of expanded
 *   CODIS core loci in the United States", Forensic Sci. Int. Genet. 17:33-34;
 *   the FBI required them for NDIS from 1 January 2017.
 * - `ESS`: Council Resolution of 30 November 2009 on the exchange of DNA
 *   analysis results (2009/C 296/01), Official Journal of the EU C 296, Annex 1;
 *   Gill P et al. (2006), Forensic Sci. Int. 156:242-244 for the ENFSI/EDNAP
 *   recommendation of the five added loci.
 * - `GlobalFiler`: Thermo Fisher Scientific, GlobalFiler PCR Amplification Kit
 *   User Guide (21 autosomal STRs, DYS391, a Y-indel, Amelogenin).
 * - `PowerPlex Fusion`: Promega, PowerPlex Fusion System Technical Manual
 *   (22 autosomal STRs, DYS391, Amelogenin).
 * - `PowerPlex Fusion 6C`: Promega, PowerPlex Fusion 6C System Technical
 *   Manual (23 autosomal STRs, DYS391, DYS570, DYS576, Amelogenin).
 * - `PowerPlex ESX 17`: Promega, PowerPlex ESX 17 System Technical Manual
 *   (16 autosomal STRs, Amelogenin).
 * - `NGM SElect`: Thermo Fisher Scientific, AmpFlSTR NGM SElect PCR
 *   Amplification Kit User Guide (16 autosomal STRs, Amelogenin).
 * - `Investigator ESSplex SE QS`: QIAGEN, Investigator ESSplex SE QS Handbook
 *   (16 autosomal STRs, Amelogenin, quality sensors).
 * - `Identifiler Plus`: Thermo Fisher Scientific, AmpFlSTR Identifiler Plus PCR
 *   Amplification Kit User Guide (15 autosomal STRs, Amelogenin).
 * - `Investigator IDplex Plus`: QIAGEN, Investigator IDplex Plus Handbook
 *   (15 autosomal STRs, Amelogenin, quality sensors).
 *
 * Kits are included only where the locus list is certain; others (for example
 * PowerPlex 21, Investigator 24plex) can be built with `customPanel`.
 * @module panels
 * @example
 * import "@deckyard/forensicgenetics/panels.j" as panels;
 *
 * def ess as panels.Panel init panels.panel("ESS");
 * def shared as list of string init panels.sharedLoci(
 *     panels.panel("Identifiler Plus"), $ess);
 * def narrowed as forensics.Profile init panels.restrictProfile($p, $ess);
 * def gaps as list of string init panels.missingInDb($db, $ess);
 */
use maps;
use strings;
use regex;

import "./forensicgenetics.j" as forensics;

/**
 * A named set of autosomal STR loci.
 * @field name {string} the panel's name, e.g. "GlobalFiler"
 * @field source {string} the document the locus list is taken from
 * @field loci {list of string} the autosomal STR loci, normalised names
 * @field otherMarkers {list of string} non-autosomal markers the kit also
 *   types (Amelogenin, Y-STRs, Y-indel); never used in a calculation
 */
export def struct Panel {
    name as string,
    source as string,
    loci as list of string,
    otherMarkers as list of string
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "panels", message: $msg, file: "", line: 0, col: 0};
}

# --- the built-in panels ----------------------------------------------------

# The original 13 CODIS loci; every later US kit is a superset.
def const CODIS13 as list of string init [
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D13S317",
    "D16S539",
    "D18S51",
    "D21S11",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA"
];

# CODIS13 plus the seven loci required from 1 January 2017 (Hares 2015).
def const CODIS20 as list of string init [
    "D1S1656",
    "D2S441",
    "D2S1338",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D13S317",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "D22S1045",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA"
];

# The seven original ESS loci plus the five added in 2009.
def const ESS as list of string init [
    "D1S1656",
    "D2S441",
    "D3S1358",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D18S51",
    "D21S11",
    "D22S1045",
    "FGA",
    "TH01",
    "vWA"
];

# CODIS20 plus SE33.
def const GLOBALFILER as list of string init [
    "D1S1656",
    "D2S441",
    "D2S1338",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D13S317",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "D22S1045",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA",
    "SE33"
];

# CODIS20 plus Penta D and Penta E.
def const FUSION as list of string init [
    "D1S1656",
    "D2S441",
    "D2S1338",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D13S317",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "D22S1045",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA",
    "Penta D",
    "Penta E"
];

# Fusion plus SE33.
def const FUSION6C as list of string init [
    "D1S1656",
    "D2S441",
    "D2S1338",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D13S317",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "D22S1045",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA",
    "Penta D",
    "Penta E",
    "SE33"
];

# ESS plus D2S1338, D16S539, D19S433 and SE33: PowerPlex ESX 17, NGM SElect
# and Investigator ESSplex SE QS type the same sixteen.
def const ESS16 as list of string init [
    "D1S1656",
    "D2S441",
    "D2S1338",
    "D3S1358",
    "D8S1179",
    "D10S1248",
    "D12S391",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "D22S1045",
    "FGA",
    "TH01",
    "vWA",
    "SE33"
];

# CODIS13 plus D2S1338 and D19S433: Identifiler Plus and IDplex Plus.
def const IDENTIFILER as list of string init [
    "D2S1338",
    "D3S1358",
    "D5S818",
    "D7S820",
    "D8S1179",
    "D13S317",
    "D16S539",
    "D18S51",
    "D19S433",
    "D21S11",
    "CSF1PO",
    "FGA",
    "TH01",
    "TPOX",
    "vWA"
];

# catalogue lists every built-in panel, in the order panelNames reports.
func catalogue() {
    return [
        Panel{
            name: "CODIS13",
            source: "Budowle et al. (1998), the original 13 CODIS core loci",
            loci: CODIS13,
            otherMarkers: []
        },
        Panel{
            name: "CODIS20",
            source: "Hares (2015), Forensic Sci. Int. Genet. 17:33-34; required by the FBI from 2017-01-01",
            loci: CODIS20,
            otherMarkers: []
        },
        Panel{
            name: "ESS",
            source: "Council Resolution 2009/C 296/01, Annex 1 (European Standard Set)",
            loci: ESS,
            otherMarkers: []
        },
        Panel{
            name: "GlobalFiler",
            source: "Thermo Fisher Scientific, GlobalFiler PCR Amplification Kit User Guide",
            loci: GLOBALFILER,
            otherMarkers: ["AMEL", "DYS391", "Y indel"]
        },
        Panel{
            name: "PowerPlex Fusion",
            source: "Promega, PowerPlex Fusion System Technical Manual",
            loci: FUSION,
            otherMarkers: ["AMEL", "DYS391"]
        },
        Panel{
            name: "PowerPlex Fusion 6C",
            source: "Promega, PowerPlex Fusion 6C System Technical Manual",
            loci: FUSION6C,
            otherMarkers: ["AMEL", "DYS391", "DYS570", "DYS576"]
        },
        Panel{
            name: "PowerPlex ESX 17",
            source: "Promega, PowerPlex ESX 17 System Technical Manual",
            loci: ESS16,
            otherMarkers: ["AMEL"]
        },
        Panel{
            name: "NGM SElect",
            source: "Thermo Fisher Scientific, AmpFlSTR NGM SElect Kit User Guide",
            loci: ESS16,
            otherMarkers: ["AMEL"]
        },
        Panel{
            name: "Investigator ESSplex SE QS",
            source: "QIAGEN, Investigator ESSplex SE QS Handbook",
            loci: ESS16,
            otherMarkers: ["AMEL"]
        },
        Panel{
            name: "Identifiler Plus",
            source: "Thermo Fisher Scientific, AmpFlSTR Identifiler Plus Kit User Guide",
            loci: IDENTIFILER,
            otherMarkers: ["AMEL"]
        },
        Panel{
            name: "Investigator IDplex Plus",
            source: "QIAGEN, Investigator IDplex Plus Handbook",
            loci: IDENTIFILER,
            otherMarkers: ["AMEL"]
        }
    ];
}

# nameKey folds a panel name for lookup: lower case, with spaces, hyphens and
# underscores removed, so "powerplex-fusion-6c" finds "PowerPlex Fusion 6C".
func nameKey(name as string) {
    def parts as list of string init [];
    for (def ch in strings.chars(strings.lower($name))) {
        if ($ch != " " and $ch != "-" and $ch != "_") {
            $parts[] = $ch;
        }
    }
    return strings.join($parts, "");
}

/**
 * The names of the built-in panels.
 * @return {list of string} every name `panel` accepts, in catalogue order
 */
export func panelNames() {
    def out as list of string init [];
    for (def p in catalogue()) {
        $out[] = $p.name;
    }
    return $out;
}

/**
 * Whether a built-in panel has this name. Case, spaces, hyphens and
 * underscores are ignored, as in `panel`.
 * @param name {string} the panel name
 * @return {bool} true if `panel` would find it
 */
export func hasPanel(name as string) {
    def key as string init nameKey($name);
    for (def p in catalogue()) {
        if (nameKey($p.name) == $key) {
            return true;
        }
    }
    return false;
}

/**
 * A built-in panel by name.
 *
 * The match ignores case, spaces, hyphens and underscores: `"globalfiler"`,
 * `"PowerPlex-ESX-17"` and `"codis20"` all resolve. `panelNames` lists the
 * names.
 * @param name {string} the panel name
 * @return {Panel} the panel
 * @throws {Error} kind "panels" if no built-in panel has that name
 */
export func panel(name as string) {
    def key as string init nameKey($name);
    for (def p in catalogue()) {
        if (nameKey($p.name) == $key) {
            return $p;
        }
    }
    fail("no built-in panel named '" + $name + "'; panelNames() lists them");
}

/**
 * A caller-defined panel, its locus names normalised.
 * @param name {string} a label for the panel
 * @param loci {list of string} the autosomal loci, in any spelling
 *   `normalizeLocus` understands
 * @return {Panel} the panel, with no `otherMarkers` and an empty `source`
 * @throws {Error} kind "panels" on an empty name or locus list, an empty
 *   locus name, or two names that normalise to the same locus
 */
export func customPanel(name as string, loci as list of string) {
    if (len($name) == 0) {
        fail("a panel needs a name");
    }
    if (len($loci) == 0) {
        fail("panel '" + $name + "' lists no loci");
    }
    def seen as map of string to bool init {};
    def out as list of string init [];
    for (def locus in $loci) {
        def canon as string init normalizeLocus($locus);
        if (maps.has($seen, $canon)) {
            fail("panel '" + $name + "' lists locus '" + $canon + "' twice");
        }
        $seen[$canon] = true;
        $out[] = $canon;
    }
    return Panel{name: $name, source: "", loci: $out, otherMarkers: []};
}

/**
 * The loci two panels have in common, in the first panel's order.
 *
 * The usual question when a profile typed with one kit is compared with a
 * database typed with another.
 * @param a {Panel} one panel
 * @param b {Panel} the other
 * @return {list of string} the shared loci; empty if none
 */
export func sharedLoci(a as Panel, b as Panel) {
    def inB as map of string to bool init {};
    for (def locus in $b.loci) {
        $inB[$locus] = true;
    }
    def out as list of string init [];
    for (def locus in $a.loci) {
        if (maps.has($inB, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}

/**
 * The intersection of two panels as a panel of its own, named
 * `"<a> & <b>"`, to pass to `restrictProfile` or `restrictDb`.
 * @param a {Panel} one panel
 * @param b {Panel} the other
 * @return {Panel} a panel of `sharedLoci(a, b)`
 * @throws {Error} kind "panels" if the two share no locus
 */
export func intersect(a as Panel, b as Panel) {
    def loci as list of string init sharedLoci($a, $b);
    if (len($loci) == 0) {
        fail("panels '" + $a.name + "' and '" + $b.name + "' share no locus");
    }
    return Panel{name: $a.name + " & " + $b.name, source: "", loci: $loci, otherMarkers: []};
}

# --- locus names ------------------------------------------------------------

# foldKey upper-cases a locus name and drops the separators people put into it
# (space, underscore, hyphen, dot, slash), so "Penta_E", "penta e" and
# "PENTA-E" all fold to "PENTAE", and "HUMVWFA31/A" to "HUMVWFA31A".
func foldKey(name as string) {
    def key as string init strings.upper(strings.trim($name));
    $key = strings.replace($key, " ", "");
    $key = strings.replace($key, "_", "");
    $key = strings.replace($key, "-", "");
    $key = strings.replace($key, ".", "");
    return strings.replace($key, "/", "");
}

# SEGMENT_NAME is the D-segment form of a locus name: "D", a chromosome
# (digits, X or Y), "S", a number - D8S1179, DXS10135, DYS391.
def const SEGMENT_NAME as string init "^D([0-9]+|X|Y)S[0-9]+$";

# isSegmentName reports whether a folded key has the D-segment form.
func isSegmentName(key as string) {
    return regex.matches(SEGMENT_NAME, $key);
}

# The locus names whose canonical spelling is not upper case, keyed by
# their folded form. HUM-prefixed GenBank names are folded before this lookup.
def const ALIASES as map of string to string init {
    "TH01": "TH01",
    "THO1": "TH01",
    "VWA": "vWA",
    "VWA31": "vWA",
    "VWA31A": "vWA",
    "VWFA31": "vWA",
    "VWFA31A": "vWA",
    "FGA": "FGA",
    "FIBRA": "FGA",
    "CSF1PO": "CSF1PO",
    "CSF1P0": "CSF1PO",
    "TPOX": "TPOX",
    "PENTAD": "Penta D",
    "PENTAE": "Penta E",
    "SE33": "SE33",
    "ACTBP2": "SE33",
    "AMEL": "AMEL",
    "AMELOGENIN": "AMEL"
};

# CANONICAL holds the canonical spellings ALIASES maps to, each its own
# normalised form; normalizeLocus returns these unchanged without folding.
def const CANONICAL as map of string to bool init {
    "TH01": true,
    "vWA": true,
    "FGA": true,
    "CSF1PO": true,
    "TPOX": true,
    "Penta D": true,
    "Penta E": true,
    "SE33": true,
    "AMEL": true
};

/**
 * A locus name in the deck's canonical spelling.
 *
 * Case, spaces, underscores, hyphens, dots and slashes are ignored, and a
 * GenBank `HUM` prefix is dropped before a known name (`HUMTH01`). Known names map to one spelling: `THO1`
 * and `TH01` to `TH01`; `VWA`, `vWA31` and `HUMVWFA31/A` to `vWA`;
 * `FIBRA` to `FGA`; `CSF1P0` to `CSF1PO`; `PentaD`, `PENTA_D` to `Penta D`
 * (and E); `ACTBP2` to `SE33`; `Amelogenin` to `AMEL`. A D-segment name
 * (`d8s1179`, `DYS391`) is upper-cased. Any other name is returned trimmed but
 * otherwise unchanged, so an unfamiliar locus is never mangled.
 * @param name {string} the locus name as written
 * @return {string} the canonical name
 * @throws {Error} kind "panels" on an empty or blank name
 */
export func normalizeLocus(name as string) {
    # Already canonical: a D-segment name in upper case, or a canonical
    # spelling, folds back to itself.
    if (isSegmentName($name) or maps.has(CANONICAL, $name)) {
        return $name;
    }
    def key as string init foldKey($name);
    if (len($key) == 0) {
        fail("locus name must not be empty");
    }
    if (strings.startsWith($key, "HUM") and len($key) > 3) {
        def bare as string init strings.substring($key, 3, len($key));
        if (maps.has(ALIASES, $bare)) {
            return ALIASES[$bare];
        }
    }
    if (maps.has(ALIASES, $key)) {
        return ALIASES[$key];
    }
    if (isSegmentName($key)) {
        return $key;
    }
    return strings.trim($name);
}

/**
 * A copy of `p` with every locus renamed by `normalizeLocus`, in the
 * profile's keys and in each `Genotype.locus`.
 * @param p {Profile} the profile to rename
 * @return {Profile} the renamed copy, loci in the same order
 * @throws {Error} kind "panels" if two of the profile's loci normalise to the
 *   same name (for example both `TH01` and `THO1`)
 */
export func normalizeProfile(p as forensics.Profile) {
    def out as forensics.Profile init forensics.Profile{id: $p.id, loci: {}};
    def origin as map of string to string init {};
    for (def locus in $p.loci) {
        def canon as string init normalizeLocus($locus);
        if (maps.has($origin, $canon)) {
            fail("profile '" + $p.id + "': loci '" + $origin[$canon] + "' and '" + $locus +
                "' both normalise to '" + $canon + "'");
        }
        $origin[$canon] = $locus;
        def g as forensics.Genotype init $p.loci[$locus];
        $out.loci[$canon] = forensics.Genotype{locus: $canon, a: $g.a, b: $g.b};
    }
    return $out;
}

/**
 * A copy of `fdb` with every locus renamed by `normalizeLocus`, in the keys
 * and in each `LocusFrequencies.locus`. Frequencies and sample sizes are
 * unchanged.
 * @param fdb {FrequencyDb} the table to rename
 * @return {FrequencyDb} the renamed copy, loci in the same order
 * @throws {Error} kind "panels" if two of the table's loci normalise to the
 *   same name
 */
export func normalizeDb(fdb as forensics.FrequencyDb) {
    def out as forensics.FrequencyDb init forensics.db($fdb.name);
    for (def locus in $fdb.loci) {
        def canon as string init normalizeLocus($locus);
        if (maps.has($out.loci, $canon)) {
            fail("frequency database '" + $fdb.name + "': two loci normalise to '" + $canon + "'");
        }
        # Written straight into this function's own value: the entry was
        # validated when it entered $fdb, and withLocus would copy the whole
        # table per locus.
        $out.loci[$canon] = forensics.LocusFrequencies{
            locus: $canon,
            n: $fdb.loci[$locus].n,
            alleles: $fdb.loci[$locus].alleles
        };
    }
    return $out;
}

# --- narrowing to a panel ---------------------------------------------------

# canonicalKeys maps each canonical locus name to the key a map actually uses,
# so a profile typed as "THO1" is found when the panel asks for "TH01". Two
# keys with one canonical name throw, as in normalizeProfile / normalizeDb:
# keeping either would hide the other's data. `owner` names the source in the
# message.
func canonicalKeys(keys as list of string, owner as string) {
    def out as map of string to string init {};
    for (def k in $keys) {
        def canon as string init normalizeLocus($k);
        if (maps.has($out, $canon)) {
            fail($owner + ": loci '" + $out[$canon] + "' and '" + $k + "' both normalise to '" +
                $canon + "'");
        }
        $out[$canon] = $k;
    }
    return $out;
}

/**
 * A copy of `p` keeping only the panel's loci, in panel order and under their
 * canonical names.
 *
 * A profile locus is matched to the panel through `normalizeLocus`, so a
 * profile typed as `THO1` contributes its genotype as `TH01`. Panel loci the
 * profile lacks are absent from the result - `missingInProfile` lists
 * them - and loci outside the panel are dropped, which is the point.
 * @param p {Profile} the profile to narrow
 * @param pan {Panel} the panel to keep
 * @return {Profile} the narrowed profile; empty if it shares no locus
 * @throws {Error} kind "panels" if two of the profile's loci normalise to the
 *   same name
 */
export func restrictProfile(p as forensics.Profile, pan as Panel) {
    def norm as forensics.Profile init normalizeProfile($p);
    def out as forensics.Profile init forensics.Profile{id: $p.id, loci: {}};
    for (def locus in $pan.loci) {
        if (maps.has($norm.loci, $locus)) {
            $out.loci[$locus] = $norm.loci[$locus];
        }
    }
    return $out;
}

/**
 * A copy of `fdb` keeping only the panel's loci, in panel order and under
 * their canonical names. Panel loci the table lacks are absent from the
 * result; `missingInDb` lists them.
 * @param fdb {FrequencyDb} the table to narrow
 * @param pan {Panel} the panel to keep
 * @return {FrequencyDb} the narrowed table
 * @throws {Error} kind "panels" if two of the table's loci normalise to the
 *   same name
 */
export func restrictDb(fdb as forensics.FrequencyDb, pan as Panel) {
    def norm as forensics.FrequencyDb init normalizeDb($fdb);
    def out as forensics.FrequencyDb init forensics.db($fdb.name);
    for (def locus in $pan.loci) {
        if (maps.has($norm.loci, $locus)) {
            $out.loci[$locus] = $norm.loci[$locus];
        }
    }
    return $out;
}

/**
 * The panel loci a profile is not typed at, in panel order. Profile loci are
 * matched through `normalizeLocus`.
 * @param p {Profile} the profile to check
 * @param pan {Panel} the panel it should cover
 * @return {list of string} the missing loci; empty if the profile covers the
 *   panel
 * @throws {Error} kind "panels" if two of the profile's loci normalise to the
 *   same name, as `restrictProfile` does
 */
export func missingInProfile(p as forensics.Profile, pan as Panel) {
    def owner as string init "profile '" + $p.id + "'";
    return missingFrom(canonicalKeys(maps.keys($p.loci), $owner), $pan);
}

/**
 * The panel loci a frequency table does not cover, in panel order. Table
 * loci are matched through `normalizeLocus`.
 *
 * Run it before a calculation over a kit's loci: `randomMatchProbability`
 * and the kinship LRs throw on an uncovered locus rather than drop it.
 * @param fdb {FrequencyDb} the table to check
 * @param pan {Panel} the panel it should cover
 * @return {list of string} the missing loci; empty if the table covers the
 *   panel
 * @throws {Error} kind "panels" if two of the table's loci normalise to the
 *   same name, as `restrictDb` does
 */
export func missingInDb(fdb as forensics.FrequencyDb, pan as Panel) {
    def owner as string init "frequency database '" + $fdb.name + "'";
    return missingFrom(canonicalKeys(forensics.loci($fdb), $owner), $pan);
}

# missingFrom lists the panel loci absent from a canonical-name map.
func missingFrom(have as map of string to string, pan as Panel) {
    def out as list of string init [];
    for (def locus in $pan.loci) {
        if (not maps.has($have, $locus)) {
            $out[] = $locus;
        }
    }
    return $out;
}
