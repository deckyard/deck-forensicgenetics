# `panels`

```jennifer
import "@deckyard/forensicgenetics/panels.j" as panels;
```

Built-in autosomal STR locus panels and a locus-name normaliser. It takes the
core's `Profile` and `FrequencyDb`. Signatures and the `Panel` fields are in
the [Cheatsheet](cheatsheet.md#panels).

## Panels

| Name | Loci | Autosomal loci | Other markers | Source |
| --- | --- | --- | --- | --- |
| `CODIS13` | 13 | CSF1PO, FGA, TH01, TPOX, vWA, D3S1358, D5S818, D7S820, D8S1179, D13S317, D16S539, D18S51, D21S11 | - | Budowle et al. (1998) |
| `CODIS20` | 20 | `CODIS13` + D1S1656, D2S441, D2S1338, D10S1248, D12S391, D19S433, D22S1045 | - | Hares (2015) |
| `ESS` | 12 | D3S1358, vWA, D8S1179, D21S11, D18S51, TH01, FGA + D1S1656, D2S441, D10S1248, D12S391, D22S1045 | - | Council Resolution 2009/C 296/01 |
| `GlobalFiler` | 21 | `CODIS20` + SE33 | AMEL, DYS391, Y indel | Thermo Fisher user guide |
| `PowerPlex Fusion` | 22 | `CODIS20` + Penta D, Penta E | AMEL, DYS391 | Promega technical manual |
| `PowerPlex Fusion 6C` | 23 | `PowerPlex Fusion` + SE33 | AMEL, DYS391, DYS570, DYS576 | Promega technical manual |
| `PowerPlex ESX 17` | 16 | `ESS` + D2S1338, D16S539, D19S433, SE33 | AMEL | Promega technical manual |
| `NGM SElect` | 16 | as `PowerPlex ESX 17` | AMEL | Thermo Fisher user guide |
| `Investigator ESSplex SE QS` | 16 | as `PowerPlex ESX 17` | AMEL | QIAGEN handbook |
| `Identifiler Plus` | 15 | `CODIS13` + D2S1338, D19S433 | AMEL | Thermo Fisher user guide |
| `Investigator IDplex Plus` | 15 | as `Identifiler Plus` | AMEL | QIAGEN handbook |

Full citations: the module docblock and
[Formulae and sources](../project/formulae.md#kit-panels).

## Conventions

- Panels hold autosomal STRs only. Amelogenin, Y-STRs and the Y-indel are in
  `Panel.otherMarkers`, never in `loci`.
- Loci are in one fixed order: D-numbered loci by chromosome, then CSF1PO,
  FGA, TH01, TPOX, vWA, Penta D, Penta E, SE33. Not the kit's dye order.
- `panel` and `hasPanel` ignore case, spaces, hyphens and underscores in the
  name.
- `normalizeLocus` folds case, drops spaces, `_`, `-`, `.` and `/`, and drops
  a `HUM` prefix when the rest is a known alias. Aliases: `THO1` -> `TH01`, `VWA` / `VWFA31A`
  -> `vWA`, `FIBRA` -> `FGA`, `CSF1P0` -> `CSF1PO`, `PENTAD` / `PENTAE` ->
  `Penta D` / `Penta E`, `ACTBP2` -> `SE33`, `Amelogenin` -> `AMEL`.
  D-segment names (`d8s1179`, `dys391`) are upper-cased; anything else is
  returned trimmed and otherwise unchanged.
- `restrictProfile` and `restrictDb` match through `normalizeLocus` and return
  canonical names in panel order. Panel loci absent from the input are absent
  from the result; `missingInProfile` / `missingInDb` list them.

## Errors

`Error{kind: "panels"}`, thrown for:

- an unknown panel name
- an empty or blank locus name
- a custom panel with no name, no loci, or a locus listed twice (after
  normalisation)
- two panels with no shared locus (`intersect`)
- two loci of one profile or table that normalise to the same name (also in
  `missingInProfile` / `missingInDb`)
