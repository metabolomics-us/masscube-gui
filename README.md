# MassCube GUI

A small Tkinter front-end for [MassCube](https://github.com/huaxuyu/masscube)
(Huaxu Yu's untargeted-metabolomics workflow). It lets you pick folders, edit
parameters with sensible defaults already filled in, and run the workflow with
a single button.

## One-time setup

```powershell
pip install -r requirements.txt
```

This installs `masscube` (which pulls in numpy, pandas, scipy, pyteomics,
ms_entropy, etc.). Python 3.9+ recommended.

## Launch

Double-click `run_gui.bat`, or:

```powershell
python masscube_gui.py
```

## Workflow

1. **Project folder** — output directory. The GUI creates `data/`,
   `single_files/`, `chromatograms/`, etc. inside it.
2. **Raw data folder** — folder containing your `.mzML` / `.mzXML` files. If
   it's not already `<project>/data`, the GUI creates a Windows directory
   junction so MassCube finds it (no admin needed).
3. **m/z + RT library (optional)** — CSV with columns `name, mz, rt` (RT in
   minutes). On Run, the GUI copies it to `<project>/mzrt_list.csv`, which
   MassCube auto-detects for targeted annotation. Tolerances reuse
   `mz_tol_alignment` and `rt_tol_alignment` from the Parameters tab.
3. **Parameters tab** — every parameter from MassCube's `Params` class,
   pre-populated with upstream defaults. Tweak as needed. *Reset parameters
   to defaults* restores the originals.
4. **Samples tab** — *Scan raw data folder* lists every detected raw file.
   Double-click a cell to toggle `is_blank` / `is_qc` or set a group label.
5. **Run workflow** — writes `parameters.csv` + `sample_table.csv` into the
   project folder and invokes
   `masscube.workflows.untargeted_metabolomics_workflow`. Live output streams
   to the *Run log* tab.
6. **MSP export** — when the workflow finishes cleanly the GUI converts the
   results into `<project>/aligned_feature_table.msp` (see below).
7. **Export for formula ID** — writes SIRIUS `.ms` input and an MS1-bearing MSP
   for MS-FINDER, for molecular-formula assignment from accurate m/z + isotope
   abundance (see below).

## MSP output

After a successful run the GUI writes an MSP spectral library of the results to
`<project>/aligned_feature_table.msp`, built by reading
`aligned_feature_table.txt` back in. Each entry carries `NAME` (the annotation,
or `Unknown`), `PRECURSORMZ`, `PRECURSORTYPE` (adduct), `RETENTIONTIME` (min),
`IONMODE`, `FORMULA` / `INCHIKEY` / `SMILES` when annotated, a `COMMENT` with
`feature_ID`, `search_mode`, `similarity`, `matched_peaks`, `database`,
`detection_rate` and the MS2 source file, then `Num Peaks` and the peak list.

Two settings on the *Output & Plots* tab control it:

- `output_msp` — run the export after the workflow (default on).
- `msp_content` — `features_with_MS2` (default, only features that actually
  have an MS2 spectrum), `annotated_only` (only annotated features; MS2-less
  ones get `Num Peaks: 0`), or `all_features` (every feature).

Both are GUI-side settings, so they are deliberately **not** written to
`parameters.csv` — that file stays a pure MassCube input.

The **Export MSP** button regenerates the file for whatever project folder is
selected without reprocessing anything, which is handy after tweaking
`msp_content` or for a project that was processed earlier. It falls back to
`project_files/aligned_feature_table_before_annotation.txt` if the annotated
table isn't there yet.

MassCube upstream also writes its own `project_files/features.msp`, but only in
the branch where feature alignment actually runs — a re-run that reuses an
existing `aligned_feature_table.txt` skips it — and it is generated before
normalization. The GUI export is independent of both.

## masscube speed patch — "Step 1: Preparing the workflow..." hang

Stock masscube can sit on `Step 1: Preparing the workflow...` for minutes or
hours with no output. This is upstream behaviour, not a GUI problem: masscube's
workflow preparation calls `get_start_time()` once per raw file (serially) to
sort samples by acquisition time, and that function does

```python
for l in f.readlines()[:300]:
```

`readlines()` materializes the **whole file** before the `[:300]` slice, so
Step 1 reads every byte of every sample just to find one timestamp in the
header. With a cloud-synced data folder (OneDrive, Dropbox, a network share)
it's worse still — the full read forces every placeholder file to download.

`patches/masscube-get-start-time.patch` replaces that with a bounded 256 KB
header read. Return values are identical; on a synthetic 406 MB mzML it went
from 717 ms to 0.5 ms per file.

If Step 1 is still slow afterwards, the bottleneck is your storage: make sure
the raw data folder is on a local disk, or right-click it in Explorer →
*Always keep on this device* so no file is a cloud-only placeholder.

## masscube isotope patch — recovering M+3 / M+4

`find_isotope_signals()` searches for isotope peaks at **integer** Da offsets.
That's deliberate upstream, and correct as far as it goes: integer steps keep the
sub-1 Da isotopes inside the tolerance window (¹⁵N +0.99703, ³⁴S +1.99580,
³⁷Cl +1.99705). The gap it leaves is pure carbon — ¹³C is +1.0033548 per step,
so the deviation grows 0.00336 each step and crosses the default
`mz_tol_feature_grouping` of 0.01 at M+3:

| isotope | offset | deviation from integer grid | |
|---|---|---|---|
| ¹³C ×1 | +1.00335 | 0.00335 | inside |
| ¹³C ×2 | +2.00671 | 0.00671 | inside |
| ¹³C ×3 | +3.01006 | 0.01006 | **outside — silently lost** |
| ¹³C ×4 | +4.01342 | 0.01342 | **outside — silently lost** |

So the `isotopes` column normally carries only M+0…M+2. Naively switching to
1.00336 spacing fixes M+3/M+4 but *breaks* ³⁴S, ³⁴S+¹³C and ³⁷Cl×2, so
`patches/masscube-isotope-spacing.patch` searches **both** grids and keeps a peak
matching either. That's a strict superset of upstream, so nothing that matched
before stops matching. Measured over the common isotopes at tol=0.01:

| grid | caught | regressions |
|---|---|---|
| integer (upstream) | 10/12 | — |
| ¹³C-spaced (naive) | 9/12 | ³⁴S, ³⁴S+¹³C, ³⁷Cl×2 |
| union (this patch) | **12/12** | none |

This one only affects features detected during a run, so it needs a **fresh
workflow run** to change any output.

## Applying the patches

```powershell
python scripts/apply_masscube_patch.py            # apply all
python scripts/apply_masscube_patch.py --check     # report status only
python scripts/apply_masscube_patch.py --revert    # restore from backups
python scripts/apply_masscube_patch.py --only isotope-spacing
```

Idempotent, and each patched file is backed up alongside the original as
`*.py.orig`. Because these edit the installed package in `site-packages`,
**re-run after every `pip install --upgrade masscube`** — an upgrade overwrites
the files and silently reverts both fixes. If masscube changes upstream such that
the expected code is no longer found, the script reports that and refuses to
guess rather than corrupting the file.

## Molecular formula assignment (SIRIUS / MS-FINDER)

Assigning a formula from accurate m/z **and isotope abundance** needs the MS1
isotope pattern, not MS2. MassCube keeps it — the feature table's `isotopes`
column holds the apex-scan signals, monoisotopic peak included — but stores it as
a numpy array repr spanning two physical lines per row, which is awkward to
consume directly.

*Export for formula ID* (button, or `scripts/export_for_formula_id.py`) writes:

- **`sirius_input/`** — one `.ms` per feature with `>parentmass`, `>charge`,
  `>ionization`, `>ms1peaks` and `>ms2peaks`. Import the folder into SIRIUS, or
  run `sirius -i sirius_input -o sirius_out formulas`. Use `--single-file` for one
  multi-compound `sirius_input.ms` instead.
- **`aligned_feature_table_msfinder.msp`** — MSP with `MSTYPE: MS1`
  (isotope pattern) and `MSTYPE: MS2` sections, for MS-FINDER.

Both deliberately **omit the formula** even when MassCube annotated the feature —
writing it would pin the answer instead of letting the tool derive it. The
annotation is carried as a `#` comment in the `.ms`, which SIRIUS ignores.

Rows that aren't compounds are filtered out, and the run log reports every count:

- `is_isotope=True` rows — these *are* the M+1/M+2 peaks of other features
- `is_in_source_fragment` rows — fragments, not intact molecules
  (`--keep-in-source-fragments` to keep them)
- multimer adducts (`[2M+H]+`, `[3M-H]-`) — the mass describes a cluster

Charge comes from the adduct string, because MassCube marks `charge` as
non-exported in its feature-table schema; with no adduct it falls back to
`ion_mode`. Adducts in `SIRIUS_IONIZATION_MAP` are passed through exactly;
anything else (formate, acetate, methanol/acetonitrile adducts, uncommon metals)
becomes the SIRIUS wildcard `[M+?]+/-` so SIRIUS considers adducts itself rather
than rejecting a string it may not parse. Add entries to that map if your SIRIUS
build accepts more, and note the plain **Export MSP** output is unaffected by any
of this.

For a one-off manual check of a single feature, [ChemCalc](https://www.chemcalc.org/)
or enviPat is quicker than either tool.

## Notes

- Defaults in `masscube_gui.py` mirror `masscube/params.py` from upstream
  `main`. If you upgrade `masscube` and a default changes, update
  `PARAM_GROUPS` at the top of `masscube_gui.py` accordingly.
- `ms2_library_path` is optional — leave it blank to skip MS2 annotation.
- **MS2 library check / fix.** Next to the MS2 library picker on the
  *Annotation* tab there's a *Check / Fix* button. Use it to scan the file
  for entries MassCube would crash on (it requires every entry to have a
  field whose name contains both `prec` and `mz`). If your `.msp` uses
  variants like `PrecursorMass`, `Precursor m/z`, `ParentMass`, or `PEPMASS`,
  the GUI offers to rename them to `PRECURSORMZ` in place. A `.msp.bak`
  backup is written next to the original on first fix. The same check also
  runs automatically before *Run workflow*, so a bad library aborts cleanly
  instead of crashing mid-way through MassCube.
- *Save parameters.csv only* writes the CSV without running the workflow, in
  case you want to inspect or edit it by hand first.
