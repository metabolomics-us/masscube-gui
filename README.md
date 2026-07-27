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
from 717 ms to 0.5 ms per file. Apply it with:

```powershell
python scripts/apply_masscube_patch.py
```

The script is idempotent, backs up the original to `utils_functions.py.orig`,
and supports `--check` and `--revert`. Because it edits the installed package in
`site-packages`, **re-run it after every `pip install --upgrade masscube`**.

If Step 1 is still slow afterwards, the bottleneck is your storage: make sure
the raw data folder is on a local disk, or right-click it in Explorer →
*Always keep on this device* so no file is a cloud-only placeholder.

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
