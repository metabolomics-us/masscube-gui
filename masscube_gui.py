"""
MassCube GUI — a tkinter front-end for Huaxu Yu's MassCube untargeted
metabolomics workflow (https://github.com/huaxuyu/masscube).

Pick a project folder, pick a folder of raw .mzML / .mzXML files, tweak the
parameters (pre-populated with MassCube's own defaults), then hit Run. The GUI
writes parameters.csv + sample_table.csv into the project folder and invokes
masscube.workflows.untargeted_metabolomics_workflow.
"""

import csv
import json
import os
import pickle
import queue
import re
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import (
    BooleanVar,
    StringVar,
    Tk,
    filedialog,
    messagebox,
    ttk,
    scrolledtext,
    END,
    N, S, E, W,
)


# ---------------------------------------------------------------------------
# Defaults pulled verbatim from masscube/params.py (Params.__init__).
# Keep these in sync with upstream if the package is updated.
# Each entry: (default_value, kind, choices_or_None, help_text)
#   kind ∈ {"float", "int", "bool", "choice", "str", "path"}
# ---------------------------------------------------------------------------

# Which features end up in the exported MSP (see write_msp_from_feature_table).
MSP_CONTENT_CHOICES = ["features_with_MS2", "annotated_only", "all_features"]

PARAM_GROUPS = {
    "Raw data": [
        ("ion_mode",          ("positive", "choice", ["positive", "negative"],
                               "MS ionization mode")),
        ("ms_type",           ("orbitrap", "choice", ["orbitrap", "qtof", "tripletof", "others"],
                               "Instrument type. Affects intensity defaults.")),
        ("is_centroid",       (True, "bool", None,
                               "Raw data is already centroided")),
        ("scan_time_unit",    ("minute", "choice", ["minute", "second"],
                               "Time unit of scan time")),
        ("mz_lower_limit",    (0.0, "float", None, "Lower m/z limit (Da)")),
        ("mz_upper_limit",    (100000.0, "float", None, "Upper m/z limit (Da)")),
        ("rt_lower_limit",    (0.0, "float", None, "Lower RT limit (min)")),
        ("rt_upper_limit",    (10000.0, "float", None, "Upper RT limit (min)")),
        ("centroid_mz_tol",   (0.002, "float", None,
                               "m/z tolerance for centroiding")),
        ("ms1_abs_int_tol",   (1000.0, "float", None,
                               "Absolute MS1 intensity threshold (Orbitrap≈30000, QTOF≈1000)")),
        ("ms2_abs_int_tol",   (500.0, "float", None,
                               "Absolute MS2 intensity threshold (Orbitrap≈10000, QTOF≈500)")),
        ("ms2_rel_int_tol",   (0.01, "float", None,
                               "Relative MS2 intensity threshold vs base peak")),
        ("precursor_mz_offset", (2.0, "float", None,
                                 "Offset for MS2 m/z range (Da)")),
    ],
    "Feature detection": [
        ("mz_tol_ms1",        (0.01, "float", None, "MS1 m/z tolerance")),
        ("mz_tol_ms2",        (0.015, "float", None, "MS2 m/z tolerance")),
        ("feature_gap_tol",   (10, "int", None,
                               "Consecutive scans without signal allowed within a feature")),
        ("batch_size",        (100, "int", None, "Parallel processing batch size")),
        ("percent_cpu_to_use", (0.8, "float", None, "Fraction of CPU cores to use")),
    ],
    "Grouping": [
        ("group_features_single_file", (False, "bool", None,
                                        "Group features within each single file")),
        ("scan_scan_cor_tol", (0.9, "float", None, "Scan-to-scan correlation tolerance")),
        ("mz_tol_feature_grouping", (0.01, "float", None, "m/z tolerance for grouping")),
        ("rt_tol_feature_grouping", (0.05, "float", None, "RT tolerance for grouping")),
        ("isotope_rel_int_limit", (1.5, "float", None,
                                   "Isotope intensity upper limit (×base peak)")),
    ],
    "Alignment": [
        ("mz_tol_alignment",   (0.01, "float", None, "m/z tolerance for alignment")),
        ("rt_tol_alignment",   (0.2, "float", None, "RT tolerance for alignment")),
        ("noise_tol",          (2.0, "float", None, "Noise score tolerance")),
        ("gaussian_similarity_tol", (0.7, "float", None, "Gaussian similarity tolerance")),
        ("rt_tol_rt_correction", (0.5, "float", None, "Expected max RT shift (min)")),
        ("correct_rt",         (True, "bool", None, "Perform RT correction")),
        ("scan_number_cutoff", (5, "int", None, "Min non-zero scans to include feature")),
        ("detection_rate_cutoff", (0.1, "float", None,
                                   "Detection rate × (qc+sample) cutoff")),
        ("merge_features",     (True, "bool", None, "Merge features with same m/z & RT")),
        ("mz_tol_merge_features", (0.005, "float", None, "m/z tol when merging")),
        ("rt_tol_merge_features", (0.03, "float", None, "RT tol when merging")),
        ("group_features_after_alignment", (True, "bool", None,
                                            "Group features after alignment")),
        ("fill_gaps",          (True, "bool", None, "Fill gaps in aligned features")),
        ("gap_filling_method", ("local_maximum", "choice", ["local_maximum"],
                                "Gap-filling strategy")),
        ("gap_filling_rt_window", (0.05, "float", None, "RT window for local maximum (min)")),
    ],
    "Annotation": [
        ("ms2_library_path",   ("", "path", None,
                                "Path to MS2 library (.msp or .pickle). Leave blank to skip.")),
        ("fuzzy_search",       (False, "bool", None, "Fuzzy search for annotation")),
        ("consider_rt",        (False, "bool", None, "Consider RT in MS2 matching")),
        ("rt_tol_annotation",  (0.2, "float", None, "RT tolerance for MS2 annotation")),
        ("ms2_sim_tol",        (0.7, "float", None, "MS2 similarity tolerance")),
        ("spectral_similarity_method",
            ("unweighted_entropy", "choice",
             ["unweighted_entropy", "entropy", "dot_product"],
             "Spectral similarity method")),
    ],
    "Normalization & Stats": [
        ("sample_normalization", (False, "bool", None,
                                  "Normalize by total sample amount")),
        ("sample_norm_method",   ("pqn", "choice", ["pqn"],
                                  "Sample normalization method. MassCube currently "
                                  "implements only PQN; total_intensity / median / "
                                  "quantile / mdfc are stubs upstream.")),
        ("signal_normalization", (False, "bool", None,
                                  "Correct systematic signal drift")),
        ("signal_norm_method",   ("lowess", "choice", ["lowess"],
                                  "Signal-drift normalization method")),
        ("run_statistics",       (False, "bool", None, "Run statistical analysis")),
    ],
    "Output & Plots": [
        ("plot_bpc",             (True, "bool", None, "Plot base peak chromatograms")),
        ("plot_ms2",             (False, "bool", None, "Plot MS2 mirror plots")),
        ("plot_normalization",   (False, "bool", None, "Plot normalization results")),
        ("output_single_file",   (True, "bool", None,
                                  "Export per-sample processed CSVs")),
        ("output_ms1_scans",     (True, "bool", None,
                                  "Export MS1 scans as pickle for fast reload")),
        ("output_aligned_file",  (True, "bool", None, "Export aligned-feature CSV")),
        ("quant_method",         ("peak_height", "choice",
                                   ["peak_height", "peak_area", "top_average"],
                                   "Quantification value to export")),
        ("output_msp",           (True, "bool", None,
                                  "GUI extra: after the run, export "
                                  "aligned_feature_table.msp from the results")),
        ("msp_content",          ("features_with_MS2", "choice", MSP_CONTENT_CHOICES,
                                  "Which features go into the MSP: only those with an "
                                  "MS2 spectrum, only annotated ones, or every feature")),
    ],
}

# Flat dict {name: spec} for quick lookup.
PARAM_SPECS = {n: spec for grp in PARAM_GROUPS.values() for n, spec in grp}

# Settings that live in the GUI only — they drive GUI-side post-processing and
# are deliberately kept out of parameters.csv, which is MassCube's own input.
GUI_ONLY_PARAMS = {"output_msp", "msp_content"}

RAW_EXTS = (".mzml", ".mzxml", ".mzjson", ".mzjson.gz")


# ---------------------------------------------------------------------------
# MS2 library validator + auto-fixer.
#
# Two failure modes are detected here:
#
# 1. Precursor m/z header. MassCube's annotation.py picks the precursor m/z
#    by searching each entry's keys for any name containing BOTH the
#    substrings "prec" and "mz" (case-insensitive). Common offenders:
#    "PrecursorMass", "Precursor m/z", "PEPMASS" (MGF style).
#
# 2. Ion_mode value. MassCube builds an ion-mode mask in annotate_aligned_
#    features (annotation.py:202+) that requires entry['ion_mode'].lower()
#    to equal params.ion_mode ("positive" or "negative"). Anything else
#    ("P", "N", "+", "POS", missing field, …) silently makes that entry
#    score 0 against every feature — so the workflow runs to completion but
#    yields zero annotations.
# ---------------------------------------------------------------------------

# Headers that almost certainly mean "precursor m/z" but won't pass MassCube's
# substring test. Used both for validation reporting and for the auto-fix.
_FIXABLE_KEY_RE = re.compile(
    r"^(precursor\s*mass|precursor\s*m\s*/?\s*z|parent\s*mass|"
    r"parent\s*m\s*/?\s*z|pepmass)$",
    re.IGNORECASE,
)
# Same pattern, but anchored to a line with a key:value (or key=value) form.
_FIXABLE_LINE_RE = re.compile(
    r"^(\s*)(precursor\s*mass|precursor\s*m\s*/?\s*z|parent\s*mass|"
    r"parent\s*m\s*/?\s*z|pepmass)(\s*[:=])",
    re.IGNORECASE | re.MULTILINE,
)


def _key_is_precursor_mz(k: str) -> bool:
    kl = k.lower()
    return "prec" in kl and "mz" in kl


def _key_looks_fixable(k: str) -> bool:
    return bool(_FIXABLE_KEY_RE.match(k.strip())) and not _key_is_precursor_mz(k)


# Ion_mode line in an .msp entry. Captures prefix / key / sep / value / trailing
# whitespace so the fixer can rewrite the value without disturbing layout.
_ION_MODE_LINE_RE = re.compile(
    r"^(\s*)(ion[_\s]?mode|ionmode)(\s*[:=]\s*)(.+?)(\s*)$",
    re.IGNORECASE | re.MULTILINE,
)
# Map of recognized-but-non-canonical Ion_mode values → MassCube-acceptable form.
# Keys are matched after .strip().lower(); any value not in this map AND not
# already "positive" / "negative" is reported as unfixable.
_ION_MODE_NORMALIZE = {
    "p": "Positive", "pos": "Positive", "positive": "Positive", "+": "Positive",
    "n": "Negative", "neg": "Negative", "negative": "Negative", "-": "Negative",
}


def _ion_mode_acceptable(value: str) -> bool:
    """True iff MassCube's mask test (value.lower() == params.ion_mode) will pass."""
    return value.strip().lower() in ("positive", "negative")


# A decimal-comma numeric token ('197,9867', '1,2e3') flanked by anything that
# isn't a digit or dot — so US-style "1,234.56" thousands separators (digit/dot
# adjacent) are left alone. Used both to detect locale-broken values and to
# rewrite them.
_COMMA_NUM_RE = re.compile(
    r"(?<![\d.])(-?\d+),(\d+(?:[eE][+-]?\d+)?)(?![\d.])"
)


def _has_decimal_comma(value: str) -> bool:
    return bool(_COMMA_NUM_RE.search(value))


def _ion_mode_fix(value: str) -> str | None:
    """Map a raw Ion_mode value to MassCube's expected casing, or None if unknown."""
    return _ION_MODE_NORMALIZE.get(value.strip().lower())


def validate_ms2_library(path: Path) -> dict:
    """
    Check a library file the way MassCube will. Returns:
        {
          "total":   int,  # entries parsed
          "bad":     list[(idx, name)],     # entries missing precursor m/z
          "fixable": list[str],             # offending header names that can be renamed
          "bad_ion_mode":          list[(idx, name, raw_value)],  # value won't match
          "bad_ion_mode_unfixable": list[str], # raw values we can't auto-map
          "missing_ion_mode":      list[(idx, name)],  # no Ion_mode line at all
          "format":  "msp" | "pickle" | "json",
          "error":   str (only on hard failure)
        }
    """
    ext = path.suffix.lower()
    if ext == ".msp":
        return _validate_msp(path)
    if ext in (".pkl", ".pickle"):
        return _validate_pickled(path)
    if ext == ".json":
        return _validate_json(path)
    return {"error": f"Unsupported library extension: {ext}"}


def _validate_msp(path: Path) -> dict:
    bad: list[tuple[int, str]] = []
    fixable: set[str] = set()
    bad_ion_mode: list[tuple[int, str, str]] = []
    missing_ion_mode: list[tuple[int, str]] = []
    bad_decimal: list[tuple[int, str, str]] = []  # precursor m/z values with ','
    bad_unfixable_count = 0   # bad entries with no renameable header either
    bad_ms1_count = 0         # subset of bad_unfixable that is Spectrum_type: MS1
    decimal_comma_peak_lines = 0
    total = 0
    cur_keys: list[str] = []
    cur_name = ""
    cur_ion_mode: str | None = None
    cur_spectrum_type: str | None = None
    cur_precursor_value: str | None = None

    def _close_entry() -> None:
        nonlocal total, bad_unfixable_count, bad_ms1_count
        if not cur_keys:
            return
        total += 1
        label = cur_name or "<no name>"
        has_pmz = any(_key_is_precursor_mz(k) for k in cur_keys)
        has_renameable = any(_key_looks_fixable(k) for k in cur_keys)
        if not has_pmz:
            bad.append((total, label))
            if not has_renameable:
                bad_unfixable_count += 1
                if cur_spectrum_type and cur_spectrum_type.strip().lower() == "ms1":
                    bad_ms1_count += 1
        for k in cur_keys:
            if _key_looks_fixable(k):
                fixable.add(k)
        if cur_ion_mode is None:
            missing_ion_mode.append((total, label))
        elif not _ion_mode_acceptable(cur_ion_mode):
            bad_ion_mode.append((total, label, cur_ion_mode))
        if cur_precursor_value is not None and _has_decimal_comma(cur_precursor_value):
            bad_decimal.append((total, label, cur_precursor_value))

    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    _close_entry()
                    cur_keys, cur_name, cur_ion_mode = [], "", None
                    cur_spectrum_type = None
                    cur_precursor_value = None
                    continue
                if ":" in stripped and not stripped[0].isdigit():
                    k, _, v = stripped.partition(":")
                    k = k.strip()
                    if k:
                        cur_keys.append(k)
                        kl = k.lower().replace(" ", "").replace("_", "")
                        if kl in ("name", "compoundname"):
                            cur_name = v.strip()
                        elif kl == "ionmode":
                            cur_ion_mode = v.strip()
                        elif kl == "spectrumtype":
                            cur_spectrum_type = v.strip()
                        if cur_precursor_value is None and (
                                _key_is_precursor_mz(k) or _key_looks_fixable(k)):
                            cur_precursor_value = v.strip()
                elif stripped[0].isdigit():
                    # Peak-list line; check for comma-decimal m/z or intensity.
                    if _has_decimal_comma(stripped):
                        decimal_comma_peak_lines += 1
            _close_entry()
    except OSError as exc:
        return {"error": f"Could not read file: {exc}"}

    bad_ion_mode_unfixable = sorted({
        v for _, _, v in bad_ion_mode if _ion_mode_fix(v) is None
    })
    return {
        "total": total,
        "bad": bad,
        "fixable": sorted(fixable),
        "bad_unfixable_count": bad_unfixable_count,
        "bad_ms1_count": bad_ms1_count,
        "bad_ion_mode": bad_ion_mode,
        "bad_ion_mode_unfixable": bad_ion_mode_unfixable,
        "missing_ion_mode": missing_ion_mode,
        "bad_decimal": bad_decimal,
        "decimal_comma_peak_lines": decimal_comma_peak_lines,
        "format": "msp",
    }


def _validate_pickled(path: Path) -> dict:
    try:
        with path.open("rb") as f:
            db = pickle.load(f)
    except Exception as exc:
        return {"error": f"Could not load pickle: {exc}"}
    return _validate_dict_list(db, "pickle")


def _validate_json(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as f:
            db = json.load(f)
    except Exception as exc:
        return {"error": f"Could not load JSON: {exc}"}
    return _validate_dict_list(db, "json")


def _validate_dict_list(db, fmt: str) -> dict:
    if not isinstance(db, list):
        return {"error": f"Library is not a list of dict entries "
                          f"(got {type(db).__name__})"}
    bad: list[tuple[int, str]] = []
    fixable: set[str] = set()
    bad_ion_mode: list[tuple[int, str, str]] = []
    missing_ion_mode: list[tuple[int, str]] = []
    bad_unfixable_count = 0
    bad_ms1_count = 0
    for i, entry in enumerate(db, start=1):
        if not isinstance(entry, dict):
            bad.append((i, "<non-dict entry>"))
            bad_unfixable_count += 1
            continue
        keys = [str(k) for k in entry.keys()]
        label = str(entry.get("name", entry.get("compound_name", "<no name>")))
        has_pmz = any(_key_is_precursor_mz(k) for k in keys)
        has_renameable = any(_key_looks_fixable(k) for k in keys)
        if not has_pmz:
            bad.append((i, label))
            if not has_renameable:
                bad_unfixable_count += 1
                spectrum_type = None
                for k in keys:
                    if k.lower().replace(" ", "").replace("_", "") == "spectrumtype":
                        spectrum_type = str(entry[k])
                        break
                if spectrum_type and spectrum_type.strip().lower() == "ms1":
                    bad_ms1_count += 1
        for k in keys:
            if _key_looks_fixable(k):
                fixable.add(k)
        # look for any ion-mode-like key (lowercased, stripped of _ and space)
        im_val = None
        for k in keys:
            if k.lower().replace(" ", "").replace("_", "") == "ionmode":
                im_val = entry[k]
                break
        if im_val is None:
            missing_ion_mode.append((i, label))
        elif not _ion_mode_acceptable(str(im_val)):
            bad_ion_mode.append((i, label, str(im_val)))
    bad_ion_mode_unfixable = sorted({
        v for _, _, v in bad_ion_mode if _ion_mode_fix(v) is None
    })
    return {
        "total": len(db),
        "bad": bad,
        "fixable": sorted(fixable),
        "bad_unfixable_count": bad_unfixable_count,
        "bad_ms1_count": bad_ms1_count,
        "bad_ion_mode": bad_ion_mode,
        "bad_ion_mode_unfixable": bad_ion_mode_unfixable,
        "missing_ion_mode": missing_ion_mode,
        "format": fmt,
    }


def _rewrite_msp_streaming(path: Path, line_sub_fn) -> int:
    """
    Stream `path` line-by-line through `line_sub_fn`, writing to a sibling
    .tmp file. `line_sub_fn(body: str) -> tuple[str, int]` receives the line
    without its trailing newline and returns (new_body, num_subs). On finish:
    if any substitutions were made, a .msp.bak backup of the original is
    written (if absent) and the temp file atomically replaces the original;
    otherwise the temp file is removed. Memory use is O(line), so multi-GB
    libraries work.
    """
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    backup = path.with_suffix(path.suffix + ".bak")
    changed = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as src, \
             tmp_path.open("w", encoding="utf-8", newline="") as dst:
            for line in src:
                if line.endswith("\r\n"):
                    body, eol = line[:-2], "\r\n"
                elif line.endswith("\n"):
                    body, eol = line[:-1], "\n"
                else:
                    body, eol = line, ""
                new_body, n = line_sub_fn(body)
                changed += n
                dst.write(new_body)
                dst.write(eol)
    except BaseException:
        try:
            tmp_path.unlink()
        except OSError:
            pass
        raise

    if changed:
        if not backup.exists():
            shutil.copyfile(path, backup)
        os.replace(tmp_path, path)
    else:
        try:
            tmp_path.unlink()
        except OSError:
            pass
    return changed


def fix_msp_headers(path: Path) -> int:
    """
    Rename precursor-m/z-like headers in an .msp file to PRECURSORMZ.
    Writes a sibling .msp.bak backup of the original. Returns # of lines changed.
    """
    def _sub_line(body: str) -> tuple[str, int]:
        n = [0]

        def _sub(m: re.Match) -> str:
            replacement = f"{m.group(1)}PRECURSORMZ{m.group(3)}"
            if replacement != m.group(0):
                n[0] += 1
            return replacement

        return _FIXABLE_LINE_RE.sub(_sub, body), n[0]

    return _rewrite_msp_streaming(path, _sub_line)


def fix_msp_ion_mode(path: Path) -> int:
    """
    Normalize Ion_mode / IONMODE values in an .msp file so MassCube's mask
    test (value.lower() == 'positive' | 'negative') succeeds. Recognized
    shorthand ('P', 'N', '+', '-', 'POS', 'NEG', etc.) is rewritten to
    'Positive' / 'Negative'. Unknown values are left alone. Writes a
    sibling .msp.bak backup. Returns # of value lines changed.
    """
    def _sub_line(body: str) -> tuple[str, int]:
        n = [0]

        def _sub(m: re.Match) -> str:
            prefix, key, sep, value, trailing = m.groups()
            # Leave already-acceptable values alone (don't re-case "positive" → "Positive").
            if _ion_mode_acceptable(value):
                return m.group(0)
            fixed = _ion_mode_fix(value)
            if fixed is None:
                return m.group(0)
            n[0] += 1
            return f"{prefix}{key}{sep}{fixed}{trailing}"

        return _ION_MODE_LINE_RE.sub(_sub, body), n[0]

    return _rewrite_msp_streaming(path, _sub_line)


def fix_msp_decimal_commas(path: Path) -> int:
    """
    Rewrite decimal-comma numbers ('197,9867' → '197.9867') in:
      - header lines whose key is a precursor m/z field (recognized or renameable)
      - peak-list lines (any line starting with a digit)
    Synonyms, SMILES, InChI, Comment, and other free-text lines are left alone
    so that legitimate commas in text are preserved. Writes a sibling .msp.bak.
    Returns # of substitutions made.
    """
    def _sub_line(body: str) -> tuple[str, int]:
        stripped = body.strip()
        if not stripped:
            return body, 0
        if stripped[0].isdigit():
            return _COMMA_NUM_RE.subn(r"\1.\2", body)
        if ":" in stripped:
            k = stripped.partition(":")[0].strip()
            if k and (_key_is_precursor_mz(k) or _key_looks_fixable(k)):
                return _COMMA_NUM_RE.subn(r"\1.\2", body)
        return body, 0

    return _rewrite_msp_streaming(path, _sub_line)


def strip_msp_bad_entries(path: Path) -> int:
    """
    Remove .msp entries that have no precursor-m/z field MassCube can recognize.
    An entry is kept if any of its keys either already matches MassCube's test
    (contains both 'prec' and 'mz') OR is a known renameable variant
    (Precursor m/z, ParentMass, PEPMASS, ...) — so this is safe to run before
    or after fix_msp_headers. Typical victims are Spectrum_type: MS1 / in-source
    CID entries from MoNA-style dumps. Writes a sibling .msp.bak backup.
    Returns # of entries removed.
    """
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    backup = path.with_suffix(path.suffix + ".bak")
    removed = 0
    entry_lines: list[str] = []
    entry_has_pmz = False
    entry_has_content = False

    try:
        with path.open("r", encoding="utf-8", errors="replace") as src, \
             tmp_path.open("w", encoding="utf-8", newline="") as dst:
            for line in src:
                entry_lines.append(line)
                stripped = line.strip()
                if not stripped:
                    if entry_has_content and not entry_has_pmz:
                        removed += 1
                    else:
                        for el in entry_lines:
                            dst.write(el)
                    entry_lines = []
                    entry_has_pmz = False
                    entry_has_content = False
                    continue
                entry_has_content = True
                if not entry_has_pmz and ":" in stripped and not stripped[0].isdigit():
                    k = stripped.partition(":")[0].strip()
                    if k and (_key_is_precursor_mz(k) or _key_looks_fixable(k)):
                        entry_has_pmz = True
            # Trailing entry without a closing blank line.
            if entry_lines:
                if entry_has_content and not entry_has_pmz:
                    removed += 1
                else:
                    for el in entry_lines:
                        dst.write(el)
    except BaseException:
        try:
            tmp_path.unlink()
        except OSError:
            pass
        raise

    if removed:
        if not backup.exists():
            shutil.copyfile(path, backup)
        os.replace(tmp_path, path)
    else:
        try:
            tmp_path.unlink()
        except OSError:
            pass
    return removed


# ---------------------------------------------------------------------------
# MSP export of the workflow results.
#
# MassCube writes project_files/features.msp on its own, but only in the branch
# where feature alignment actually runs — a re-run that reuses an existing
# aligned_feature_table.txt skips it — and it lands in project_files/ rather
# than next to the results. This exporter reads the aligned feature table back
# and writes <project>/aligned_feature_table.msp, so it works either way and
# can also be re-run on an old project without reprocessing anything.
#
# MS2 spectra are stored in the table's "MS2" column as
#   "mz1;intensity1|mz2;intensity2|…"
# (masscube.utils_functions.convert_signals_to_string).
# ---------------------------------------------------------------------------
MSP_FILE_NAME = "aligned_feature_table.msp"

# Feature-table candidates, best first. The normalized table carries the same
# annotation/MS2 columns and only differs in the per-sample intensities, which
# an MSP doesn't use, so the aligned table is always the right source.
FEATURE_TABLE_CANDIDATES = (
    ("aligned_feature_table.txt",),
    ("project_files", "aligned_feature_table_before_annotation.txt"),
)

# Values pandas/MassCube can leave in a text column to mean "nothing here".
_EMPTY_CELL_VALUES = {"", "nan", "none", "null"}

# Extra fields folded into the COMMENT line: (column, label).
_MSP_COMMENT_FIELDS = (
    ("feature_ID", "feature_ID"),
    ("search_mode", "search_mode"),
    ("similarity", "similarity"),
    ("matched_peak_number", "matched_peaks"),
    ("database", "database"),
    ("detection_rate", "detection_rate"),
    ("MS2_reference_file", "MS2_file"),
)

# Comment fields where 0 means "no match", not a real measurement — MassCube
# fills these with zeros for unannotated features, which is just noise here.
_MSP_COMMENT_SKIP_IF_ZERO = {"similarity", "matched_peaks"}


def _cell(row: dict, col: str) -> str:
    """Read a feature-table cell, normalizing MassCube's empties to ''."""
    value = row.get(col)
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in _EMPTY_CELL_VALUES else text


def _fmt_num(x: float, digits: int = 4) -> str:
    s = f"{x:.{digits}f}".rstrip("0").rstrip(".")
    return s or "0"


def _comment_value(row: dict, col: str, label: str) -> str | None:
    """Format one COMMENT field, or None if it carries no information."""
    value = _cell(row, col)
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return value
    if number == 0.0 and label in _MSP_COMMENT_SKIP_IF_ZERO:
        return None
    # pandas turns integer columns into floats on round-trip: 3.0 -> 3
    return _fmt_num(number, 6)


def _parse_ms2_string(raw: str) -> list[tuple[float, float]]:
    """Parse "mz;int|mz;int|…" into [(mz, intensity), …], skipping junk."""
    peaks: list[tuple[float, float]] = []
    for chunk in raw.split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(";")
        if len(parts) < 2:
            continue
        try:
            peaks.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return peaks


def find_feature_table(project: Path) -> Path | None:
    for parts in FEATURE_TABLE_CANDIDATES:
        candidate = project.joinpath(*parts)
        if candidate.is_file():
            return candidate
    return None


def write_msp_from_feature_table(
    table_path: Path,
    msp_path: Path,
    ion_mode: str = "",
    content: str = "features_with_MS2",
    isotope_ms1: bool = False,
) -> dict:
    """
    Convert a MassCube aligned feature table into an MSP spectral library.

    Parameters
    ----------
    table_path : Path
        aligned_feature_table.txt (tab-separated, as MassCube writes it).
    msp_path : Path
        Destination .msp file (overwritten).
    ion_mode : str
        "positive" / "negative" — written as IONMODE. Skipped when blank.
    content : str
        One of MSP_CONTENT_CHOICES:
          features_with_MS2 — only features that carry an MS2 spectrum (default)
          annotated_only    — only features with an annotation
          all_features      — every row; MS2-less features get "Num Peaks: 0"
    isotope_ms1 : bool
        Write the MS1 isotope pattern alongside MS2 using MS-FINDER's
        "MSTYPE: MS1" / "MSTYPE: MS2" sections, so a formula search can use
        isotope abundance. Off by default — a plain MSP stays exactly as it was.
        Also drops isotope-peak rows and in-source fragments, which are not
        compounds and would pollute a formula search.

    Returns
    -------
    dict with keys: rows, written, with_peaks, with_isotopes, isotope_peaks,
    skipped_no_ms2, skipped_unannotated, skipped_no_mz, skipped_isotope_row,
    skipped_in_source.
    """

    if content not in MSP_CONTENT_CHOICES:
        raise ValueError(f"Unknown MSP content option: {content!r}")

    stats = {"rows": 0, "written": 0, "with_peaks": 0,
             "with_isotopes": 0, "isotope_peaks": 0,
             "skipped_no_ms2": 0, "skipped_unannotated": 0, "skipped_no_mz": 0,
             "skipped_isotope_row": 0, "skipped_in_source": 0,
             "skipped_multimer": 0}

    ion_mode_label = ion_mode.strip().capitalize() if ion_mode.strip() else ""

    with table_path.open("r", newline="", encoding="utf-8-sig") as src:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [c for c in ("m/z", "RT", "MS2") if c not in fields]
        if isotope_ms1 and "isotopes" not in fields:
            missing.append("isotopes")
        if missing:
            raise ValueError(
                f"{table_path.name} is missing expected column(s): "
                f"{', '.join(missing)}. Is this a MassCube aligned feature table?")

        tmp_path = msp_path.with_suffix(msp_path.suffix + ".tmp")
        try:
            with tmp_path.open("w", encoding="utf-8", newline="\n") as out:
                for row in reader:
                    stats["rows"] += 1

                    mz = _cell(row, "m/z")
                    if not mz:
                        stats["skipped_no_mz"] += 1
                        continue

                    isotopes: list[tuple[float, float]] = []
                    if isotope_ms1:
                        # These rows are not compounds — see the module comment.
                        if _cell_is_true(row, "is_isotope"):
                            stats["skipped_isotope_row"] += 1
                            continue
                        if _cell_is_true(row, "is_in_source_fragment"):
                            stats["skipped_in_source"] += 1
                            continue
                        # [2M+H]+ measures a dimer, so a formula search on it
                        # describes the cluster rather than the compound.
                        adduct_cell = _cell(row, "adduct")
                        if adduct_cell and _MULTIMER_RE.match(adduct_cell):
                            stats["skipped_multimer"] += 1
                            continue
                        isotopes = _parse_isotope_cell(_cell(row, "isotopes"))

                    ms2_raw = _cell(row, "MS2")
                    peaks = _parse_ms2_string(ms2_raw) if ms2_raw else []
                    annotation = _cell(row, "annotation")

                    if content == "features_with_MS2" and not peaks:
                        stats["skipped_no_ms2"] += 1
                        continue
                    if content == "annotated_only" and not annotation:
                        stats["skipped_unannotated"] += 1
                        continue

                    # NAME must stay single-line for every MSP reader.
                    name = " ".join(annotation.split()) if annotation else "Unknown"

                    out.write(f"NAME: {name}\n")
                    out.write(f"PRECURSORMZ: {mz}\n")
                    adduct = _cell(row, "adduct")
                    if adduct:
                        out.write(f"PRECURSORTYPE: {adduct}\n")
                    rt = _cell(row, "RT")
                    if rt:
                        out.write(f"RETENTIONTIME: {rt}\n")
                    if ion_mode_label:
                        out.write(f"IONMODE: {ion_mode_label}\n")
                    for col, key in (("formula", "FORMULA"),
                                     ("InChIKey", "INCHIKEY"),
                                     ("SMILES", "SMILES")):
                        value = _cell(row, col)
                        if value:
                            out.write(f"{key}: {value}\n")

                    comment_parts = []
                    for col, label in _MSP_COMMENT_FIELDS:
                        value = _comment_value(row, col, label)
                        if value is not None:
                            comment_parts.append(f"{label}={value}")
                    comment = "; ".join(comment_parts)
                    if comment:
                        out.write(f"COMMENT: {comment}\n")

                    if isotope_ms1:
                        # MS-FINDER's MSP dialect: explicit MS1/MS2 sections, so
                        # the isotope pattern is read as MS1 rather than as
                        # fragments. MS1 comes first, matching MS-FINDER's own
                        # exports.
                        out.write("MSTYPE: MS1\n")
                        out.write(f"Num Peaks: {len(isotopes)}\n")
                        for peak_mz, intensity in isotopes:
                            out.write(f"{_fmt_num(peak_mz)}\t{_fmt_num(intensity)}\n")
                        out.write("MSTYPE: MS2\n")
                        if isotopes:
                            stats["with_isotopes"] += 1
                            stats["isotope_peaks"] += len(isotopes)

                    out.write(f"Num Peaks: {len(peaks)}\n")
                    for peak_mz, intensity in peaks:
                        out.write(f"{_fmt_num(peak_mz)}\t{_fmt_num(intensity)}\n")
                    out.write("\n")

                    stats["written"] += 1
                    if peaks:
                        stats["with_peaks"] += 1
        except BaseException:
            try:
                tmp_path.unlink()
            except OSError:
                pass
            raise

    os.replace(tmp_path, msp_path)
    return stats


# ---------------------------------------------------------------------------
# Formula-ID export (SIRIUS .ms, and MS1-bearing MSP for MS-FINDER).
#
# Assigning a molecular formula from accurate m/z + isotope abundance needs the
# MS1 isotope pattern, not the MS2 spectrum. MassCube keeps it: the feature
# table's "isotopes" column holds the apex-scan signals as a numpy array repr,
#   "[[300.1234 1000000.] \n [301.1268 220000.]]"
# and it *does* include the monoisotopic peak (find_isotope_signals keeps any
# signal below base_intensity * isotope_rel_int_limit, and that limit defaults
# to 1.5, i.e. above the base peak).
#
# Two caveats that shape the code below:
#
#  * Only M+0..M+2 are present. MassCube searches for isotopes at integer Da
#    offsets (mz + arange(n)) while real 13C spacing is 1.00336 Da, so by M+3
#    the drift (0.0101) exceeds mz_tol_feature_grouping (0.01) and the peak is
#    missed. patches/masscube-isotope-spacing.patch fixes this upstream; the
#    exporter works either way and reports how many isotope peaks it saw.
#
#  * Rows with is_isotope=True *are* the M+1/M+2 peaks of other features, and
#    in-source fragments are not intact molecules. Both are excluded by default
#    — submitting them yields formulas for things that aren't compounds.
# ---------------------------------------------------------------------------
SIRIUS_DIR_NAME = "sirius_input"
SIRIUS_SINGLE_FILE_NAME = "sirius_input.ms"
MSFINDER_MSP_FILE_NAME = "aligned_feature_table_msfinder.msp"

# MassCube adduct -> SIRIUS ionization string, for adducts SIRIUS is known to
# parse. Anything else (formate/acetate/methanol/acetonitrile adducts, metals
# beyond the common ones) is exported as the SIRIUS wildcard "[M+?]+/-", which
# tells SIRIUS to consider adducts itself rather than trusting a string it may
# reject. Add entries here if your SIRIUS build accepts more.
SIRIUS_IONIZATION_MAP = {
    "[M+H]+": "[M+H]+",
    "[M+Na]+": "[M+Na]+",
    "[M+K]+": "[M+K]+",
    "[M+NH4]+": "[M+NH4]+",
    "[M]+": "[M]+",
    "[M+Li]+": "[M+Li]+",
    "[M+2H]2+": "[M+2H]2+",
    "[M+3H]3+": "[M+3H]3+",
    "[M+H-H2O]+": "[M-H2O+H]+",
    "[M-H]-": "[M-H]-",
    "[M+Cl]-": "[M+Cl]-",
    "[M-2H]2-": "[M-2H]2-",
    "[M-3H]3-": "[M-3H]3-",
    "[M-H-H2O]-": "[M-H2O-H]-",
}

# [2M+H]+, [3M-H]- and friends: the measured mass is a dimer/trimer, so a
# formula search on it describes the cluster, not the compound. Skipped.
_MULTIMER_RE = re.compile(r"^\[\s*([2-9]\d*)\s*M")

# Trailing charge of an adduct string: "[M+H]+" -> 1, "[M-2H]2-" -> -2.
_ADDUCT_CHARGE_RE = re.compile(r"\]\s*(\d*)\s*([+-])\s*$")

# Any float in a numpy array repr, including "3.0012340e+02" scientific form.
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

_TRUE_CELL_VALUES = {"true", "1", "yes", "t"}


def _cell_is_true(row: dict, col: str) -> bool:
    """Read a boolean feature-table cell. Missing/blank counts as False."""
    return _cell(row, col).lower() in _TRUE_CELL_VALUES


def _parse_isotope_cell(raw: str) -> list[tuple[float, float]]:
    """
    Parse MassCube's "isotopes" cell into [(m/z, intensity), …], m/z ascending.

    The cell is str() of an (n, 2) numpy array, so it spans two physical lines
    per row and may use scientific notation. numpy elides very large arrays with
    "...", which would silently drop peaks — refuse those rather than guess.
    """
    if not raw or "..." in raw:
        return []
    tokens = _FLOAT_RE.findall(raw)
    if not tokens or len(tokens) % 2:
        return []
    peaks: list[tuple[float, float]] = []
    for i in range(0, len(tokens), 2):
        try:
            peaks.append((float(tokens[i]), float(tokens[i + 1])))
        except ValueError:
            return []
    return sorted(peaks)


def _adduct_charge(adduct: str, ion_mode: str) -> int:
    """
    Charge of an adduct string, falling back to ion mode.

    MassCube drops the `charge` column from the exported table (its export
    schema marks it export=False), so it has to come from the adduct.
    """
    match = _ADDUCT_CHARGE_RE.search(adduct) if adduct else None
    if match:
        magnitude = int(match.group(1)) if match.group(1) else 1
        return magnitude if match.group(2) == "+" else -magnitude
    return -1 if ion_mode.strip().lower().startswith("neg") else 1


def _sirius_ionization(adduct: str, charge: int) -> tuple[str, bool]:
    """Map an adduct to a SIRIUS ionization string. Returns (value, is_exact)."""
    mapped = SIRIUS_IONIZATION_MAP.get(adduct)
    if mapped:
        return mapped, True
    return ("[M+?]+" if charge > 0 else "[M+?]-"), False


def _safe_stem(text: str) -> str:
    """Filesystem- and SIRIUS-safe compound id."""
    cleaned = re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("_")
    return cleaned[:80] or "feature"


def write_sirius_ms_from_feature_table(
    table_path: Path,
    out_path: Path,
    ion_mode: str = "",
    single_file: bool = False,
    include_ms2: bool = True,
    require_isotopes: bool = False,
    skip_in_source_fragments: bool = True,
) -> dict:
    """
    Convert a MassCube aligned feature table into SIRIUS .ms input.

    Parameters
    ----------
    table_path : Path
        aligned_feature_table.txt (tab-separated, as MassCube writes it).
    out_path : Path
        Directory to fill with one .ms per feature, or the .ms file to write
        when single_file is True. Either way it is replaced, not merged.
    ion_mode : str
        "positive" / "negative", used only when a row has no adduct.
    single_file : bool
        Write every compound into one multi-compound .ms instead of a directory.
        A directory of per-compound files is the most portable SIRIUS input.
    include_ms2 : bool
        Emit >ms2peaks when the feature has an MS2 spectrum. MS2 sharply
        improves SIRIUS's formula ranking, so this is on by default.
    require_isotopes : bool
        Skip features whose isotope pattern has fewer than two peaks. Off by
        default: SIRIUS can still rank formulas from accurate mass + MS2.
    skip_in_source_fragments : bool
        Skip rows flagged is_in_source_fragment.

    Notes
    -----
    No >formula line is written even when MassCube annotated the feature —
    that would pin the answer instead of letting SIRIUS derive it. The
    annotation is emitted as a "#" comment, which SIRIUS ignores.

    Returns
    -------
    dict of counters (see summarize_formula_export_stats).
    """

    stats = {"rows": 0, "written": 0, "with_isotopes": 0, "with_ms2": 0,
             "isotope_peaks": 0, "wildcard_adduct": 0, "skipped_isotope_row": 0,
             "skipped_in_source": 0, "skipped_multimer": 0,
             "skipped_no_isotopes": 0, "skipped_no_mz": 0}

    with table_path.open("r", newline="", encoding="utf-8-sig") as src:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [c for c in ("m/z", "RT", "isotopes") if c not in fields]
        if missing:
            raise ValueError(
                f"{table_path.name} is missing expected column(s): "
                f"{', '.join(missing)}. Is this a MassCube aligned feature table?")

        entries: list[tuple[str, str]] = []  # (compound id, .ms body)

        for row in reader:
            stats["rows"] += 1

            mz_text = _cell(row, "m/z")
            try:
                mz = float(mz_text)
            except ValueError:
                stats["skipped_no_mz"] += 1
                continue

            # M+1/M+2 rows describe other features' isotopes, not compounds.
            if _cell_is_true(row, "is_isotope"):
                stats["skipped_isotope_row"] += 1
                continue
            if skip_in_source_fragments and _cell_is_true(row, "is_in_source_fragment"):
                stats["skipped_in_source"] += 1
                continue

            adduct = _cell(row, "adduct")
            if adduct and _MULTIMER_RE.match(adduct):
                stats["skipped_multimer"] += 1
                continue

            isotopes = _parse_isotope_cell(_cell(row, "isotopes"))
            if require_isotopes and len(isotopes) < 2:
                stats["skipped_no_isotopes"] += 1
                continue

            ms2 = _parse_ms2_string(_cell(row, "MS2")) if include_ms2 else []
            charge = _adduct_charge(adduct, ion_mode)
            ionization, exact = _sirius_ionization(adduct, charge)
            if not exact:
                stats["wildcard_adduct"] += 1

            feature_id = _cell(row, "feature_ID") or str(stats["rows"])
            # pandas round-trips integer columns as floats: "12.0" -> "12".
            # (_fmt_num is no help here — with digits=0 it would strip the
            # trailing zero of a real value and turn 10 into 1.)
            if feature_id.endswith(".0"):
                feature_id = feature_id[:-2]
            annotation = " ".join(_cell(row, "annotation").split())
            compound = _safe_stem(f"FT{feature_id}")

            lines = [f">compound {compound}",
                     f">parentmass {_fmt_num(mz, 5)}",
                     f">charge {charge}",
                     f">ionization {ionization}"]

            rt_text = _cell(row, "RT")
            try:
                # MassCube reports RT in minutes; SIRIUS .ms expects seconds.
                lines.append(f">rt {_fmt_num(float(rt_text) * 60.0, 2)}s")
            except ValueError:
                pass

            if annotation:
                lines.append(f"#annotation {annotation}")
            if adduct:
                lines.append(f"#masscube_adduct {adduct}")
            if not exact and adduct:
                lines.append("#note adduct not mapped to SIRIUS, using [M+?]")

            if isotopes:
                lines.append("")
                lines.append(">ms1peaks")
                lines.extend(f"{_fmt_num(p, 5)} {_fmt_num(i, 2)}" for p, i in isotopes)
                stats["with_isotopes"] += 1
                stats["isotope_peaks"] += len(isotopes)
            if ms2:
                lines.append("")
                lines.append(">ms2peaks")
                lines.extend(f"{_fmt_num(p, 5)} {_fmt_num(i, 2)}" for p, i in ms2)
                stats["with_ms2"] += 1

            entries.append((compound, "\n".join(lines) + "\n"))
            stats["written"] += 1

    _write_sirius_entries(entries, out_path, single_file)
    return stats


def _write_sirius_entries(entries: list[tuple[str, str]], out_path: Path,
                          single_file: bool) -> None:
    """Write .ms entries either as one file or one file per compound."""
    if single_file:
        tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
        try:
            with tmp_path.open("w", encoding="utf-8", newline="\n") as out:
                for _, body in entries:
                    out.write(body)
                    out.write("\n")
        except BaseException:
            tmp_path.unlink(missing_ok=True)
            raise
        os.replace(tmp_path, out_path)
        return

    out_path.mkdir(parents=True, exist_ok=True)
    # Clear previous .ms files so a re-export can't leave stale compounds behind
    # for SIRIUS to pick up. Only our own extension, never the whole directory.
    for stale in out_path.glob("*.ms"):
        try:
            stale.unlink()
        except OSError:
            pass
    used: set[str] = set()
    for compound, body in entries:
        name = compound
        suffix = 2
        while name.lower() in used:
            name = f"{compound}_{suffix}"
            suffix += 1
        used.add(name.lower())
        (out_path / f"{name}.ms").write_text(body, encoding="utf-8", newline="\n")


def summarize_formula_export_stats(stats: dict, out_path: Path,
                                   table_path: Path) -> str:
    written = stats["written"]
    lines = [
        f"Wrote {written} compound{'' if written == 1 else 's'} to:",
        f"  {out_path}",
        f"Source: {table_path.name} ({stats['rows']} rows)",
        f"With MS1 isotope pattern: {stats['with_isotopes']}"
        f" ({stats['isotope_peaks']} isotope peaks total)",
        f"With MS2 spectrum: {stats['with_ms2']}",
    ]
    if stats["wildcard_adduct"]:
        lines.append(f"Used [M+?] wildcard ionization for "
                     f"{stats['wildcard_adduct']} feature(s) with an adduct "
                     f"SIRIUS may not parse.")
    for key, label in (("skipped_isotope_row", "isotope peak row(s) of other features"),
                       ("skipped_in_source", "in-source fragment(s)"),
                       ("skipped_multimer", "multimer adduct(s) ([2M+H]+ etc.)"),
                       ("skipped_no_isotopes", "feature(s) without an isotope pattern"),
                       ("skipped_no_mz", "row(s) without an m/z")):
        if stats.get(key):
            lines.append(f"Skipped {stats[key]} {label}.")
    if written and not stats["with_isotopes"]:
        lines.append("WARNING: no isotope patterns found. Was the workflow run "
                     "with feature grouping enabled?")
    return "\n".join(lines)


def summarize_msp_stats(stats: dict, msp_path: Path, table_path: Path) -> str:
    lines = [
        f"Wrote {stats['written']} entr{'y' if stats['written'] == 1 else 'ies'} "
        f"({stats['with_peaks']} with MS2 peaks) to:",
        f"  {msp_path}",
        f"Source: {table_path.name} ({stats['rows']} features)",
    ]
    if stats.get("with_isotopes"):
        lines.append(f"With MS1 isotope pattern: {stats['with_isotopes']}"
                     f" ({stats['isotope_peaks']} isotope peaks total)")
    if stats["skipped_no_ms2"]:
        lines.append(f"Skipped {stats['skipped_no_ms2']} feature(s) without MS2.")
    if stats["skipped_unannotated"]:
        lines.append(f"Skipped {stats['skipped_unannotated']} unannotated feature(s).")
    if stats["skipped_no_mz"]:
        lines.append(f"Skipped {stats['skipped_no_mz']} row(s) without an m/z.")
    if stats.get("skipped_isotope_row"):
        lines.append(f"Skipped {stats['skipped_isotope_row']} isotope peak row(s) "
                     "of other features.")
    if stats.get("skipped_in_source"):
        lines.append(f"Skipped {stats['skipped_in_source']} in-source fragment(s).")
    if stats.get("skipped_multimer"):
        lines.append(f"Skipped {stats['skipped_multimer']} multimer adduct(s) "
                     "([2M+H]+ etc.).")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class MassCubeGUI:
    def __init__(self, root: Tk):
        self.root = root
        root.title("MassCube GUI")
        root.geometry("1000x780")

        self.project_dir = StringVar()
        self.data_dir = StringVar()
        self.mzrt_list_path = StringVar()  # optional in-house m/z+RT library CSV
        self.param_vars: dict[str, StringVar | BooleanVar] = {}
        self.sample_rows: list[dict] = []   # populated when data folder is scanned
        self.proc: subprocess.Popen | None = None
        self.log_queue: queue.Queue[str] = queue.Queue()

        self._build_layout()
        self._reset_to_defaults()
        self.root.after(100, self._drain_log_queue)

    # ---- layout -----------------------------------------------------------
    def _build_layout(self) -> None:
        root = self.root
        root.columnconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)

        # Top: folder pickers
        top = ttk.LabelFrame(root, text="Folders", padding=8)
        top.grid(row=0, column=0, sticky=(N, S, E, W), padx=8, pady=(8, 4))
        top.columnconfigure(1, weight=1)

        ttk.Label(top, text="Project folder:").grid(row=0, column=0, sticky=W)
        ttk.Entry(top, textvariable=self.project_dir).grid(
            row=0, column=1, sticky=(E, W), padx=4)
        ttk.Button(top, text="Browse…",
                   command=self._pick_project_dir).grid(row=0, column=2)

        ttk.Label(top, text="Raw data folder:").grid(row=1, column=0, sticky=W, pady=(4, 0))
        ttk.Entry(top, textvariable=self.data_dir).grid(
            row=1, column=1, sticky=(E, W), padx=4, pady=(4, 0))
        ttk.Button(top, text="Browse…",
                   command=self._pick_data_dir).grid(row=1, column=2, pady=(4, 0))

        ttk.Label(top, text="m/z + RT library (optional):").grid(
            row=2, column=0, sticky=W, pady=(4, 0))
        ttk.Entry(top, textvariable=self.mzrt_list_path).grid(
            row=2, column=1, sticky=(E, W), padx=4, pady=(4, 0))
        ttk.Button(top, text="Browse…",
                   command=self._pick_mzrt_list).grid(row=2, column=2, pady=(4, 0))

        ttk.Label(top,
                  text="Project folder receives outputs. Raw data folder should hold "
                       ".mzML / .mzXML files; if it isn't <project>/data, the GUI will "
                       "create a junction. The m/z+RT library is a CSV with columns "
                       "name, mz, rt (min) — copied to <project>/mzrt_list.csv on run.",
                  foreground="gray", wraplength=900, justify="left").grid(
            row=3, column=0, columnspan=3, sticky=W, pady=(4, 0))

        # Middle: notebook with Parameters / Samples / Run
        nb = ttk.Notebook(root)
        nb.grid(row=1, column=0, sticky=(N, S, E, W), padx=8, pady=4)
        self._build_parameters_tab(nb)
        self._build_samples_tab(nb)
        self._build_run_tab(nb)

        # Bottom: action buttons
        btns = ttk.Frame(root)
        btns.grid(row=2, column=0, sticky=(E, W), padx=8, pady=(0, 8))
        ttk.Button(btns, text="Reset parameters to defaults",
                   command=self._reset_to_defaults).pack(side="left")
        ttk.Button(btns, text="Load parameters.csv",
                   command=self._load_parameters_csv).pack(side="left", padx=(8, 0))
        ttk.Button(btns, text="Save parameters.csv only",
                   command=self._save_parameters_only).pack(side="left", padx=(8, 0))
        ttk.Button(btns, text="Export MSP",
                   command=self._export_msp_now).pack(side="left", padx=(8, 0))
        ttk.Button(btns, text="Export for formula ID",
                   command=self._export_formula_id_now).pack(side="left", padx=(8, 0))
        ttk.Button(btns, text="Run workflow",
                   command=self._run_workflow).pack(side="right")

    def _build_parameters_tab(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb, padding=4)
        nb.add(frame, text="Parameters")
        inner_nb = ttk.Notebook(frame)
        inner_nb.pack(fill="both", expand=True)

        for group_name, items in PARAM_GROUPS.items():
            tab = ttk.Frame(inner_nb, padding=8)
            inner_nb.add(tab, text=group_name)
            tab.columnconfigure(1, weight=1)
            for row, (pname, (default, kind, choices, helptext)) in enumerate(items):
                ttk.Label(tab, text=pname).grid(row=row, column=0, sticky=W,
                                                pady=2, padx=(0, 8))
                self._make_param_widget(tab, pname, kind, choices, row)
                ttk.Label(tab, text=helptext, foreground="gray").grid(
                    row=row, column=2, sticky=W, padx=(8, 0))

    def _make_param_widget(self, parent, pname, kind, choices, row):
        if kind == "bool":
            var = BooleanVar()
            ttk.Checkbutton(parent, variable=var).grid(row=row, column=1, sticky=W)
        elif kind == "choice":
            var = StringVar()
            ttk.Combobox(parent, textvariable=var, values=choices,
                         state="readonly", width=24).grid(
                row=row, column=1, sticky=W)
        elif kind == "path":
            var = StringVar()
            holder = ttk.Frame(parent)
            holder.grid(row=row, column=1, sticky=(E, W))
            holder.columnconfigure(0, weight=1)
            ttk.Entry(holder, textvariable=var).grid(row=0, column=0, sticky=(E, W))
            ttk.Button(holder, text="…", width=3,
                       command=lambda v=var: self._pick_file_into(v)).grid(
                row=0, column=1, padx=(4, 0))
            if pname == "ms2_library_path":
                ttk.Button(holder, text="Check / Fix", width=11,
                           command=lambda v=var: self._check_ms2_library(v.get())
                           ).grid(row=0, column=2, padx=(4, 0))
        else:  # float, int, str
            var = StringVar()
            ttk.Entry(parent, textvariable=var, width=24).grid(
                row=row, column=1, sticky=W)
        self.param_vars[pname] = var

    def _build_samples_tab(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb, padding=8)
        nb.add(frame, text="Samples")
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(1, weight=1)

        bar = ttk.Frame(frame)
        bar.grid(row=0, column=0, sticky=(E, W))
        ttk.Button(bar, text="Scan raw data folder",
                   command=self._scan_data_folder).pack(side="left")
        ttk.Label(bar, text="  Double-click is_blank / is_qc / group to edit.",
                  foreground="gray").pack(side="left", padx=(8, 0))

        cols = ("sample_name", "is_blank", "is_qc", "group")
        self.tree = ttk.Treeview(frame, columns=cols, show="headings", height=12)
        for c, w in zip(cols, (380, 90, 90, 160)):
            self.tree.heading(c, text=c)
            self.tree.column(c, width=w, anchor=W)
        self.tree.grid(row=1, column=0, sticky=(N, S, E, W), pady=(8, 0))
        sb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        sb.grid(row=1, column=1, sticky=(N, S))
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.bind("<Double-1>", self._edit_cell)

    def _build_run_tab(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb, padding=8)
        nb.add(frame, text="Run log")
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.log = scrolledtext.ScrolledText(frame, wrap="word", height=24,
                                             font=("Consolas", 9))
        self.log.grid(row=0, column=0, sticky=(N, S, E, W))
        self.log.configure(state="disabled")

    # ---- folder pickers ---------------------------------------------------
    def _pick_project_dir(self):
        d = filedialog.askdirectory(title="Select project folder (outputs go here)")
        if d:
            self.project_dir.set(d)

    def _pick_data_dir(self):
        d = filedialog.askdirectory(title="Select folder with raw .mzML / .mzXML files")
        if d:
            self.data_dir.set(d)
            self._scan_data_folder()

    def _pick_file_into(self, var: StringVar):
        f = filedialog.askopenfilename(
            title="Select MS2 library",
            filetypes=[("MS2 library", "*.msp *.pickle"), ("All files", "*.*")])
        if f:
            var.set(f)

    # ---- MS2 library check / auto-fix ------------------------------------
    def _check_ms2_library(self, path_str: str, silent_when_ok: bool = False) -> bool:
        """
        Validate the MS2 library and offer to fix it. Returns True if the
        library is OK after this call (or empty/path missing was acceptable),
        False otherwise. Two classes of issues are detected:

          - Entries with no precursor-m/z field MassCube recognizes (hard
            failure — annotation crashes).
          - Entries with an Ion_mode value MassCube won't match (silent
            failure — annotation runs but yields nothing).
        """
        path_str = (path_str or "").strip()
        if not path_str:
            if not silent_when_ok:
                messagebox.showinfo("MS2 library",
                    "No MS2 library set — annotation step will be skipped.")
            return True
        p = Path(path_str)
        if not p.is_file():
            messagebox.showerror("MS2 library", f"File not found:\n{p}")
            return False

        result = validate_ms2_library(p)
        if "error" in result:
            messagebox.showerror("MS2 library", result["error"])
            return False

        total = result["total"]
        is_msp = result["format"] == "msp"
        bad = result["bad"]
        fixable_hdrs = result["fixable"]
        bad_ion = result.get("bad_ion_mode", [])
        bad_ion_unfixable = result.get("bad_ion_mode_unfixable", [])
        missing_ion = result.get("missing_ion_mode", [])
        bad_decimal = result.get("bad_decimal", [])
        peak_comma_lines = result.get("decimal_comma_peak_lines", 0)

        if (not bad and not bad_ion and not missing_ion
                and not bad_decimal and not peak_comma_lines):
            if not silent_when_ok:
                messagebox.showinfo(
                    "MS2 library OK",
                    f"{total} entries — all have precursor m/z, a "
                    f"MassCube-recognized Ion_mode value, and dot-decimal numbers.")
            return True

        # Build issue report + fix plan.
        issues: list[str] = []
        fixes: list[str] = []

        bad_unfixable = result.get("bad_unfixable_count", 0)
        bad_ms1 = result.get("bad_ms1_count", 0)

        if bad:
            sample = "\n".join(f"    #{i}: {name}" for i, name in bad[:3])
            if len(bad) > 3:
                sample += f"\n    …and {len(bad) - 3} more"
            note = ""
            if bad_unfixable:
                ms1_part = (f" ({bad_ms1} tagged Spectrum_type: MS1 — in-source "
                            f"CID scans MassCube can't use)") if bad_ms1 else ""
                note = (f"\n  {bad_unfixable} of these have no renameable "
                        f"header either — only fix is to drop them"
                        f"{ms1_part}.")
            issues.append(
                f"• {len(bad)} of {total} entries have no precursor m/z field "
                f"MassCube recognizes (needs a key containing both 'prec' and "
                f"'mz'). MassCube will crash on these.{note}\n{sample}")
            if fixable_hdrs and is_msp:
                fixes.append(
                    f"  - rename header(s) {', '.join(fixable_hdrs)} → PRECURSORMZ")
            if bad_unfixable and is_msp:
                fixes.append(
                    f"  - remove {bad_unfixable} entries with no precursor m/z")

        if bad_ion:
            unique_vals = sorted({v for _, _, v in bad_ion})
            issues.append(
                f"• {len(bad_ion)} of {total} entries have an Ion_mode value "
                f"MassCube doesn't recognize. Values found: "
                f"{', '.join(repr(v) for v in unique_vals)}\n"
                f"  (MassCube requires the value to lowercase to 'positive' "
                f"or 'negative'; everything else silently excludes the entry "
                f"from MS2 matching — so the workflow runs but no annotations "
                f"appear.)")
            if is_msp and not bad_ion_unfixable:
                mapping = ", ".join(
                    f"{v!r}→{_ion_mode_fix(v)!r}" for v in unique_vals)
                fixes.append(f"  - normalize Ion_mode values: {mapping}")

        if missing_ion:
            issues.append(
                f"• {len(missing_ion)} of {total} entries have no Ion_mode "
                f"field at all. These will also be excluded from matching "
                f"and can't be auto-fixed (the original polarity is unknown).")

        if bad_decimal or peak_comma_lines:
            parts = []
            if bad_decimal:
                sample = "; ".join(
                    f"#{i} {name}: '{v}'" for i, name, v in bad_decimal[:2])
                parts.append(
                    f"{len(bad_decimal)} entries have a comma-decimal precursor "
                    f"m/z value (e.g. {sample}). MassCube parses this with "
                    f"float() and will crash on the first such entry.")
            if peak_comma_lines:
                parts.append(
                    f"{peak_comma_lines} peak-list line(s) also contain "
                    f"comma-decimal numbers — these would crash MS2 matching "
                    f"the same way.")
            issues.append("• " + " ".join(parts) + " (Locale issue: the library "
                          "was exported with ',' as the decimal separator "
                          "instead of '.'.)")
            if is_msp:
                fixes.append(
                    "  - swap ',' → '.' in precursor m/z headers and peak-list "
                    "lines (free-text fields like Synon/SMILES/InChI/Comment "
                    "are left alone)")

        can_fix_precursor = bool(bad) and bool(fixable_hdrs) and is_msp
        can_fix_ion = bool(bad_ion) and not bad_ion_unfixable and is_msp
        can_strip_bad = bool(bad_unfixable) and is_msp
        can_fix_decimal = bool(bad_decimal or peak_comma_lines) and is_msp

        # Hard precursor failures crash MassCube — block if we can't fix them.
        # Ion-mode issues (missing field, or unrecognized value we can't map)
        # only silently exclude entries from matching, so warn + let the user
        # decide whether to continue.
        if not (can_fix_precursor or can_fix_ion or can_strip_bad or can_fix_decimal):
            extra = ""
            if not is_msp:
                extra = ("\n\nAuto-fix is only supported for .msp files. "
                         "Regenerate or convert the library to fix these.")
            hard_block = bool(bad) or bool(bad_decimal) or bool(peak_comma_lines)
            if hard_block:
                messagebox.showerror(
                    "MS2 library has issues",
                    "Issues found:\n\n" + "\n\n".join(issues) + extra)
                return False
            return messagebox.askokcancel(
                "MS2 library has issues",
                "Issues found:\n\n" + "\n\n".join(issues) + extra
                + "\n\nThese only reduce annotation coverage — MassCube "
                  "will still run, but the affected entries won't match. "
                  "Continue?")

        # Offer the fix.
        msg = (
            f"Issues in\n{p}:\n\n"
            + "\n\n".join(issues)
            + "\n\nProposed fixes:\n"
            + "\n".join(fixes)
            + f"\n\nA .msp.bak backup will be written next to the file."
            + "\n\nApply these fixes now?"
        )
        if not messagebox.askyesno("MS2 library needs fixing", msg):
            return False

        nh = ni = ns = nd = 0
        if can_fix_precursor:
            nh = fix_msp_headers(p)
        if can_fix_ion:
            ni = fix_msp_ion_mode(p)
        if can_strip_bad:
            ns = strip_msp_bad_entries(p)
        if can_fix_decimal:
            nd = fix_msp_decimal_commas(p)

        after = validate_ms2_library(p)
        if after.get("error"):
            messagebox.showerror("MS2 library", after["error"])
            return False

        remaining: list[str] = []
        if after["bad"]:
            remaining.append(
                f"{len(after['bad'])} entries still missing precursor m/z")
        if after.get("bad_ion_mode"):
            remaining.append(
                f"{len(after['bad_ion_mode'])} entries still have an "
                f"unrecognized Ion_mode value")
        if after.get("missing_ion_mode"):
            remaining.append(
                f"{len(after['missing_ion_mode'])} entries have no Ion_mode "
                f"field (silently excluded from matching)")
        if after.get("bad_decimal"):
            remaining.append(
                f"{len(after['bad_decimal'])} entries still have a "
                f"comma-decimal precursor m/z value")
        if after.get("decimal_comma_peak_lines"):
            remaining.append(
                f"{after['decimal_comma_peak_lines']} peak line(s) still "
                f"contain comma-decimal numbers")

        applied = []
        if nh:
            applied.append(f"renamed {nh} header line(s)")
        if ni:
            applied.append(f"normalized {ni} Ion_mode value(s)")
        if ns:
            applied.append(f"removed {ns} entries without precursor m/z")
        if nd:
            applied.append(f"converted {nd} decimal-comma number(s) to dot")
        applied_str = ", ".join(applied) if applied else "no changes were needed"
        backup_note = (
            f" Backup saved to {p.with_suffix(p.suffix + '.bak').name}."
            if (nh or ni or ns or nd) else "")

        if remaining:
            messagebox.showwarning(
                "MS2 library still has issues",
                f"Applied: {applied_str}.{backup_note}\n\n"
                f"Remaining:\n  - " + "\n  - ".join(remaining)
                + "\n\nYou'll need to edit the file by hand for those.")
            # Hard precursor / decimal-comma failures block the run; ion_mode
            # issues don't crash MassCube but reduce annotation coverage —
            # treat those as warnings.
            return not (after["bad"] or after.get("bad_decimal")
                        or after.get("decimal_comma_peak_lines"))

        messagebox.showinfo(
            "MS2 library fixed",
            f"Applied: {applied_str}.{backup_note}")
        return True

    def _pick_mzrt_list(self):
        f = filedialog.askopenfilename(
            title="Select m/z + RT library (CSV: name, mz, rt)",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if f:
            self.mzrt_list_path.set(f)

    # ---- defaults / params I/O -------------------------------------------
    def _reset_to_defaults(self):
        for pname, (default, kind, _, _) in PARAM_SPECS.items():
            var = self.param_vars[pname]
            if kind == "bool":
                var.set(bool(default))
            elif kind == "path":
                var.set("" if default in (None, "") else str(default))
            else:
                var.set("" if default is None else str(default))

    def _collect_params(self) -> list[tuple[str, str]]:
        rows: list[tuple[str, str]] = []
        for pname, (default, kind, _, _) in PARAM_SPECS.items():
            if pname in GUI_ONLY_PARAMS:
                continue  # GUI-side post-processing, not a MassCube parameter
            var = self.param_vars[pname]
            if kind == "bool":
                rows.append((pname, "True" if var.get() else "False"))
                continue
            raw = var.get().strip()
            if kind == "path" and not raw:
                # leave ms2_library_path blank — workflow will skip annotation
                continue
            if kind in ("float", "int"):
                if raw == "":
                    raise ValueError(f"Parameter '{pname}' is required and was left blank.")
                try:
                    if kind == "int":
                        int(float(raw))
                    else:
                        float(raw)
                except ValueError:
                    raise ValueError(f"Parameter '{pname}' must be numeric (got '{raw}').")
            rows.append((pname, raw))
        return rows

    def _write_parameters_csv(self, project: Path) -> Path:
        rows = self._collect_params()
        out = project / "parameters.csv"
        with out.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["parameter", "value"])
            w.writerows(rows)
        return out

    def _save_parameters_only(self):
        try:
            project = self._validated_project_dir()
            out = self._write_parameters_csv(project)
            messagebox.showinfo("Saved", f"Wrote {out}")
        except Exception as exc:
            messagebox.showerror("Error", str(exc))

    def _load_parameters_csv(self):
        proj = self.project_dir.get().strip()
        initial_dir = proj if proj and os.path.isdir(proj) else None
        path = filedialog.askopenfilename(
            title="Load parameters.csv",
            initialdir=initial_dir,
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", newline="", encoding="utf-8-sig") as f:
                reader = csv.reader(f)
                rows = [r for r in reader if r and any(c.strip() for c in r)]
        except OSError as exc:
            messagebox.showerror("Load parameters.csv", f"Could not read file:\n{exc}")
            return

        if not rows:
            messagebox.showerror("Load parameters.csv", "File is empty.")
            return

        # Skip header if present.
        header = [c.strip().lower() for c in rows[0]]
        if header[:2] == ["parameter", "value"]:
            rows = rows[1:]

        applied: list[str] = []
        unknown: list[str] = []
        invalid: list[str] = []

        for row in rows:
            if len(row) < 2:
                continue
            pname = row[0].strip()
            raw = row[1].strip()
            if not pname:
                continue
            spec = PARAM_SPECS.get(pname)
            if spec is None:
                unknown.append(pname)
                continue
            _default, kind, choices, _help = spec
            var = self.param_vars[pname]
            try:
                if kind == "bool":
                    low = raw.lower()
                    if low in ("true", "1", "yes", "y"):
                        var.set(True)
                    elif low in ("false", "0", "no", "n", ""):
                        var.set(False)
                    else:
                        raise ValueError(f"expected True/False, got '{raw}'")
                elif kind == "choice":
                    if raw and choices and raw not in choices:
                        raise ValueError(
                            f"'{raw}' not in allowed values {choices}")
                    var.set(raw)
                elif kind == "int":
                    if raw != "":
                        int(float(raw))  # validate
                    var.set(raw)
                elif kind == "float":
                    if raw != "":
                        float(raw)  # validate
                    var.set(raw)
                else:  # str, path
                    var.set(raw)
                applied.append(pname)
            except ValueError as exc:
                invalid.append(f"{pname}: {exc}")

        msg_parts = [f"Loaded {len(applied)} parameter(s) from:\n{path}"]
        if unknown:
            msg_parts.append(
                f"\nIgnored {len(unknown)} unknown parameter(s):\n  - "
                + "\n  - ".join(unknown))
        if invalid:
            msg_parts.append(
                f"\nSkipped {len(invalid)} invalid value(s):\n  - "
                + "\n  - ".join(invalid))
        if unknown or invalid:
            messagebox.showwarning("Load parameters.csv", "\n".join(msg_parts))
        else:
            messagebox.showinfo("Load parameters.csv", "\n".join(msg_parts))

    # ---- sample table -----------------------------------------------------
    def _scan_data_folder(self):
        d = self.data_dir.get().strip()
        if not d or not os.path.isdir(d):
            messagebox.showwarning("No data folder",
                                   "Pick the raw data folder first.")
            return
        files = sorted(
            f for f in os.listdir(d)
            if not f.startswith(".") and f.lower().endswith(RAW_EXTS)
        )
        if not files:
            messagebox.showwarning("No raw files",
                                   f"No .mzML / .mzXML files in {d}.")
            return
        self.tree.delete(*self.tree.get_children())
        self.sample_rows.clear()
        for fname in files:
            name = os.path.splitext(fname)[0]
            row = {"sample_name": name, "is_blank": "no", "is_qc": "no", "group": ""}
            self.sample_rows.append(row)
            self.tree.insert("", END, values=tuple(row.values()))

    def _edit_cell(self, event):
        region = self.tree.identify("region", event.x, event.y)
        if region != "cell":
            return
        col = self.tree.identify_column(event.x)
        item = self.tree.identify_row(event.y)
        if not item:
            return
        col_idx = int(col.replace("#", "")) - 1
        col_name = ("sample_name", "is_blank", "is_qc", "group")[col_idx]
        cur = self.tree.set(item, col_name)
        if col_name in ("is_blank", "is_qc"):
            new = "yes" if cur != "yes" else "no"
        elif col_name == "group":
            new = _prompt_string(self.root, "Group", "Group label:", cur)
            if new is None:
                return
        else:
            return  # sample_name is read-only
        self.tree.set(item, col_name, new)
        idx = self.tree.index(item)
        self.sample_rows[idx][col_name] = new

    def _sync_mzrt_list(self, project: Path) -> None:
        """Copy the chosen m/z+RT library to <project>/mzrt_list.csv (or clean up)."""
        target = project / "mzrt_list.csv"
        src = self.mzrt_list_path.get().strip()
        if not src:
            return  # leave existing file alone
        src_path = Path(src)
        if not src_path.is_file():
            raise ValueError(f"m/z + RT library not found: {src}")
        try:
            same = target.exists() and src_path.resolve() == target.resolve()
        except OSError:
            same = False
        if not same:
            shutil.copyfile(src_path, target)

    def _write_sample_table(self, project: Path) -> Path | None:
        if not self.sample_rows:
            return None
        out = project / "sample_table.csv"
        with out.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["sample_name", "is_blank", "is_qc", "group"])
            for r in self.sample_rows:
                w.writerow([r["sample_name"], r["is_blank"], r["is_qc"], r["group"]])
        return out

    # ---- data folder linking ---------------------------------------------
    def _ensure_data_folder(self, project: Path, data: Path) -> None:
        """Make sure <project>/data points at the raw-data folder."""
        target = project / "data"
        try:
            if target.resolve() == data.resolve():
                return
        except OSError:
            pass
        if target.exists():
            if not messagebox.askyesno(
                "Replace existing data folder?",
                f"{target} already exists. Replace it with a link to {data}?",
            ):
                raise RuntimeError("Aborted: existing data folder kept.")
            if target.is_symlink() or _is_junction(target):
                target.unlink()
            else:
                shutil.rmtree(target)
        # Junction works on Windows without admin rights
        if sys.platform.startswith("win"):
            subprocess.check_call(
                ["cmd", "/c", "mklink", "/J", str(target), str(data)],
                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
            )
        else:
            os.symlink(data, target, target_is_directory=True)

    # ---- run --------------------------------------------------------------
    def _validated_project_dir(self) -> Path:
        p = self.project_dir.get().strip()
        if not p:
            raise ValueError("Pick a project folder first.")
        project = Path(p)
        project.mkdir(parents=True, exist_ok=True)
        return project

    def _run_workflow(self):
        if self.proc and self.proc.poll() is None:
            messagebox.showinfo("Already running", "The workflow is still running.")
            return
        try:
            project = self._validated_project_dir()
            d = self.data_dir.get().strip()
            if not d or not os.path.isdir(d):
                raise ValueError("Pick the raw data folder.")
            ms2_path = self.param_vars["ms2_library_path"].get().strip()
            if ms2_path and not self._check_ms2_library(ms2_path, silent_when_ok=True):
                return  # validator already showed an error/cancellation dialog
            self._ensure_data_folder(project, Path(d))
            self._write_parameters_csv(project)
            self._write_sample_table(project)
            self._sync_mzrt_list(project)
        except Exception as exc:
            messagebox.showerror("Cannot start", str(exc))
            return

        # Read the MSP settings here, on the main thread — Tk variables must not
        # be touched from the reader thread that does the post-run export.
        msp_opts = self._msp_options()

        self._log_clear()
        self._log(f"[GUI] project: {project}\n")
        self._log(f"[GUI] data:    {d}\n")
        self._log("[GUI] launching MassCube…\n\n")

        script = (
            "import sys\n"
            "try:\n"
            "    from masscube.workflows import untargeted_metabolomics_workflow\n"
            "except ImportError:\n"
            "    sys.stderr.write('masscube is not installed. Run: pip install masscube\\n')\n"
            "    sys.exit(2)\n"
            f"untargeted_metabolomics_workflow(path=r'{project}')\n"
        )
        # MassCube 1.2.x prints box-drawing characters (e.g. U+2550) in its
        # banner. On Windows the subprocess stdout defaults to cp1252 and
        # crashes on those chars — force UTF-8 in the child + on our reader.
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"
        env["PYTHONUTF8"] = "1"
        self.proc = subprocess.Popen(
            [sys.executable, "-u", "-c", script],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
            encoding="utf-8", errors="replace",
            env=env,
        )
        threading.Thread(target=self._reader_thread,
                         args=(self.proc, project, msp_opts), daemon=True).start()

    def _reader_thread(self, proc: subprocess.Popen, project: Path, msp_opts: dict):
        assert proc.stdout is not None
        for line in proc.stdout:
            self.log_queue.put(line)
        proc.wait()
        self.log_queue.put(f"\n[GUI] workflow exited with code {proc.returncode}\n")
        if proc.returncode == 0 and msp_opts["enabled"]:
            self._export_msp_after_run(project, msp_opts)

    # ---- MSP export -------------------------------------------------------
    def _msp_options(self) -> dict:
        """Snapshot the MSP-export settings (must run on the main thread)."""
        content = self.param_vars["msp_content"].get().strip() or MSP_CONTENT_CHOICES[0]
        return {
            "enabled": bool(self.param_vars["output_msp"].get()),
            "content": content,
            "ion_mode": self.param_vars["ion_mode"].get().strip(),
        }

    def _export_msp_after_run(self, project: Path, msp_opts: dict) -> None:
        """Post-run MSP export. Runs on the reader thread — log only, no dialogs."""
        self.log_queue.put(f"\n[GUI] exporting MSP ({msp_opts['content']})…\n")
        table = find_feature_table(project)
        if table is None:
            self.log_queue.put(
                "[GUI] no aligned feature table found — MSP export skipped.\n")
            return
        msp_path = project / MSP_FILE_NAME
        try:
            stats = write_msp_from_feature_table(
                table, msp_path,
                ion_mode=msp_opts["ion_mode"], content=msp_opts["content"])
        except Exception as exc:
            self.log_queue.put(f"[GUI] MSP export failed: {exc}\n")
            return
        for line in summarize_msp_stats(stats, msp_path, table).splitlines():
            self.log_queue.put(f"[GUI] {line}\n")

    def _export_msp_now(self):
        """Export the MSP for the selected project without re-running MassCube."""
        proj = self.project_dir.get().strip()
        if not proj or not os.path.isdir(proj):
            messagebox.showwarning(
                "Export MSP", "Pick an existing project folder first.")
            return
        project = Path(proj)
        table = find_feature_table(project)
        if table is None:
            messagebox.showerror(
                "Export MSP",
                "No aligned feature table found in:\n"
                f"{project}\n\nRun the workflow first (the exporter reads "
                "aligned_feature_table.txt).")
            return
        opts = self._msp_options()
        msp_path = project / MSP_FILE_NAME
        try:
            stats = write_msp_from_feature_table(
                table, msp_path, ion_mode=opts["ion_mode"], content=opts["content"])
        except Exception as exc:
            messagebox.showerror("Export MSP", f"Export failed:\n{exc}")
            return
        summary = summarize_msp_stats(stats, msp_path, table)
        self._log(f"\n[GUI] MSP export ({opts['content']})\n")
        for line in summary.splitlines():
            self._log(f"[GUI] {line}\n")
        if stats["written"]:
            messagebox.showinfo("Export MSP", summary)
        else:
            messagebox.showwarning(
                "Export MSP",
                "No entries matched the selected msp_content option "
                f"({opts['content']}).\n\n{summary}")

    # ---- formula-ID export ------------------------------------------------
    def _export_formula_id_now(self):
        """Write SIRIUS .ms input + an MS1-bearing MSP for MS-FINDER."""
        proj = self.project_dir.get().strip()
        if not proj or not os.path.isdir(proj):
            messagebox.showwarning(
                "Export for formula ID", "Pick an existing project folder first.")
            return
        project = Path(proj)
        table = find_feature_table(project)
        if table is None:
            messagebox.showerror(
                "Export for formula ID",
                "No aligned feature table found in:\n"
                f"{project}\n\nRun the workflow first (the exporter reads "
                "aligned_feature_table.txt).")
            return

        ion_mode = self.param_vars["ion_mode"].get().strip()
        sirius_dir = project / SIRIUS_DIR_NAME
        msp_path = project / MSFINDER_MSP_FILE_NAME

        self._log("\n[GUI] formula-ID export\n")
        try:
            sirius_stats = write_sirius_ms_from_feature_table(
                table, sirius_dir, ion_mode=ion_mode)
            msp_stats = write_msp_from_feature_table(
                table, msp_path, ion_mode=ion_mode,
                content="all_features", isotope_ms1=True)
        except Exception as exc:
            messagebox.showerror("Export for formula ID", f"Export failed:\n{exc}")
            return

        summary = "\n".join([
            "SIRIUS (.ms):",
            summarize_formula_export_stats(sirius_stats, sirius_dir, table),
            "",
            "MS-FINDER (.msp, MS1 + MS2):",
            summarize_msp_stats(msp_stats, msp_path, table),
        ])
        for line in summary.splitlines():
            self._log(f"[GUI] {line}\n")

        if sirius_stats["written"]:
            messagebox.showinfo("Export for formula ID", summary)
        else:
            messagebox.showwarning(
                "Export for formula ID",
                "No features were exported — every row was filtered out.\n\n"
                + summary)

    # ---- log helpers ------------------------------------------------------
    def _drain_log_queue(self):
        try:
            while True:
                line = self.log_queue.get_nowait()
                self._log(line)
        except queue.Empty:
            pass
        self.root.after(100, self._drain_log_queue)

    def _log(self, text: str):
        self.log.configure(state="normal")
        self.log.insert(END, text)
        self.log.see(END)
        self.log.configure(state="disabled")

    def _log_clear(self):
        self.log.configure(state="normal")
        self.log.delete("1.0", END)
        self.log.configure(state="disabled")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _is_junction(p: Path) -> bool:
    if not sys.platform.startswith("win") or not p.exists():
        return False
    try:
        import ctypes
        FILE_ATTRIBUTE_REPARSE_POINT = 0x0400
        attrs = ctypes.windll.kernel32.GetFileAttributesW(str(p))
        return attrs != -1 and bool(attrs & FILE_ATTRIBUTE_REPARSE_POINT)
    except Exception:
        return False


def _prompt_string(parent, title, prompt, initial=""):
    """Tiny inline replacement for simpledialog.askstring (keeps imports tidy)."""
    from tkinter import simpledialog
    return simpledialog.askstring(title, prompt, initialvalue=initial, parent=parent)


def main():
    root = Tk()
    MassCubeGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
