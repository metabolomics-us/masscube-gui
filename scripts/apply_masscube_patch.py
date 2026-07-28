"""Apply local fixes to the installed masscube package.

Two patches, both editing masscube in site-packages:

  get-start-time    STEP 1 ("Preparing the workflow...") calls get_start_time()
                    once per raw file and upstream reads the WHOLE file via
                    f.readlines()[:300] to find a timestamp in the header. On
                    GB-sized mzML -- especially in a cloud-synced folder, where
                    the read forces a full download -- STEP 1 becomes a silent
                    multi-minute or multi-hour hang. Reads a bounded 256 KB
                    header chunk instead. Return values unchanged.

  isotope-spacing   find_isotope_signals() searches for isotopes at integer Da
                    offsets. That is deliberate upstream (it keeps 15N +0.99703,
                    34S +1.99580 and 37Cl +1.99705 inside mz_tol), but pure 13C
                    drifts +0.00336 per step, so from M+3 (+3.01006) on it falls
                    outside the default 0.01 tolerance and is silently lost.
                    Searches the integer grid AND the 13C grid -- a superset, so
                    nothing that matched before stops matching. Recovers M+3/M+4
                    for formula assignment from isotope abundance.

Usage:

    python scripts/apply_masscube_patch.py                  # apply all
    python scripts/apply_masscube_patch.py --check           # report only
    python scripts/apply_masscube_patch.py --revert          # restore backups
    python scripts/apply_masscube_patch.py --only get-start-time

Safe to run repeatedly -- already-patched files are detected and skipped. Re-run
after every `pip install --upgrade masscube`, which overwrites site-packages and
silently reverts these fixes.

isotope-spacing only affects newly detected features, so it needs a fresh
workflow run to change anything. get-start-time takes effect immediately.

See patches/ for the corresponding unified diffs and rationale.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

BACKUP_SUFFIX = ".orig"

PATCHES = [
    {
        "name": "get-start-time",
        "module_file": "utils_functions.py",
        "summary": "bounded 256 KB header read instead of reading whole files",
        "original": '''    if os.path.exists(str(file_name)):
        with open(file_name, "rb") as f:
            # check the first 300 rows for the start time
            for l in f.readlines()[:300]:
                l = str(l)
                if "startTimeStamp" in str(l):
                    t = l.split("startTimeStamp")[1].split('"')[1]
                    return datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ")
    return None
''',
        "patched": '''    if os.path.exists(str(file_name)):
        with open(file_name, "rb") as f:
            # Read a bounded header chunk instead of f.readlines(), which pulls the
            # whole file into memory before the [:300] slice — on GB-sized mzML
            # files (or cloud-synced folders) that alone can take minutes per file.
            # 256 KB covers far more than 300 lines of mzML/mzXML header.
            head = f.read(262144)
        # check the first 300 rows for the start time
        for l in head.splitlines()[:300]:
            l = str(l)
            if "startTimeStamp" in str(l):
                t = l.split("startTimeStamp")[1].split('"')[1]
                return datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ")
    return None
''',
    },
    {
        "name": "isotope-spacing",
        "module_file": "feature_grouping.py",
        "summary": "search integer AND 13C offset grids, recovering M+3/M+4",
        "original": "    targets = mz + np.arange(num, dtype=mzs.dtype)\n",
        "patched": '''    # Search two offset grids, because no single spacing covers every isotope:
    #   integer steps  -- keeps the sub-1 Da isotopes in range (15N +0.99703,
    #                     34S +1.99580, 37Cl +1.99705), which is why upstream
    #                     uses them; see the note at the top of this module.
    #   13C steps      -- pure carbon drifts +0.00336 per step, so by M+3
    #                     (+3.01006) it falls outside mz_tol of the integer
    #                     grid and the peak is silently lost.
    # The union is a superset of the integer grid, so nothing that matched
    # before stops matching.
    steps = np.arange(num, dtype=mzs.dtype)
    targets = mz + np.concatenate([steps, 1.0033548 * steps])
''',
    },
]


def masscube_dir() -> Path:
    try:
        import masscube
    except ImportError:
        sys.exit("masscube is not installed in this interpreter. "
                 "Run: pip install masscube")
    return Path(masscube.__file__).parent


def process(patch: dict, root: Path, mode: str) -> str:
    """Apply/check/revert one patch. Returns 'ok', 'skip', or 'fail'."""
    path = root / patch["module_file"]
    label = f"{patch['name']:16}"

    if not path.is_file():
        print(f"{label} FAIL   missing file: {path}", file=sys.stderr)
        return "fail"

    backup = path.with_name(path.name + BACKUP_SUFFIX)
    text = path.read_text(encoding="utf-8")

    if mode == "revert":
        if not backup.is_file():
            print(f"{label} SKIP   no backup at {backup.name}")
            return "skip"
        shutil.copyfile(backup, path)
        print(f"{label} OK     reverted from {backup.name}")
        return "ok"

    if patch["patched"] in text:
        print(f"{label} SKIP   already patched")
        return "skip"
    if patch["original"] not in text:
        print(f"{label} FAIL   upstream code not found — masscube has probably "
              f"changed.\n{'':17}Check whether the fix is still needed before "
              f"editing {patch['module_file']} by hand.", file=sys.stderr)
        return "fail"

    if mode == "check":
        print(f"{label} TODO   unpatched, upstream code matched "
              f"({patch['summary']})")
        return "ok"

    if not backup.exists():
        shutil.copyfile(path, backup)
    path.write_text(text.replace(patch["original"], patch["patched"]),
                    encoding="utf-8")
    print(f"{label} OK     patched — {patch['summary']}")
    return "ok"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Patches: " + ", ".join(p["name"] for p in PATCHES))
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true",
                      help="report status without writing")
    mode.add_argument("--revert", action="store_true",
                      help="restore the originals from backup")
    ap.add_argument("--only", action="append", metavar="NAME",
                    choices=[p["name"] for p in PATCHES],
                    help="limit to one patch (repeatable)")
    args = ap.parse_args(argv)

    mode_name = "check" if args.check else "revert" if args.revert else "apply"
    root = masscube_dir()
    print(f"masscube: {root}")
    print(f"mode:     {mode_name}\n")

    selected = [p for p in PATCHES if not args.only or p["name"] in args.only]
    results = [process(p, root, mode_name) for p in selected]

    failed = results.count("fail")
    print(f"\n{results.count('ok')} ok, {results.count('skip')} skipped, "
          f"{failed} failed")
    if mode_name == "apply" and not failed:
        print("Re-run this after any `pip install --upgrade masscube`.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
