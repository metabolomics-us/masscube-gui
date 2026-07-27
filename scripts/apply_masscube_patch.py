"""Apply the get_start_time() speed fix to the installed masscube package.

masscube's STEP 1 ("Preparing the workflow...") calls get_start_time() once per
raw file to sort samples by acquisition time, and upstream reads the whole file
via f.readlines()[:300] to find a timestamp in the header. On GB-sized mzML
files -- especially in a cloud-synced folder, where the read forces a full
download -- that turns STEP 1 into a silent multi-minute (or multi-hour) hang.

This script rewrites that one function to read a bounded 256 KB header chunk
instead. Return values are unchanged.

Run it after installing or upgrading masscube:

    python scripts/apply_masscube_patch.py           # apply
    python scripts/apply_masscube_patch.py --check    # report status only
    python scripts/apply_masscube_patch.py --revert   # restore from backup

Safe to run repeatedly -- it detects an already-patched file and does nothing.
See patches/masscube-get-start-time.patch for the diff and rationale.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

BACKUP_SUFFIX = ".orig"

ORIGINAL = '''    if os.path.exists(str(file_name)):
        with open(file_name, "rb") as f:
            # check the first 300 rows for the start time
            for l in f.readlines()[:300]:
                l = str(l)
                if "startTimeStamp" in str(l):
                    t = l.split("startTimeStamp")[1].split('"')[1]
                    return datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ")
    return None
'''

PATCHED = '''    if os.path.exists(str(file_name)):
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
'''


def target_file() -> Path:
    """Locate utils_functions.py inside the installed masscube package."""
    try:
        import masscube
    except ImportError:
        sys.exit("masscube is not installed in this interpreter. Run: pip install masscube")
    path = Path(masscube.__file__).parent / "utils_functions.py"
    if not path.is_file():
        sys.exit(f"Expected file not found: {path}")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="report status without writing")
    mode.add_argument("--revert", action="store_true", help="restore the original from backup")
    args = ap.parse_args()

    path = target_file()
    backup = path.with_name(path.name + BACKUP_SUFFIX)
    text = path.read_text(encoding="utf-8")
    patched = PATCHED in text
    original = ORIGINAL in text

    print(f"target: {path}")

    if args.revert:
        if not backup.is_file():
            return err(f"no backup at {backup} — nothing to revert to")
        shutil.copyfile(backup, path)
        print(f"reverted from {backup.name}")
        return 0

    if patched:
        print("status: already patched — nothing to do")
        return 0
    if not original:
        return err(
            "status: NOT patched, but the expected upstream code was not found either.\n"
            "        masscube's get_start_time() has probably changed upstream — check\n"
            "        whether the fix is still needed before editing by hand.\n"
            "        See patches/masscube-get-start-time.patch"
        )

    print("status: unpatched, upstream code matched")
    if args.check:
        print("(--check given, no changes written)")
        return 0

    if not backup.exists():
        shutil.copyfile(path, backup)
        print(f"backed up original to {backup.name}")
    path.write_text(text.replace(ORIGINAL, PATCHED), encoding="utf-8")
    print("patched successfully")
    return 0


def err(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
