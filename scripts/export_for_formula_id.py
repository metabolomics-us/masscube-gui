"""Export a MassCube project for molecular-formula assignment.

Turns aligned_feature_table.txt into input for a formula-ID tool that scores
candidate formulas against accurate m/z *and* MS1 isotope abundance:

  SIRIUS      one .ms per feature (default) or a single multi-compound .ms,
              with >parentmass, >charge, >ionization, >ms1peaks, >ms2peaks
  MS-FINDER   an .msp carrying MSTYPE: MS1 (isotope pattern) + MSTYPE: MS2

Isotope-peak rows and in-source fragments are dropped, multimer adducts
([2M+H]+ etc.) are skipped, and no formula is ever written into the output --
the point is for the tool to derive it.

Examples
--------
    python scripts/export_for_formula_id.py "D:/projects/run7"
    python scripts/export_for_formula_id.py "D:/projects/run7" --format sirius --single-file
    python scripts/export_for_formula_id.py "D:/projects/run7" --ion-mode negative

Then in SIRIUS: File > Import > select the sirius_input folder, or on the CLI
    sirius -i sirius_input -o sirius_out formulas
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Reuse the exporters and feature-table handling that back the GUI button,
# rather than reimplementing the parsing here.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from masscube_gui import (  # noqa: E402
    MSFINDER_MSP_FILE_NAME,
    SIRIUS_DIR_NAME,
    SIRIUS_SINGLE_FILE_NAME,
    find_feature_table,
    summarize_formula_export_stats,
    summarize_msp_stats,
    write_msp_from_feature_table,
    write_sirius_ms_from_feature_table,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Reads aligned_feature_table.txt, or falls back to "
               "project_files/aligned_feature_table_before_annotation.txt.")
    ap.add_argument("project", type=Path, help="MassCube project folder")
    ap.add_argument("--format", choices=("both", "sirius", "msfinder"),
                    default="both", help="output to write (default: both)")
    ap.add_argument("--ion-mode", default="",
                    help="positive/negative; used only for rows with no adduct")
    ap.add_argument("--single-file", action="store_true",
                    help="write one multi-compound sirius_input.ms instead of a "
                         "directory of per-compound files")
    ap.add_argument("--require-isotopes", action="store_true",
                    help="skip features with fewer than 2 isotope peaks")
    ap.add_argument("--keep-in-source-fragments", action="store_true",
                    help="keep rows flagged is_in_source_fragment (off by "
                         "default: fragments are not intact molecules)")
    ap.add_argument("--no-ms2", action="store_true",
                    help="omit MS2 from the SIRIUS export (mass + isotopes only)")
    args = ap.parse_args(argv)

    project: Path = args.project
    if not project.is_dir():
        return err(f"Not a directory: {project}")
    table = find_feature_table(project)
    if table is None:
        return err(f"No aligned feature table found in {project}\n"
                   "Run the MassCube workflow first.")

    print(f"Source: {table}")

    if args.format in ("both", "sirius"):
        out_path = (project / SIRIUS_SINGLE_FILE_NAME if args.single_file
                    else project / SIRIUS_DIR_NAME)
        stats = write_sirius_ms_from_feature_table(
            table, out_path,
            ion_mode=args.ion_mode,
            single_file=args.single_file,
            include_ms2=not args.no_ms2,
            require_isotopes=args.require_isotopes,
            skip_in_source_fragments=not args.keep_in_source_fragments,
        )
        print("\nSIRIUS (.ms):")
        print(summarize_formula_export_stats(stats, out_path, table))

    if args.format in ("both", "msfinder"):
        msp_path = project / MSFINDER_MSP_FILE_NAME
        stats = write_msp_from_feature_table(
            table, msp_path, ion_mode=args.ion_mode,
            content="all_features", isotope_ms1=True)
        print("\nMS-FINDER (.msp, MS1 + MS2):")
        print(summarize_msp_stats(stats, msp_path, table))

    return 0


def err(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
