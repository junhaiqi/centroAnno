#!/usr/bin/env python3
"""
continue_pipeline.py — Generic post-processing pipeline for centroAnno outputs.

Given existing centroAnno Stage-1 results (*_decomposedResult.csv), this script
runs Stages 2–5 automatically:

  Stage 2: detectHORs.py       → HOR patterns & statistics
  Stage 3: centroFinder.py     → Centromere candidate prediction
  Stage 4: visualizeMonomerArray.py (basic) → Repeat / HOR region maps
  Stage 5: Refine top-1 BED → anno-sat-asm → detectHORs.py → composite viz

Usage modes
-----------
1. Batch mode (one output root with chromosome subdirectories):
   python3 continue_pipeline.py -d out/ -f ref.fa -t 16

2. Single-result mode (one CSV file):
   python3 continue_pipeline.py -c chr1_decomposedResult.csv -f ref.fa -t 16

3. Batch with explicit FASTA directory (auto-match by basename):
   python3 continue_pipeline.py -d out/ -F fastas/ -t 16

Output structure
----------------
For each input CSV, results are placed alongside it:

  <input_dir>/
  ├── *_decomposedResult.csv           (existing)
  ├── *_monomerTemplates.fa            (existing)
  ├── 02_HORs/
  ├── 03_cenSatArray/
  ├── 04_visualization/
  └── 05_cenSatArray_refined/

Skip options
------------
  --skip-hors          Skip Stage 2 (HOR detection)
  --skip-centrofinder  Skip Stage 3 (CenSatArray detection)
  --skip-visualize     Skip Stage 4 (basic visualization)
  --skip-refine        Skip Stage 5 (centromere refinement)
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg: str):
    import datetime
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [pipeline] {msg}", flush=True)


def run_cmd(cmd: list, cwd=None):
    """Run a command, log it, and return (rc, stdout, stderr)."""
    log(f"  RUN: {' '.join(str(c) for c in cmd)}")
    try:
        proc = subprocess.run(
            [str(c) for c in cmd],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=None,
        )
        if proc.stdout:
            for line in proc.stdout.strip().splitlines():
                log(f"    OUT: {line}")
        if proc.returncode != 0 and proc.stderr:
            for line in proc.stderr.strip().splitlines():
                log(f"    ERR: {line}")
        return proc.returncode, proc.stdout, proc.stderr
    except Exception as e:
        log(f"  EXCEPTION: {e}")
        return 1, "", str(e)


def resolve_scripts(script_dir: Path):
    """Verify that all required helper scripts exist."""
    scripts = {
        "detectHORs": script_dir / "detectHORs.py",
        "centroFinder": script_dir / "centroFinder.py",
        "visualize": script_dir / "visualizeMonomerArray.py",
        "extractBed": script_dir / "extractBedRegions.py",
    }
    for name, path in scripts.items():
        if not path.is_file():
            log(f"[ERROR] Missing required script: {path}")
            sys.exit(1)
    return scripts


def resolve_binary(bin_path: Path):
    if not bin_path.is_file():
        log(f"[ERROR] centroAnno binary not found: {bin_path}")
        sys.exit(1)
    return bin_path


def find_csv_files(root: Path) -> list:
    """Recursively find all *_decomposedResult.csv under root."""
    return sorted(root.rglob("*_decomposedResult.csv"))


def infer_fasta(csv_path: Path, fasta: Path = None, fasta_dir: Path = None) -> Path:
    """
    Infer the source FASTA for a given CSV.

    Priority:
      1. Explicit -f/--fasta (single mode)
      2. Explicit -F/--fasta-dir with basename matching
      3. Parent directory name matching (heuristic)
    """
    if fasta is not None and fasta.is_file():
        return fasta

    base = csv_path.stem.replace("_decomposedResult", "")
    # Try fasta_dir
    if fasta_dir is not None and fasta_dir.is_dir():
        for ext in (".fasta", ".fa", ".fna"):
            candidate = fasta_dir / f"{base}{ext}"
            if candidate.is_file():
                return candidate

    # Heuristic: parent dir name
    parent_name = csv_path.parent.name
    if fasta_dir is not None and fasta_dir.is_dir():
        for ext in (".fasta", ".fa", ".fna"):
            candidate = fasta_dir / f"{parent_name}{ext}"
            if candidate.is_file():
                return candidate

    return None


# ---------------------------------------------------------------------------
# Stage runners
# ---------------------------------------------------------------------------

def stage2_detect_hors(csv_path: Path, scripts: dict, args, out_dir: Path) -> Path:
    """Run detectHORs.py. Returns path to _allHORs.csv or None."""
    hor_dir = out_dir / "02_HORs"
    hor_dir.mkdir(parents=True, exist_ok=True)
    basename = csv_path.stem.replace("_decomposedResult", "")
    prefix = hor_dir / basename

    cmd = [
        args.python,
        scripts["detectHORs"],
        "-i", csv_path,
        "-o", prefix,
        "-m", str(args.max_hor_len),
        "-c", str(args.min_hor_copies),
        "--min-identity", str(args.min_hor_identity),
        "--min-max-copies", str(args.min_max_copies),
    ]
    rc, _, _ = run_cmd(cmd)
    hor_csv = prefix.parent / f"{basename}_allHORs.csv"
    if rc == 0 and hor_csv.is_file():
        log(f"  [OK] HORs → {hor_csv}")
        return hor_csv
    log(f"  [WARN] HOR detection may have failed (rc={rc})")
    return hor_csv if hor_csv.is_file() else None


def stage3_centro_finder(csv_path: Path, fasta: Path, scripts: dict, args, out_dir: Path) -> Path:
    """Run centroFinder.py. Returns path to BED or None."""
    cf_dir = out_dir / "03_cenSatArray"
    cf_dir.mkdir(parents=True, exist_ok=True)
    basename = csv_path.stem.replace("_decomposedResult", "")
    prefix = cf_dir / basename

    cmd = [
        args.python,
        scripts["centroFinder"],
        "--mono-csv", csv_path,
        "--fasta", fasta,
        "--out-prefix", prefix,
        "--window", "50000",
        "--step", "10000",
        "--density-thr", str(args.density_thr),
        "--gap-factor", str(args.gap_factor),
        "--min-size", str(args.min_size),
        "--max-size", str(args.max_size),
        "--min-mono-len", str(args.min_mono_len),
        "--mono-gap-bp", str(args.mono_gap_bp),
        "--identity-boundary-thr", str(args.identity_boundary_thr),
    ]
    if args.no_offset:
        cmd.append("--no-offset")

    rc, _, _ = run_cmd(cmd)
    bed = prefix.parent / f"{basename}_censatarray_candidates.bed"
    if rc == 0 and bed.is_file():
        log(f"  [OK] Centromere BED → {bed}")
        return bed
    log(f"  [WARN] centroFinder may have failed (rc={rc})")
    return bed if bed.is_file() else None


def stage4_visualize(csv_path: Path, hor_csv: Path, bed: Path, scripts: dict, args, out_dir: Path):
    """Run visualizeMonomerArray.py (basic, no composite)."""
    viz_dir = out_dir / "04_visualization"
    viz_dir.mkdir(parents=True, exist_ok=True)
    basename = csv_path.stem.replace("_decomposedResult", "")
    prefix = viz_dir / basename

    cmd = [
        args.python,
        scripts["visualize"],
        "--mono-csv", csv_path,
        "--out-prefix", prefix,
        "--mode", args.viz_mode,
    ]
    if hor_csv is not None and hor_csv.is_file():
        cmd += ["--hor-csv", hor_csv]
    if bed is not None and bed.is_file():
        cmd += ["--centro-bed", bed]

    rc, _, _ = run_cmd(cmd)
    if rc == 0:
        log(f"  [OK] Visualization → {viz_dir}")
    else:
        log(f"  [WARN] Visualization failed (rc={rc})")


def stage5_refine(bed: Path, fasta: Path, csv_path: Path, scripts: dict, args, out_dir: Path):
    """Refine top-1 centromere region and re-annotate with anno-sat-asm."""
    refine_dir = out_dir / "05_cenSatArray_refined"
    refine_dir.mkdir(parents=True, exist_ok=True)
    basename = csv_path.stem.replace("_decomposedResult", "")

    # Read top-1 line from BED
    top1_bed = refine_dir / f".{basename}_top1.bed"
    with open(bed) as fh, open(top1_bed, "w") as out:
        first = fh.readline()
        if first:
            out.write(first)

    if top1_bed.stat().st_size == 0:
        log("  [SKIP] BED is empty, no region to refine")
        top1_bed.unlink(missing_ok=True)
        return

    # Extract FASTA region
    cmd = [
        args.python, scripts["extractBed"],
        "--bed", top1_bed,
        "--fasta", fasta,
        "--out-prefix", refine_dir / basename,
    ]
    rc, _, _ = run_cmd(cmd)
    top1_bed.unlink(missing_ok=True)

    if rc != 0:
        log(f"  [WARN] extractBedRegions failed (rc={rc})")
        return

    # Find extracted FASTA(s)
    for region_fa in sorted(refine_dir.glob(f"{basename}_*.fa")):
        region_name = region_fa.stem
        region_out = refine_dir / region_name
        region_out.mkdir(parents=True, exist_ok=True)

        log(f"    [refine] Re-annotating {region_name} with anno-sat-asm ...")
        cmd = [
            args.centroanno_bin,
            "-o", region_out,
            "-x", "anno-sat-asm",
            "-t", str(args.threads),
            "-M", str(args.max_hor_len),
            str(region_fa),
        ]
        rc, _, _ = run_cmd(cmd)
        if rc != 0:
            log(f"    [WARN] anno-sat-asm failed for {region_name} (rc={rc})")
            continue

        refined_csvs = list(region_out.rglob("*_decomposedResult.csv"))
        if not refined_csvs:
            log(f"    [WARN] No refined CSV found for {region_name}")
            continue
        refined_csv = refined_csvs[0]

        # Re-detect HORs
        cmd = [
            args.python, scripts["detectHORs"],
            "-i", refined_csv,
            "-o", region_out / region_name,
            "--min-max-copies", str(args.min_max_copies),
        ]
        rc, _, _ = run_cmd(cmd)

        refined_hor = region_out / f"{region_name}_allHORs.csv"
        if rc != 0 or not refined_hor.is_file():
            log(f"    [WARN] HOR detection failed for refined {region_name}")
            continue

        # Composite visualization
        cmd = [
            args.python, scripts["visualize"],
            "--mono-csv", refined_csv,
            "--hor-csv", refined_hor,
            "--fasta", region_fa,
            "--out-prefix", region_out / region_name,
            "--mode", args.viz_mode,
            "--composite",
        ]
        rc, _, _ = run_cmd(cmd)
        if rc == 0:
            log(f"    [OK] Composite viz for {region_name}")
        else:
            log(f"    [WARN] Composite viz failed for {region_name}")


# ---------------------------------------------------------------------------
# Per-sample worker
# ---------------------------------------------------------------------------

def process_one(args_dict, csv_path_str, fasta_path_str):
    """Worker function for parallel execution."""
    args = argparse.Namespace(**args_dict)
    csv_path = Path(csv_path_str)
    fasta = Path(fasta_path_str) if fasta_path_str else None
    out_dir = csv_path.parent
    basename = csv_path.stem.replace("_decomposedResult", "")

    log(f"{'='*70}")
    log(f"Processing: {basename}")
    log(f"  CSV : {csv_path}")
    log(f"  FASTA: {fasta}")
    log(f"{'='*70}")

    scripts = resolve_scripts(Path(args.script_dir))

    hor_csv = None
    bed = None

    if not args.skip_hors:
        log("  Stage 2/5: HOR detection")
        hor_csv = stage2_detect_hors(csv_path, scripts, args, out_dir)

    if not args.skip_centrofinder and fasta is not None and fasta.is_file():
        log("  Stage 3/5: CenSatArray detection")
        bed = stage3_centro_finder(csv_path, fasta, scripts, args, out_dir)

    if not args.skip_visualize:
        log("  Stage 4/5: Visualization")
        stage4_visualize(csv_path, hor_csv, bed, scripts, args, out_dir)

    if not args.skip_refine and bed is not None and bed.is_file():
        log("  Stage 5/5: Refinement")
        stage5_refine(bed, fasta, csv_path, scripts, args, out_dir)

    log(f"[DONE] {basename}")
    return basename


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Continue centroAnno pipeline from Stage-1 results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Batch mode — auto-discover CSVs and match FASTAs
  python3 continue_pipeline.py -d benchmark_CHM13/ -F Data/CHM13_2.0_Human/ -t 16

  # Single sample
  python3 continue_pipeline.py -c chr1_decomposedResult.csv -f chr1.fasta -t 16

  # Skip refinement
  python3 continue_pipeline.py -d out/ -F fastas/ --skip-refine
        """,
    )

    # Input / output
    parser.add_argument("-d", "--dir", type=Path, default=None,
                        help="Root directory containing chromosome subdirectories with CSVs")
    parser.add_argument("-c", "--csv", type=Path, default=None,
                        help="Single *_decomposedResult.csv file to process")
    parser.add_argument("-f", "--fasta", type=Path, default=None,
                        help="Source FASTA for the single CSV (use with -c)")
    parser.add_argument("-F", "--fasta-dir", type=Path, default=None,
                        help="Directory with source FASTAs (auto-matched by basename)")

    # Paths
    parser.add_argument("--script-dir", type=Path,
                        default=Path(__file__).parent,
                        help="Directory containing Python helper scripts [default: same dir as this script]")
    parser.add_argument("--centroanno-bin", type=Path,
                        default=Path(__file__).parent.parent / "centroAnno",
                        help="Path to centroAnno binary [default: ../centroAnno]")
    parser.add_argument("--python", default="python3",
                        help="Python interpreter [default: python3]")

    # Parallelism
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Threads for anno-sat-asm refinement [default: 8]")
    parser.add_argument("-j", "--jobs", type=int, default=1,
                        help="Number of parallel samples to process [default: 1]")

    # Visualization
    parser.add_argument("--mode", dest="viz_mode", default="png", choices=["png", "svg"],
                        help="Visualization format [default: png]")

    # HOR parameters
    parser.add_argument("--max-hor-len", type=int, default=50)
    parser.add_argument("--min-hor-copies", type=int, default=2)
    parser.add_argument("--min-hor-identity", type=float, default=0.90)
    parser.add_argument("--min-max-copies", type=int, default=3)

    # Centromere parameters
    parser.add_argument("--density-thr", type=float, default=0.8)
    parser.add_argument("--gap-factor", type=int, default=3)
    parser.add_argument("--min-size", type=int, default=50000)
    parser.add_argument("--max-size", type=int, default=10000000)
    parser.add_argument("--min-mono-len", type=int, default=80)
    parser.add_argument("--mono-gap-bp", type=int, default=100000)
    parser.add_argument("--identity-boundary-thr", type=float, default=0.8)
    parser.add_argument("--no-offset", action="store_true",
                        help="Disable auto coordinate offset for centroFinder")

    # Skip stages
    parser.add_argument("--skip-hors", action="store_true",
                        help="Skip Stage 2 (HOR detection)")
    parser.add_argument("--skip-centrofinder", action="store_true",
                        help="Skip Stage 3 (CenSatArray detection)")
    parser.add_argument("--skip-visualize", action="store_true",
                        help="Skip Stage 4 (basic visualization)")
    parser.add_argument("--skip-refine", action="store_true",
                        help="Skip Stage 5 (centromere refinement)")

    args = parser.parse_args()

    # Validate inputs
    if args.csv is not None:
        # Single mode
        if not args.csv.is_file():
            log(f"[ERROR] CSV not found: {args.csv}")
            sys.exit(1)
        if args.fasta is None and args.fasta_dir is None:
            log("[ERROR] -f or -F required when using -c (single CSV mode)")
            sys.exit(1)
        csv_files = [args.csv]
    elif args.dir is not None:
        if not args.dir.is_dir():
            log(f"[ERROR] Directory not found: {args.dir}")
            sys.exit(1)
        csv_files = find_csv_files(args.dir)
        if not csv_files:
            log(f"[ERROR] No *_decomposedResult.csv found under {args.dir}")
            sys.exit(1)
        log(f"Discovered {len(csv_files)} decomposition result(s) under {args.dir}")
    else:
        log("[ERROR] Either -d (directory) or -c (single CSV) is required.")
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Resolve binary & scripts
    args.centroanno_bin = resolve_binary(args.centroanno_bin)
    resolve_scripts(args.script_dir)

    # Build per-sample (csv, fasta) pairs
    tasks = []
    for csv_path in csv_files:
        fasta = infer_fasta(csv_path, args.fasta, args.fasta_dir)
        if fasta is None or not fasta.is_file():
            log(f"[SKIP] Cannot infer FASTA for {csv_path}")
            continue
        tasks.append((csv_path, fasta))

    if not tasks:
        log("[ERROR] No valid (CSV, FASTA) pairs found.")
        sys.exit(1)

    log(f"Ready to process {len(tasks)} sample(s) with {args.jobs} parallel job(s)")

    # Convert Namespace to dict for pickling in multiprocessing
    args_dict = vars(args)

    if args.jobs <= 1:
        for csv_path, fasta in tasks:
            process_one(args_dict, str(csv_path), str(fasta))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(process_one, args_dict, str(csv_path), str(fasta)): csv_path
                for csv_path, fasta in tasks
            }
            for future in as_completed(futures):
                csv_path = futures[future]
                try:
                    future.result()
                except Exception as e:
                    log(f"[ERROR] Exception processing {csv_path}: {e}")

    log("=" * 70)
    log("PIPELINE COMPLETE")
    log("=" * 70)


if __name__ == "__main__":
    main()
