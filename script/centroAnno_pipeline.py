#!/usr/bin/env python3
"""
centroAnno Pipeline
===================
A unified pipeline that chains:
  1. centroAnno    – de novo monomer decomposition
  2. detectHORs.py – HOR detection from monomer blocks
  3. centroFinder.py – CenSatArray region detection
  4. visualizeMonomerArray.py – monomer-array + HOR + centromere visualization

All results are organized under a single output folder.
"""

import argparse
import os
import sys
import subprocess
import glob
import time
from pathlib import Path


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [pipeline] {msg}", flush=True)


def run_cmd(cmd_list, step_name, check=True):
    """Run a command and log its execution."""
    cmd_str = " ".join(str(c) for c in cmd_list)
    log(f"[{step_name}] Running: {cmd_str}")
    start = time.time()
    result = subprocess.run(cmd_list, capture_output=False, text=True)
    elapsed = time.time() - start
    if result.returncode != 0:
        log(f"[{step_name}] FAILED (exit={result.returncode}, {elapsed:.1f}s)")
        if check:
            sys.exit(result.returncode)
    else:
        log(f"[{step_name}] Completed in {elapsed:.1f}s")
    return result


def resolve_executable(name, hint_dir):
    """Find an executable; prefer hint_dir, then PATH."""
    if hint_dir:
        candidate = os.path.join(hint_dir, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    # Search PATH
    for path in os.environ.get("PATH", "").split(os.pathsep):
        candidate = os.path.join(path, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    # Try CWD
    if os.path.isfile(name) and os.access(name, os.X_OK):
        return os.path.abspath(name)
    return name  # let shell resolve


def main():
    parser = argparse.ArgumentParser(
        description="centroAnno unified pipeline: monomer decomposition → HOR detection → CenSatArray detection → visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline on a centromeric sequence
  python script/centroAnno_pipeline.py -i example/cen21.fa -o output/

  # Full pipeline on a whole chromosome (anno-asm mode)
  python script/centroAnno_pipeline.py -i chr1.fa -o output/ -x anno-asm -t 16

  # Skip HOR detection or CenSatArray detection
  python script/centroAnno_pipeline.py -i chr1.fa -o output/ --skip-hors
  python script/centroAnno_pipeline.py -i chr1.fa -o output/ --skip-centrofinder

  # Use custom centroAnno binary and script directory
  python script/centroAnno_pipeline.py -i chr1.fa -o output/ \
      --centroanno-path ./centroAnno --script-dir ./script
""")
    # ------------------------------------------------------------------
    # I/O
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA/FASTQ file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory (will be created)")
    parser.add_argument("--centroanno-path", default="",
                        help="Path to centroAnno executable [default: search PATH/CWD]")
    parser.add_argument("--script-dir", default="",
                        help="Directory containing detectHORs.py, centroFinder.py, visualizeMonomerArray.py [default: same dir as this script]")

    # ------------------------------------------------------------------
    # centroAnno parameters
    parser.add_argument("-x", "--mode", default="anno-sat-asm",
                        choices=["anno-sat-asm", "anno-asm", "anno-read"],
                        help="centroAnno annotation mode [default: anno-sat-asm]")
    parser.add_argument("-m", "--mono-template", default="",
                        help="Monomer template FASTA (optional)")
    parser.add_argument("-k", "--kmer", type=int, default=11,
                        help="k-mer size [default: 11]")
    parser.add_argument("-f", "--fps-cutoff", type=float, default=0.6,
                        help="FPS cutoff [default: 0.6]")
    parser.add_argument("-r", "--rep-cutoff", type=float, default=0.2,
                        help="Repeat ratio cutoff [default: 0.2]")
    parser.add_argument("-w", "--window", type=int, default=500000,
                        help="Window size for template inference [default: 500000]")
    parser.add_argument("-c", "--hpc", type=int, default=1,
                        help="Homopolymer compression (1=yes, 0=no) [default: 1]")
    parser.add_argument("-e", "--epsilon", type=float, default=0.95,
                        help="DBSCAN identity cutoff [default: 0.95]")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads [default: 8]")
    parser.add_argument("-M", "--max-hor-len", type=int, default=50,
                        help="Maximum monomers in a HOR [default: 50]")
    parser.add_argument("-L", "--length-cutoff", type=int, default=5000,
                        help="Minimum sequence length to annotate [default: 5000]")
    parser.add_argument("-A", "--max-region-len", type=int, default=6000000,
                        help="Maximum region length for genome annotation [default: 6000000]")
    parser.add_argument("-N", "--min-region-len", type=int, default=100,
                        help="Minimum region length for genome annotation [default: 100]")
    parser.add_argument("-F", "--genome-cutoff", type=float, default=0.8,
                        help="Identity cutoff for genome annotation [default: 0.8]")

    # ------------------------------------------------------------------
    # detectHORs.py parameters
    parser.add_argument("--skip-hors", action="store_true",
                        help="Skip HOR detection")
    parser.add_argument("--min-hor-copies", type=int, default=2,
                        help="Minimum HOR copies per occurrence [default: 2]")
    parser.add_argument("--min-hor-identity", type=float, default=0.90,
                        help="Minimum mean identity of HOR region [default: 0.90]")
    parser.add_argument("--min-max-copies", type=int, default=2,
                        help="Minimum max copies for a (region, pattern) to be reported [default: 2]")

    # ------------------------------------------------------------------
    # centroFinder.py parameters
    parser.add_argument("--skip-centrofinder", action="store_true",
                        help="Skip CenSatArray region detection")
    parser.add_argument("--density-thr", type=float, default=0.8,
                        help="TR density threshold [default: 0.8]")
    parser.add_argument("--gap-factor", type=int, default=3,
                        help="Gap merge tolerance [default: 3]")
    parser.add_argument("--min-size", type=int, default=50000,
                        help="Minimum candidate size (bp) [default: 50000]")
    parser.add_argument("--max-size", type=int, default=10000000,
                        help="Maximum candidate size (bp) [default: 10000000]")
    parser.add_argument("--no-offset", action="store_true",
                        help="Disable auto-detected coordinate offset")

    # ------------------------------------------------------------------
    # Visualization
    parser.add_argument("--skip-visualize", action="store_true",
                        help="Skip visualization")
    parser.add_argument("--viz-format", default="png",
                        choices=["png", "svg", "pdf", "jpg"],
                        help="Visualization output format [default: png]")

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Resolve paths
    input_fa = os.path.abspath(args.input)
    outdir = os.path.abspath(args.output)
    os.makedirs(outdir, exist_ok=True)

    centroanno_exe = resolve_executable(
        args.centroanno_path if args.centroanno_path else "centroAnno",
        os.path.dirname(os.path.abspath(__file__))
    )
    if not os.path.isfile(centroanno_exe) or not os.access(centroanno_exe, os.X_OK):
        # Try project root relative to script dir
        script_dir = os.path.dirname(os.path.abspath(__file__))
        proj_root = os.path.dirname(script_dir)
        alt = os.path.join(proj_root, "centroAnno")
        if os.path.isfile(alt) and os.access(alt, os.X_OK):
            centroanno_exe = alt

    if args.script_dir:
        script_dir = os.path.abspath(args.script_dir)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))

    detect_hors_py = os.path.join(script_dir, "detectHORs.py")
    centro_finder_py = os.path.join(script_dir, "centroFinder.py")
    viz_py = os.path.join(script_dir, "visualizeMonomerArray.py")

    log("=" * 70)
    log("centroAnno Pipeline")
    log(f"  Input        : {input_fa}")
    log(f"  Output dir   : {outdir}")
    log(f"  Mode         : {args.mode}")
    log(f"  centroAnno   : {centroanno_exe}")
    log(f"  Script dir   : {script_dir}")
    log("=" * 70)

    # ------------------------------------------------------------------
    # Stage 1: centroAnno
    log("STAGE 1/4: centroAnno – monomer decomposition")
    ca_outdir = os.path.join(outdir, "01_centroAnno")
    os.makedirs(ca_outdir, exist_ok=True)

    ca_cmd = [
        centroanno_exe,
        "-o", ca_outdir,
        "-x", args.mode,
        "-k", str(args.kmer),
        "-f", str(args.fps_cutoff),
        "-r", str(args.rep_cutoff),
        "-w", str(args.window),
        "-c", str(args.hpc),
        "-e", str(args.epsilon),
        "-t", str(args.threads),
        "-M", str(args.max_hor_len),
        "-L", str(args.length_cutoff),
        "-A", str(args.max_region_len),
        "-N", str(args.min_region_len),
        "-F", str(args.genome_cutoff),
        input_fa
    ]
    if args.mono_template:
        ca_cmd.insert(-1, "-m")
        ca_cmd.insert(-1, os.path.abspath(args.mono_template))

    run_cmd(ca_cmd, "centroAnno")

    # Collect centroAnno outputs
    decomp_csvs = sorted(glob.glob(os.path.join(ca_outdir, "*_decomposedResult.csv")))
    if not decomp_csvs:
        log("ERROR: No *_decomposedResult.csv found in centroAnno output!")
        sys.exit(1)
    log(f"  Found {len(decomp_csvs)} decomposition result file(s)")

    # ------------------------------------------------------------------
    # Stage 2: HOR detection
    if not args.skip_hors:
        log("STAGE 2/4: detectHORs.py – HOR detection")
        hor_outdir = os.path.join(outdir, "02_HORs")
        os.makedirs(hor_outdir, exist_ok=True)
        for csv_path in decomp_csvs:
            basename = os.path.basename(csv_path).replace("_decomposedResult.csv", "")
            prefix = os.path.join(hor_outdir, basename)
            hor_cmd = [
                sys.executable, detect_hors_py,
                "-i", csv_path,
                "-o", prefix,
                "-m", str(args.max_hor_len),
                "-c", str(args.min_hor_copies),
                "--min-identity", str(args.min_hor_identity),
                "--min-max-copies", str(args.min_max_copies),
            ]
            run_cmd(hor_cmd, f"detectHORs({basename})", check=False)
    else:
        log("STAGE 2/4: HOR detection SKIPPED")

    # ------------------------------------------------------------------
    # Stage 3: CenSatArray detection
    if not args.skip_centrofinder:
        log("STAGE 3/4: centroFinder.py – CenSatArray detection")
        cf_outdir = os.path.join(outdir, "03_cenSatArray")
        os.makedirs(cf_outdir, exist_ok=True)
        for csv_path in decomp_csvs:
            basename = os.path.basename(csv_path).replace("_decomposedResult.csv", "")
            prefix = os.path.join(cf_outdir, basename)
            cf_cmd = [
                sys.executable, centro_finder_py,
                "--mono-csv", csv_path,
                "--fasta", input_fa,
                "--out-prefix", prefix,
                "--window", "50000",
                "--step", "10000",
                "--density-thr", str(args.density_thr),
                "--gap-factor", str(args.gap_factor),
                "--min-size", str(args.min_size),
                "--max-size", str(args.max_size),
            ]
            # For anno-sat-asm / anno-read the input FASTA is already the
            # subregion/reads; no coordinate offset should be applied.
            if args.no_offset or args.mode in ("anno-sat-asm", "anno-read"):
                cf_cmd.append("--no-offset")
            run_cmd(cf_cmd, f"centroFinder({basename})", check=False)
    else:
        log("STAGE 3/4: CenSatArray detection SKIPPED")

    # ------------------------------------------------------------------
    # Stage 4: Visualization
    if not args.skip_visualize:
        log("STAGE 4/4: visualizeMonomerArray.py – visualization")
        viz_outdir = os.path.join(outdir, "04_visualization")
        os.makedirs(viz_outdir, exist_ok=True)
        for csv_path in decomp_csvs:
            basename = os.path.basename(csv_path).replace("_decomposedResult.csv", "")
            viz_prefix = os.path.join(viz_outdir, basename)
            viz_cmd = [
                sys.executable, viz_py,
                "--mono-csv", csv_path,
                "--out-prefix", viz_prefix,
                "--format", args.viz_format,
            ]
            # Add HOR file if available
            hor_csv = os.path.join(outdir, "02_HORs", f"{basename}_allHORs.csv")
            if os.path.isfile(hor_csv) and not args.skip_hors:
                viz_cmd.extend(["--hor-csv", hor_csv])
            # Add centromere BED if available
            cf_bed = os.path.join(outdir, "03_cenSatArray", f"{basename}_censatarray_candidates.bed")
            if os.path.isfile(cf_bed) and not args.skip_centrofinder:
                viz_cmd.extend(["--centro-bed", cf_bed])
            run_cmd(viz_cmd, f"visualize({basename})", check=False)
    else:
        log("STAGE 4/4: Visualization SKIPPED")

    # ------------------------------------------------------------------
    # Summary report
    log("=" * 70)
    log("PIPELINE COMPLETE")
    log(f"  Results written to: {outdir}")
    log("")
    log("Output structure:")
    log(f"  {outdir}/")
    log(f"  ├── 01_centroAnno/           – Monomer decomposition & templates")
    log(f"  │     └── *_decomposedResult.csv")
    log(f"  │     └── *_monomerTemplates.fa")
    if not args.skip_hors:
        log(f"  ├── 02_HORs/                 – HOR patterns & statistics")
        log(f"  │     └── *_allHORs.csv")
        log(f"  │     └── *_allHORStats.txt")
        log(f"  │     └── *_allHORPatterns.txt")
        log(f"  │     └── *_allHORPatterns_compressed.txt")
    if not args.skip_centrofinder:
        log(f"  ├── 03_cenSatArray/           – CenSatArray detections & signal plots")
        log(f"  │     └── *_censatarray_candidates.bed")
        log(f"  │     └── *_censatarray_report.txt")
        log(f"  │     └── *_censatarray_signals.png")
    if not args.skip_visualize:
        log(f"  ├── 04_visualization/        – Monomer-array & HOR figures")
        log(f"  │     └── *_monomerArray.{args.viz_format}")
    log("=" * 70)


if __name__ == "__main__":
    main()
