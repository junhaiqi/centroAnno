#!/usr/bin/env bash
# =============================================================================
# run_centroAnno.sh — Fully automated centroAnno pipeline
# =============================================================================
# Chains: centroAnno → detectHORs.py → centroFinder.py → visualizeMonomerArray.py
# All results are organized into a single output folder.
#
# Usage:
#   ./run_centroAnno.sh -i <input.fa> -o <output_dir> [options]
#
# Examples:
#   ./run_centroAnno.sh -i example/cen21.fa -o results/cen21
#   ./run_centroAnno.sh -i chr1.fa -o results/chr1 -t 16 -A 2000000
#   ./run_centroAnno.sh -i reads.fa -o results/reads -x anno-read
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Determine script / binary paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CENTROANNO_BIN="${SCRIPT_DIR}/centroAnno"
SCRIPT_PY_DIR="${SCRIPT_DIR}/script"

if [[ ! -x "${CENTROANNO_BIN}" ]]; then
    echo "[ERROR] centroAnno binary not found or not executable: ${CENTROANNO_BIN}" >&2
    echo "        Please compile first:  cd ${SCRIPT_DIR} && make -j8" >&2
    exit 1
fi

for py in detectHORs.py centroFinder.py visualizeMonomerArray.py; do
    if [[ ! -f "${SCRIPT_PY_DIR}/${py}" ]]; then
        echo "[ERROR] Missing Python script: ${SCRIPT_PY_DIR}/${py}" >&2
        exit 1
    fi
done

PYTHON="${PYTHON:-python3}"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INPUT=""
OUTDIR=""
MODE="anno-asm"
MONO_TEMPLATE=""
THREADS="8"
KMER="10"
FPS_CUTOFF="0.6"
REP_CUTOFF="0.3"
WINDOW="500000"
HPC="1"
EPSILON="0.95"
MAX_HOR_LEN="50"
LENGTH_CUTOFF="5000"
MAX_REGION_LEN="1000000"
MIN_REGION_LEN="100"
GENOME_CUTOFF="0.8"

MIN_HOR_COPIES="2"
MIN_HOR_IDENTITY="0.90"
MIN_MAX_COPIES="2"

DENSITY_THR="0.8"
GAP_FACTOR="3"
MIN_SIZE="50000"
MAX_SIZE="30000000"
MIN_MONO_LEN="80"
MONO_GAP_BP="100000"
IDENTITY_BOUNDARY_THR="0.8"

SKIP_HORS=0
SKIP_CENTRO=0
SKIP_VIZ=0
NO_OFFSET=0
VIZ_MODE="png"
FORCE_CENTRO=0
REFINE_CENTRO=1
SKIP_CENTROANNO=0

# ---------------------------------------------------------------------------
# Helper: timestamped log
# ---------------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [pipeline] $*"
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input.fa> -o <output_dir> [options]

Required:
  -i FILE    Input FASTA/FASTQ
  -o DIR     Output directory (will be created)

Main options:
  -x MODE    Analysis mode: anno-asm | anno-sat-asm | anno-read [default: anno-asm]
  -m FILE    Monomer template FASTA (optional)
  -t INT     Threads [default: 8]
  -k INT     k-mer size [default: 10]
  -f FLOAT   FPS cutoff [default: 0.6]
  -r FLOAT   Repeat ratio cutoff [default: 0.3]
  -w INT     Window size for template inference [default: 500000]
  -c INT     Homopolymer compression (1=yes, 0=no) [default: 1]
  -e FLOAT   DBSCAN identity cutoff [default: 0.95]
  -M INT     Max monomers in a HOR [default: 50]
  -L INT     Minimum sequence length to annotate [default: 5000]
  -A INT     Max region length (genome mode) [default: 1000000]
  -N INT     Min region length (genome mode) [default: 100]
  -F FLOAT   Genome annotation identity cutoff [default: 0.8]

HOR options:
  --min-hor-copies INT      [default: 2]
  --min-hor-identity FLOAT  [default: 0.90]
  --min-max-copies INT      [default: 2]

Centromere options:
  --density-thr FLOAT       [default: 0.8]
  --gap-factor INT          [default: 3]
  --min-size INT            [default: 50000]
  --max-size INT            [default: 30000000]
  --min-mono-len INT        [default: 80]  Minimum mean monomer length for centromere scoring
  --mono-gap-bp INT         [default: 100000] Max tolerated low-monomer gap inside candidate (bp)
  --identity-boundary-thr FLOAT [default: 0.8] Trim boundary where mean monomer identity < threshold
  --no-offset               Disable auto coordinate offset for centroFinder

Visualization:
  --mode FORMAT             png | svg [default: png]

Skip / force stages:
  --skip-centroanno         Skip centroAnno (reuse existing 01_centroAnno/ outputs)
  --skip-hors               Skip HOR detection
  --skip-centrofinder       Skip centromere prediction
  --centrofinder            Force centromere prediction even in non-anno-asm modes
  --skip-visualize          Skip visualization
  --refine-centromere       After anno-asm, re-annotate predicted centromere
                            regions with anno-sat-asm for higher resolution
                            [default: enabled]
  --no-refine-centromere    Skip centromere region refinement

Examples:
  $(basename "$0") -i example/cen21.fa -o results/cen21
  $(basename "$0") -i chr1.fa -o results/chr1 -t 16 -A 2000000
  $(basename "$0") -i reads.fa -o results/reads -x anno-read
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)        INPUT="$2"; shift 2 ;;
        -o|--output)       OUTDIR="$2"; shift 2 ;;
        -x|--mode)         MODE="$2"; shift 2 ;;
        -m|--mono-template) MONO_TEMPLATE="$2"; shift 2 ;;
        -t|--threads)      THREADS="$2"; shift 2 ;;
        -k|--kmer)         KMER="$2"; shift 2 ;;
        -f|--fps-cutoff)   FPS_CUTOFF="$2"; shift 2 ;;
        -r|--rep-cutoff)   REP_CUTOFF="$2"; shift 2 ;;
        -w|--window)       WINDOW="$2"; shift 2 ;;
        -c|--hpc)          HPC="$2"; shift 2 ;;
        -e|--epsilon)      EPSILON="$2"; shift 2 ;;
        -M|--max-hor-len)  MAX_HOR_LEN="$2"; shift 2 ;;
        -L|--length-cutoff) LENGTH_CUTOFF="$2"; shift 2 ;;
        -A|--max-region-len) MAX_REGION_LEN="$2"; shift 2 ;;
        -N|--min-region-len) MIN_REGION_LEN="$2"; shift 2 ;;
        -F|--genome-cutoff) GENOME_CUTOFF="$2"; shift 2 ;;
        --min-hor-copies)   MIN_HOR_COPIES="$2"; shift 2 ;;
        --min-hor-identity) MIN_HOR_IDENTITY="$2"; shift 2 ;;
        --min-max-copies)   MIN_MAX_COPIES="$2"; shift 2 ;;
        --density-thr)      DENSITY_THR="$2"; shift 2 ;;
        --gap-factor)       GAP_FACTOR="$2"; shift 2 ;;
        --min-size)         MIN_SIZE="$2"; shift 2 ;;
        --max-size)         MAX_SIZE="$2"; shift 2 ;;
        --min-mono-len)     MIN_MONO_LEN="$2"; shift 2 ;;
        --mono-gap-bp)      MONO_GAP_BP="$2"; shift 2 ;;
        --identity-boundary-thr) IDENTITY_BOUNDARY_THR="$2"; shift 2 ;;
        --mode)             VIZ_MODE="$2"; shift 2 ;;
        --skip-hors)        SKIP_HORS=1; shift ;;
        --skip-centrofinder) SKIP_CENTRO=1; shift ;;
        --centrofinder)      FORCE_CENTRO=1; shift ;;
        --skip-visualize)   SKIP_VIZ=1; shift ;;
        --refine-centromere) REFINE_CENTRO=1; shift ;;
        --no-refine-centromere) REFINE_CENTRO=0; shift ;;
        --skip-centroanno)  SKIP_CENTROANNO=1; shift ;;
        --no-offset)        NO_OFFSET=1; shift ;;
        -h|--help)          usage; exit 0 ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${INPUT}" || -z "${OUTDIR}" ]]; then
    echo "[ERROR] -i (input) and -o (output) are required." >&2
    usage >&2
    exit 1
fi

# Default: only anno-asm runs centromere prediction.
# anno-sat-asm / anno-read skip it unless --centrofinder is given.
if [[ "${MODE}" != "anno-asm" && ${FORCE_CENTRO} -eq 0 ]]; then
    SKIP_CENTRO=1
fi

INPUT_ABS="$(cd "$(dirname "${INPUT}")" && pwd)/$(basename "${INPUT}")"
if [[ ! -f "${INPUT_ABS}" ]]; then
    echo "[ERROR] Input file not found: ${INPUT_ABS}" >&2
    exit 1
fi
OUTDIR_ABS="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"

# ---------------------------------------------------------------------------
# Stage 0: Print banner
# ---------------------------------------------------------------------------
log "========================================================================"
log "centroAnno Automated Pipeline"
log "  Input      : ${INPUT_ABS}"
log "  Output dir : ${OUTDIR_ABS}"
log "  Mode       : ${MODE}"
log "========================================================================"

# ---------------------------------------------------------------------------
# Stage 1: centroAnno — monomer decomposition
# ---------------------------------------------------------------------------
CA_OUT="${OUTDIR_ABS}/01_centroAnno"
mkdir -p "${CA_OUT}"

if [[ ${SKIP_CENTROANNO} -eq 1 ]]; then
    log "STAGE 1/4: centroAnno SKIPPED (using existing outputs in ${CA_OUT})"
else
    log "STAGE 1/4: centroAnno — monomer decomposition"
    CA_ARGS=(
        -o "${CA_OUT}"
        -x "${MODE}"
        -k "${KMER}"
        -f "${FPS_CUTOFF}"
        -r "${REP_CUTOFF}"
        -w "${WINDOW}"
        -c "${HPC}"
        -e "${EPSILON}"
        -t "${THREADS}"
        -M "${MAX_HOR_LEN}"
        -L "${LENGTH_CUTOFF}"
        -A "${MAX_REGION_LEN}"
        -N "${MIN_REGION_LEN}"
        -F "${GENOME_CUTOFF}"
    )

    if [[ -n "${MONO_TEMPLATE}" ]]; then
        CA_ARGS+=(-m "$(cd "$(dirname "${MONO_TEMPLATE}")" && pwd)/$(basename "${MONO_TEMPLATE}")")
    fi

    CA_ARGS+=("${INPUT_ABS}")
    /usr/bin/time -v --output="${CA_OUT}/centroAnno_time.log" "${CENTROANNO_BIN}" "${CA_ARGS[@]}"
fi

# Collect decomposition CSV files
mapfile -t CSV_FILES < <(find "${CA_OUT}" -maxdepth 1 -name '*_decomposedResult.csv' | sort)
if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
    log "ERROR: No *_decomposedResult.csv found in ${CA_OUT}"
    exit 1
fi
log "  Found ${#CSV_FILES[@]} decomposition result file(s)"

# ---------------------------------------------------------------------------
# Stage 2: detectHORs.py — HOR detection
# ---------------------------------------------------------------------------
if [[ ${SKIP_HORS} -eq 0 ]]; then
    log "STAGE 2/4: detectHORs.py — HOR detection"
    HOR_OUT="${OUTDIR_ABS}/02_HORs"
    mkdir -p "${HOR_OUT}"

    for csv in "${CSV_FILES[@]}"; do
        BASENAME=$(basename "${csv}" _decomposedResult.csv)
        PREFIX="${HOR_OUT}/${BASENAME}"
        log "  [detectHORs] Processing: ${BASENAME}"
        ${PYTHON} "${SCRIPT_PY_DIR}/detectHORs.py" \
            -i "${csv}" \
            -o "${PREFIX}" \
            -m "${MAX_HOR_LEN}" \
            -c "${MIN_HOR_COPIES}" \
            --min-identity "${MIN_HOR_IDENTITY}" \
            --min-max-copies "${MIN_MAX_COPIES}" || true
    done
else
    log "STAGE 2/4: HOR detection SKIPPED"
fi

# ---------------------------------------------------------------------------
# Stage 3: centroFinder.py — CenSatArray detection
# ---------------------------------------------------------------------------
if [[ ${SKIP_CENTRO} -eq 0 ]]; then
    log "STAGE 3/4: centroFinder.py — CenSatArray detection"
    CF_OUT="${OUTDIR_ABS}/03_cenSatArray"
    mkdir -p "${CF_OUT}"

    # For anno-sat-asm / anno-read the input FASTA IS the subregion/reads,
    # so we should disable coordinate offset.
    # Also, if the sequence name itself looks like chr:start-end, it is already
    # a subregion FASTA and offset should be disabled.
    AUTO_NO_OFFSET=0
    if [[ "${MODE}" == "anno-sat-asm" || "${MODE}" == "anno-read" ]]; then
        AUTO_NO_OFFSET=1
    fi

    for csv in "${CSV_FILES[@]}"; do
        BASENAME=$(basename "${csv}" _decomposedResult.csv)
        PREFIX="${CF_OUT}/${BASENAME}"
        log "  [centroFinder] Processing: ${BASENAME}"

        # Detect subregion-style sequence names (e.g. CP068257.1:11699867-12031015)
        CSV_NO_OFFSET=${AUTO_NO_OFFSET}
        if [[ "${BASENAME}" == *":"*"-"* && ${CSV_NO_OFFSET} -eq 0 ]]; then
            CSV_NO_OFFSET=1
            log "    Auto-detected subregion name → --no-offset"
        fi

        CF_ARGS=(
            --mono-csv "${csv}"
            --fasta "${INPUT_ABS}"
            --out-prefix "${PREFIX}"
            --window 50000
            --step 10000
            --density-thr "${DENSITY_THR}"
            --gap-factor "${GAP_FACTOR}"
            --min-size "${MIN_SIZE}"
            --max-size "${MAX_SIZE}"
            --min-mono-len "${MIN_MONO_LEN}"
            --mono-gap-bp "${MONO_GAP_BP}"
            --identity-boundary-thr "${IDENTITY_BOUNDARY_THR}"
        )
        if [[ ${NO_OFFSET} -eq 1 || ${CSV_NO_OFFSET} -eq 1 ]]; then
            CF_ARGS+=(--no-offset)
        fi
        ${PYTHON} "${SCRIPT_PY_DIR}/centroFinder.py" "${CF_ARGS[@]}" || true
    done
else
    log "STAGE 3/4: CenSatArray detection SKIPPED"
fi

# ---------------------------------------------------------------------------
# Stage 4: visualizeMonomerArray.py — visualization
# ---------------------------------------------------------------------------
if [[ ${SKIP_VIZ} -eq 0 ]]; then
    log "STAGE 4/4: visualizeMonomerArray.py — visualization"
    VIZ_OUT="${OUTDIR_ABS}/04_visualization"
    mkdir -p "${VIZ_OUT}"

    for csv in "${CSV_FILES[@]}"; do
        BASENAME=$(basename "${csv}" _decomposedResult.csv)
        PREFIX="${VIZ_OUT}/${BASENAME}"
        log "  [visualize] Processing: ${BASENAME}"
        VIZ_ARGS=(
            --mono-csv "${csv}"
            --out-prefix "${PREFIX}"
            --mode "${VIZ_MODE}"
        )
        HOR_CSV="${OUTDIR_ABS}/02_HORs/${BASENAME}_allHORs.csv"
        if [[ -f "${HOR_CSV}" && ${SKIP_HORS} -eq 0 ]]; then
            VIZ_ARGS+=(--hor-csv "${HOR_CSV}")
        fi
        CF_BED="${OUTDIR_ABS}/03_cenSatArray/${BASENAME}_censatarray_candidates.bed"
        if [[ -f "${CF_BED}" && ${SKIP_CENTRO} -eq 0 ]]; then
            VIZ_ARGS+=(--centro-bed "${CF_BED}")
        fi
        ${PYTHON} "${SCRIPT_PY_DIR}/visualizeMonomerArray.py" "${VIZ_ARGS[@]}" || true
    done
else
    log "STAGE 4/4: Visualization SKIPPED"
fi

# ---------------------------------------------------------------------------
# Stage 5: Refine centromere regions (optional, anno-asm only)
# ---------------------------------------------------------------------------
if [[ ${REFINE_CENTRO} -eq 1 ]]; then
    if [[ "${MODE}" != "anno-asm" ]]; then
        log "STAGE 5/5: Refinement SKIPPED (only available in anno-asm mode)"
    else
        log "STAGE 5/5: Refining CenSatArray regions with anno-sat-asm"
        REFINE_OUT="${OUTDIR_ABS}/05_cenSatArray_refined"
        mkdir -p "${REFINE_OUT}"

        EXTRACT_PY="${SCRIPT_PY_DIR}/extractBedRegions.py"
        if [[ ! -f "${EXTRACT_PY}" ]]; then
            log "  [refine] ERROR: extractBedRegions.py not found at ${EXTRACT_PY}"
            exit 1
        fi

        for bed in "${OUTDIR_ABS}"/03_cenSatArray/*_censatarray_candidates.bed; do
            [[ -f "${bed}" ]] || continue
            BASENAME=$(basename "${bed}" _censatarray_candidates.bed)
            log "  [refine] Processing BED: ${BASENAME}"

            # Extract only the top-1 (best) centromere region
            top1_bed="${REFINE_OUT}/.${BASENAME}_top1.bed"
            head -1 "${bed}" > "${top1_bed}"

            ${PYTHON} "${EXTRACT_PY}" \
                --bed "${top1_bed}" \
                --fasta "${INPUT_ABS}" \
                --out-prefix "${REFINE_OUT}/${BASENAME}" \
                || true

            rm -f "${top1_bed}"

            # Re-annotate the extracted region with anno-sat-asm
            for region_fa in "${REFINE_OUT}/${BASENAME}"_*.fa; do
                [[ -f "${region_fa}" ]] || continue
                region_name=$(basename "${region_fa}" .fa)
                region_out="${REFINE_OUT}/${region_name}"
                mkdir -p "${region_out}"

                log "    [refine] Annotating ${region_name} with anno-sat-asm ..."
                "${CENTROANNO_BIN}" -o "${region_out}" \
                    -x anno-sat-asm \
                    -t "${THREADS}" \
                    -k "${KMER}" \
                    -f "${FPS_CUTOFF}" \
                    -r "${REP_CUTOFF}" \
                    -w "${WINDOW}" \
                    -c "${HPC}" \
                    -e "${EPSILON}" \
                    -M "${MAX_HOR_LEN}" \
                    -L "${LENGTH_CUTOFF}" \
                    "${region_fa}" || true

                # Find the refined decomposition CSV
                refined_csv=$(find "${region_out}" -maxdepth 1 -name '*_decomposedResult.csv' | head -1)
                if [[ -f "${refined_csv}" ]]; then
                    log "    [refine] Detecting HORs for ${region_name} ..."
                    "${PYTHON}" "${SCRIPT_PY_DIR}/detectHORs.py" \
                        -i "${refined_csv}" \
                        -o "${region_out}/${region_name}" \
                        --min-max-copies "${MIN_MAX_COPIES}" \
                        || true

                    refined_hor="${region_out}/${region_name}_allHORs.csv"
                    if [[ -f "${refined_hor}" ]]; then
                        log "    [refine] Visualizing composite for ${region_name} ..."
                        "${PYTHON}" "${SCRIPT_PY_DIR}/visualizeMonomerArray.py" \
                            --mono-csv "${refined_csv}" \
                            --hor-csv "${refined_hor}" \
                            --fasta "${region_fa}" \
                            --out-prefix "${region_out}/${region_name}" \
                            --mode "${VIZ_MODE}" \
                            --composite \
                            || true
                    fi
                fi
            done
        done
    fi
else
    log "STAGE 5/5: Refinement SKIPPED"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log "========================================================================"
log "PIPELINE COMPLETE"
log "  Results: ${OUTDIR_ABS}"
log ""
log "Output structure:"
log "  ${OUTDIR_ABS}/"
log "  ├── 01_centroAnno/           – Monomer decomposition & templates"
log "  ├── 02_HORs/                 – HOR patterns & statistics"
log "  ├── 03_cenSatArray/           – CenSatArray candidates & signal plots"
log "  ├── 04_visualization/        – Figures & maps"
log "  └── 05_cenSatArray_refined/   – Refined CenSatArray annotations (optional)"
log "========================================================================"
