#!/usr/bin/env bash
# =============================================================================
# continue_pipeline.sh — Resume centroAnno Stages 2–5 from existing CSVs.
#
# Purpose:
#   You already ran Stage 1 (centroAnno monomer decomposition) and have one or
#   more *_decomposedResult.csv files. This script picks up from there and runs
#   Stages 2–5 (HOR detection, CenSatArray detection, visualization, optional
#   refinement) for every CSV it finds.
#
# FASTA dependency by stage:
#   - Stage 2 (HOR detection)    : only needs the CSV.
#   - Stage 4 (basic visualization): only needs the CSV (plus Stage 2 output).
#   - Stage 3 (CenSatArray)      : requires the original FASTA used for Stage 1.
#   - Stage 5 (refinement)       : requires the original FASTA.
#
#   Therefore -F / -f is optional IF you skip Stage 3 and Stage 5.
#
# Usage:
#   ./continue_pipeline.sh -d <bench_dir> [-F <fasta_dir> | -f <fasta>] [options]
#
# Examples:
#   # HOR detection + visualization only (no FASTA needed)
#   ./continue_pipeline.sh -d results/ --skip-centrofinder --skip-refine
#
#   # Full Stages 2–5 (needs FASTA)
#   ./continue_pipeline.sh -d results/ -F Data/ -j 8 -t 16
#
#   # Process only chr1-22,chrX,chrY
#   ./continue_pipeline.sh -d results/ -F Data/ --target-regex '^chr([0-9]+|X|Y)$'
# =============================================================================

# If sourced in worker mode, only define functions and exit.
if [[ "${_CP_WORKER_MODE:-}" == "1" ]]; then
    _define_functions_only=1
fi

set -uo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PY_DIR="${SCRIPT_DIR}/script"
CENTROANNO_BIN="${SCRIPT_DIR}/centroAnno"
PYTHON="${PYTHON:-python3}"

for py in detectHORs.py centroFinder.py visualizeMonomerArray.py extractBedRegions.py; do
    if [[ ! -f "${SCRIPT_PY_DIR}/${py}" ]]; then
        echo "[ERROR] Missing Python script: ${SCRIPT_PY_DIR}/${py}" >&2
        exit 1
    fi
done

if [[ ! -x "${CENTROANNO_BIN}" ]]; then
    echo "[ERROR] centroAnno binary not found: ${CENTROANNO_BIN}" >&2
    echo "        Please compile first: cd ${SCRIPT_DIR} && make -j8" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BENCH_DIR=""
FASTA_DIR=""
FASTA_FILE=""
THREADS=16
JOBS=1
VIZ_MODE="png"

MAX_HOR_LEN=50
MIN_HOR_COPIES=2
MIN_HOR_IDENTITY=0.90
MIN_MAX_COPIES=2

DENSITY_THR=0.8
GAP_FACTOR=3
MIN_SIZE=50000
MAX_SIZE=10000000
MIN_MONO_LEN=80
MONO_GAP_BP=100000
IDENTITY_BOUNDARY_THR=0.8
NO_OFFSET=0

SKIP_HORS=0
SKIP_CENTRO=0
SKIP_VIZ=0
SKIP_REFINE=0
STRICT=0
DRY_RUN=0
TARGET_REGEX=".*"
MAX_DEPTH=2

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [continue-pipeline] $*"
}

warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [continue-pipeline] [WARN] $*" >&2
}

error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [continue-pipeline] [ERROR] $*" >&2
}

usage() {
    cat <<EOF
Usage: $(basename "$0") -d <bench_dir> [options]

Required:
  -d DIR          Benchmark output directory containing *_decomposedResult.csv files

FASTA (only needed for Stage 3/Stage 5):
  -F DIR          Directory with source FASTA files (*.fasta, *.fa, *.fna)
  -f FILE         Single FASTA file to use for ALL CSVs (e.g. a multi-sequence file)

Options:
  -j INT          Parallel samples [default: 1]
  -t INT          Threads for anno-sat-asm refinement [default: 16]
  --mode          Visualization format: png | svg [default: png]
  --max-depth INT Max depth to search for CSVs under BENCH_DIR [default: 2]
  --target-regex  Only process CSVs whose basename matches this regex [default: .*]
  --min-max-copies INT        [default: 3]
  --min-mono-len INT          [default: 80]
  --mono-gap-bp INT           [default: 100000]
  --identity-boundary-thr FLOAT [default: 0.8]
  --no-offset                 Disable auto coordinate offset for centroFinder
  --skip-hors                 Skip Stage 2
  --skip-centrofinder         Skip Stage 3
  --skip-visualize            Skip Stage 4
  --skip-refine               Skip Stage 5
  --strict                    Exit immediately if any stage command fails
  --dry-run                   Show what would be run, then exit without running
  -h, --help                  Show this help

Examples:
  # Just HOR detection + visualization, no FASTA needed
  $(basename "$0") -d results/ --skip-centrofinder --skip-refine

  # Full Stages 2-5 using FASTA directory
  $(basename "$0") -d results/ -F Data/ -j 8 -t 16

  # Single FASTA file (multi-sequence) for all CSVs
  $(basename "$0") -d results/ -f genome.fa -j 4

  # Only canonical human chromosomes
  $(basename "$0") -d results/ -F Data/ --target-regex '^chr([0-9]+|X|Y)$'
EOF
}

# ---------------------------------------------------------------------------
# Parse args (only in main mode)
# ---------------------------------------------------------------------------
if [[ "${_define_functions_only:-}" != "1" ]]; then
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -d|--bench-dir)      BENCH_DIR="$2"; shift 2 ;;
            -F|--fasta-dir)      FASTA_DIR="$2"; shift 2 ;;
            -f|--fasta)          FASTA_FILE="$2"; shift 2 ;;
            -j|--jobs)           JOBS="$2"; shift 2 ;;
            -t|--threads)        THREADS="$2"; shift 2 ;;
            --mode)              VIZ_MODE="$2"; shift 2 ;;
            --max-depth)         MAX_DEPTH="$2"; shift 2 ;;
            --target-regex)      TARGET_REGEX="$2"; shift 2 ;;
            --min-max-copies)    MIN_MAX_COPIES="$2"; shift 2 ;;
            --min-mono-len)      MIN_MONO_LEN="$2"; shift 2 ;;
            --mono-gap-bp)       MONO_GAP_BP="$2"; shift 2 ;;
            --identity-boundary-thr) IDENTITY_BOUNDARY_THR="$2"; shift 2 ;;
            --no-offset)         NO_OFFSET=1; shift ;;
            --skip-hors)         SKIP_HORS=1; shift ;;
            --skip-centrofinder) SKIP_CENTRO=1; shift ;;
            --skip-visualize)    SKIP_VIZ=1; shift ;;
            --skip-refine)       SKIP_REFINE=1; shift ;;
            --strict)            STRICT=1; shift ;;
            --dry-run)           DRY_RUN=1; shift ;;
            -h|--help)           usage; exit 0 ;;
            *)
                error "Unknown option: $1"
                usage >&2
                exit 1
                ;;
        esac
    done

    if [[ -z "${BENCH_DIR}" ]]; then
        error "-d is required."
        usage >&2
        exit 1
    fi

    if [[ ! -d "${BENCH_DIR}" ]]; then
        error "Benchmark directory not found: ${BENCH_DIR}"
        exit 1
    fi
    BENCH_DIR="$(cd "${BENCH_DIR}" && pwd)"

    if [[ -n "${FASTA_FILE}" && -n "${FASTA_DIR}" ]]; then
        error "Use either -f (single FASTA) or -F (FASTA directory), not both."
        exit 1
    fi

    if [[ -n "${FASTA_FILE}" ]]; then
        if [[ ! -f "${FASTA_FILE}" ]]; then
            error "FASTA file not found: ${FASTA_FILE}"
            exit 1
        fi
        FASTA_FILE="$(cd "$(dirname "${FASTA_FILE}")" && pwd)/$(basename "${FASTA_FILE}")"
    fi

    if [[ -n "${FASTA_DIR}" ]]; then
        if [[ ! -d "${FASTA_DIR}" ]]; then
            error "FASTA directory not found: ${FASTA_DIR}"
            exit 1
        fi
        FASTA_DIR="$(cd "${FASTA_DIR}" && pwd)"
    fi

    NEEDS_FASTA=0
    [[ ${SKIP_CENTRO} -eq 0 ]] && NEEDS_FASTA=1
    [[ ${SKIP_REFINE} -eq 0 ]] && NEEDS_FASTA=1

    if [[ ${NEEDS_FASTA} -eq 1 && -z "${FASTA_FILE}" && -z "${FASTA_DIR}" ]]; then
        warn "Stage 3 (CenSatArray) and/or Stage 5 (refinement) are enabled but no FASTA was provided."
        warn "Those stages will be skipped. Provide -F or -f to enable them."
        SKIP_CENTRO=1
        SKIP_REFINE=1
    fi
fi

# ---------------------------------------------------------------------------
# Resolve output directory
# ---------------------------------------------------------------------------
resolve_out_dir() {
    local csv_path="$1"
    local parent
    parent=$(basename "$(dirname "${csv_path}")")
    if [[ "${parent}" == "01_centroAnno" ]]; then
        dirname "$(dirname "${csv_path}")"
    else
        dirname "${csv_path}"
    fi
}

# ---------------------------------------------------------------------------
# Infer FASTA for a given CSV
# ---------------------------------------------------------------------------
infer_fasta() {
    local csv_path="$1"
    local base parent
    base=$(basename "${csv_path}" _decomposedResult.csv)
    parent=$(basename "$(dirname "${csv_path}")")

    if [[ -n "${FASTA_FILE}" ]]; then
        echo "${FASTA_FILE}"
        return 0
    fi

    if [[ -z "${FASTA_DIR}" ]]; then
        echo ""
        return 1
    fi

    local candidates=()
    # Search recursively under FASTA_DIR
    while IFS= read -r -d '' fa; do
        candidates+=("${fa}")
    done < <(find "${FASTA_DIR}" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) -print0)

    # Priority 1: exact basename match
    for ext in .fasta .fa .fna; do
        local c
        for c in "${candidates[@]}"; do
            if [[ "$(basename "${c}")" == "${base}${ext}" || "$(basename "${c}")" == "${base}${ext}.gz" ]]; then
                echo "${c}"
                return 0
            fi
        done
    done

    # Priority 2: parent directory name match
    for ext in .fasta .fa .fna; do
        local c
        for c in "${candidates[@]}"; do
            if [[ "$(basename "${c}")" == "${parent}${ext}" || "$(basename "${c}")" == "${parent}${ext}.gz" ]]; then
                echo "${c}"
                return 0
            fi
        done
    done

    # Priority 3: basename appears as prefix in FASTA name
    local c
    for c in "${candidates[@]}"; do
        local bn
        bn=$(basename "${c}")
        bn="${bn%.gz}"
        bn="${bn%.fasta}"
        bn="${bn%.fa}"
        bn="${bn%.fna}"
        if [[ "${base}" == *"${bn}"* || "${bn}" == *"${base}"* ]]; then
            echo "${c}"
            return 0
        fi
    done

    echo ""
    return 1
}

# ---------------------------------------------------------------------------
# Per-sample worker
# ---------------------------------------------------------------------------
process_sample() {
    local csv_path="$1"
    local fasta="$2"
    local out_dir
    out_dir=$(resolve_out_dir "${csv_path}")
    local basename
    basename=$(basename "${csv_path}" _decomposedResult.csv)

    log "================================================================================"
    log "Processing: ${basename}"
    log "  CSV : ${csv_path}"
    if [[ -n "${fasta}" ]]; then
        log "  FASTA: ${fasta}"
    else
        log "  FASTA: (not provided)"
    fi
    log "================================================================================"

    local hor_csv=""
    local bed=""

    # --- Stage 2: HOR detection ---
    if [[ ${SKIP_HORS} -eq 0 ]]; then
        log "  [${basename}] Stage 2/5: HOR detection"
        local hor_dir="${out_dir}/02_HORs"
        mkdir -p "${hor_dir}"
        set +e
        ${PYTHON} "${SCRIPT_PY_DIR}/detectHORs.py" \
            -i "${csv_path}" \
            -o "${hor_dir}/${basename}" \
            -m "${MAX_HOR_LEN}" \
            -c "${MIN_HOR_COPIES}" \
            --min-identity "${MIN_HOR_IDENTITY}" \
            --min-max-copies "${MIN_MAX_COPIES}" \
            > "${hor_dir}/${basename}_detectHORs.log" 2>&1
        local rc=$?
        set -e
        if [[ ${rc} -ne 0 ]]; then
            if [[ ${STRICT} -eq 1 ]]; then
                error "${basename}: Stage 2 (HOR detection) failed (exit ${rc}). Aborting."
                return 1
            fi
            warn "${basename}: Stage 2 (HOR detection) failed (exit ${rc}). Continuing."
        fi
        hor_csv="${hor_dir}/${basename}_allHORs.csv"
        [[ -f "${hor_csv}" ]] || hor_csv=""
    fi

    # --- Stage 3: CenSatArray detection ---
    if [[ ${SKIP_CENTRO} -eq 0 && -n "${fasta}" && -f "${fasta}" ]]; then
        log "  [${basename}] Stage 3/5: CenSatArray detection"
        local cf_dir="${out_dir}/03_cenSatArray"
        mkdir -p "${cf_dir}"

        # Auto-detect subregion-style sequence names (e.g. chrY:10500000-10950000)
        local csv_no_offset=${NO_OFFSET}
        if [[ "${basename}" =~ ^[^:]+:[0-9]+-[0-9]+$ && ${csv_no_offset} -eq 0 ]]; then
            csv_no_offset=1
            log "    [${basename}] Auto-detected subregion name → --no-offset"
        fi

        local cf_args=(
            --mono-csv "${csv_path}"
            --fasta "${fasta}"
            --out-prefix "${cf_dir}/${basename}"
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
        [[ ${csv_no_offset} -eq 1 ]] && cf_args+=(--no-offset)
        set +e
        ${PYTHON} "${SCRIPT_PY_DIR}/centroFinder.py" "${cf_args[@]}" \
            > "${cf_dir}/${basename}_centroFinder.log" 2>&1
        local rc=$?
        set -e
        if [[ ${rc} -ne 0 ]]; then
            if [[ ${STRICT} -eq 1 ]]; then
                error "${basename}: Stage 3 (CenSatArray detection) failed (exit ${rc}). Aborting."
                return 1
            fi
            warn "${basename}: Stage 3 (CenSatArray detection) failed (exit ${rc}). Continuing."
        fi
        bed="${cf_dir}/${basename}_censatarray_candidates.bed"
        [[ -f "${bed}" ]] || bed=""
    fi

    # --- Stage 4: Visualization (basic) ---
    if [[ ${SKIP_VIZ} -eq 0 ]]; then
        log "  [${basename}] Stage 4/5: Visualization"
        local viz_dir="${out_dir}/04_visualization"
        mkdir -p "${viz_dir}"
        local viz_args=(
            --mono-csv "${csv_path}"
            --out-prefix "${viz_dir}/${basename}"
            --mode "${VIZ_MODE}"
        )
        [[ -n "${hor_csv}" && -f "${hor_csv}" ]] && viz_args+=(--hor-csv "${hor_csv}")
        [[ -n "${bed}" && -f "${bed}" ]] && viz_args+=(--centro-bed "${bed}")
        set +e
        ${PYTHON} "${SCRIPT_PY_DIR}/visualizeMonomerArray.py" "${viz_args[@]}" \
            > "${viz_dir}/${basename}_visualize.log" 2>&1
        local rc=$?
        set -e
        if [[ ${rc} -ne 0 ]]; then
            if [[ ${STRICT} -eq 1 ]]; then
                error "${basename}: Stage 4 (visualization) failed (exit ${rc}). Aborting."
                return 1
            fi
            warn "${basename}: Stage 4 (visualization) failed (exit ${rc}). Continuing."
        fi
    fi

    # --- Stage 5: Refinement ---
    if [[ ${SKIP_REFINE} -eq 0 && -n "${bed}" && -f "${bed}" && -n "${fasta}" && -f "${fasta}" ]]; then
        log "  [${basename}] Stage 5/5: Refinement"
        local refine_dir="${out_dir}/05_cenSatArray_refined"
        mkdir -p "${refine_dir}"

        local top1_bed="${refine_dir}/.${basename}_top1.bed"
        head -1 "${bed}" > "${top1_bed}"

        if [[ ! -s "${top1_bed}" ]]; then
            log "    [${basename}] BED empty, skipping refinement"
            rm -f "${top1_bed}"
        else
            set +e
            ${PYTHON} "${SCRIPT_PY_DIR}/extractBedRegions.py" \
                --bed "${top1_bed}" \
                --fasta "${fasta}" \
                --out-prefix "${refine_dir}/${basename}" \
                > "${refine_dir}/${basename}_extract.log" 2>&1
            local rc=$?
            set -e
            rm -f "${top1_bed}"
            if [[ ${rc} -ne 0 ]]; then
                if [[ ${STRICT} -eq 1 ]]; then
                    error "${basename}: Stage 5 (extract BED region) failed (exit ${rc}). Aborting."
                    return 1
                fi
                warn "${basename}: Stage 5 (extract BED region) failed (exit ${rc}). Continuing."
            fi

            for region_fa in "${refine_dir}/${basename}"_*.fa; do
                [[ -f "${region_fa}" ]] || continue
                local region_name
                region_name=$(basename "${region_fa}" .fa)
                local region_out="${refine_dir}/${region_name}"
                mkdir -p "${region_out}"

                log "    [${basename}] Refining ${region_name} with anno-sat-asm ..."
                set +e
                "${CENTROANNO_BIN}" -o "${region_out}" \
                    -x anno-sat-asm \
                    -t "${THREADS}" \
                    -M "${MAX_HOR_LEN}" \
                    "${region_fa}" \
                    > "${region_out}/${region_name}_centroAnno.log" 2>&1
                rc=$?
                set -e
                if [[ ${rc} -ne 0 ]]; then
                    if [[ ${STRICT} -eq 1 ]]; then
                        error "${basename}: Stage 5 (anno-sat-asm on ${region_name}) failed (exit ${rc}). Aborting."
                        return 1
                    fi
                    warn "${basename}: Stage 5 (anno-sat-asm on ${region_name}) failed (exit ${rc}). Continuing."
                    continue
                fi

                local refined_csv
                refined_csv=$(find "${region_out}" -maxdepth 1 -name '*_decomposedResult.csv' | head -1)
                if [[ -z "${refined_csv}" || ! -f "${refined_csv}" ]]; then
                    warn "    [${basename}] No refined CSV for ${region_name}"
                    continue
                fi

                log "    [${basename}] Detecting HORs for ${region_name} ..."
                set +e
                ${PYTHON} "${SCRIPT_PY_DIR}/detectHORs.py" \
                    -i "${refined_csv}" \
                    -o "${region_out}/${region_name}" \
                    --min-max-copies "${MIN_MAX_COPIES}" \
                    > "${region_out}/${region_name}_detectHORs.log" 2>&1
                rc=$?
                set -e
                if [[ ${rc} -ne 0 ]]; then
                    warn "    [${basename}] HOR detection failed for refined ${region_name} (exit ${rc})"
                    continue
                fi

                local refined_hor="${region_out}/${region_name}_allHORs.csv"
                if [[ ! -f "${refined_hor}" ]]; then
                    warn "    [${basename}] No HORs for refined ${region_name}"
                    continue
                fi

                log "    [${basename}] Composite viz for ${region_name} ..."
                set +e
                ${PYTHON} "${SCRIPT_PY_DIR}/visualizeMonomerArray.py" \
                    --mono-csv "${refined_csv}" \
                    --hor-csv "${refined_hor}" \
                    --fasta "${region_fa}" \
                    --out-prefix "${region_out}/${region_name}" \
                    --mode "${VIZ_MODE}" \
                    --composite \
                    > "${region_out}/${region_name}_visualize.log" 2>&1
                rc=$?
                set -e
                if [[ ${rc} -ne 0 ]]; then
                    warn "    [${basename}] Composite visualization failed for refined ${region_name} (exit ${rc})"
                fi
            done
        fi
    fi

    log "[DONE] ${basename}"
}

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------
if [[ "${_define_functions_only:-}" == "1" ]]; then
    return 0
fi

log "========================================================================"
log "centroAnno Continue Pipeline"
log "  Benchmark dir  : ${BENCH_DIR}"
if [[ -n "${FASTA_FILE}" ]]; then
    log "  FASTA file    : ${FASTA_FILE}"
elif [[ -n "${FASTA_DIR}" ]]; then
    log "  FASTA dir     : ${FASTA_DIR}"
else
    log "  FASTA         : not provided (Stage 3 & 5 will be skipped)"
fi
log "  Parallel jobs : ${JOBS}"
log "  Threads/refine: ${THREADS}"
log "  Target regex  : ${TARGET_REGEX}"
log "  Skip stages   : HORS=${SKIP_HORS} CENTRO=${SKIP_CENTRO} VIZ=${SKIP_VIZ} REFINE=${SKIP_REFINE}"
log "========================================================================"

# Gather all CSV files
mapfile -t CSV_FILES < <(find "${BENCH_DIR}" -maxdepth "${MAX_DEPTH}" -name '*_decomposedResult.csv' | sort)
if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
    error "No *_decomposedResult.csv found under ${BENCH_DIR} (max depth ${MAX_DEPTH})"
    exit 1
fi
log "Found ${#CSV_FILES[@]} decomposition result(s)"

# Build task list
SKIPPED_NOTARGET=0
SKIPPED_NOFASTA=0
declare -a TASKS=()
for csv in "${CSV_FILES[@]}"; do
    base=$(basename "${csv}" _decomposedResult.csv)
    if [[ ! "${base}" =~ ${TARGET_REGEX} ]]; then
        ((SKIPPED_NOTARGET++))
        log "[SKIP] ${base}: does not match target regex"
        continue
    fi

    fa=$(infer_fasta "${csv}")
    if [[ ${NEEDS_FASTA} -eq 1 && ( -z "${fa}" || ! -f "${fa}" ) ]]; then
        ((SKIPPED_NOFASTA++))
        warn "[SKIP] ${base}: cannot infer FASTA (Stage 3/5 enabled)"
        continue
    fi

    TASKS+=("${csv}|${fa}")
done

log "Ready to process ${#TASKS[@]} sample(s) with ${JOBS} parallel job(s)"
log "Skipped by regex: ${SKIPPED_NOTARGET}, skipped by missing FASTA: ${SKIPPED_NOFASTA}"

if [[ ${DRY_RUN} -eq 1 ]]; then
    log "DRY RUN — would process the following samples:"
    for task in "${TASKS[@]}"; do
        csv="${task%%|*}"
        fa="${task##*|}"
        log "  CSV: ${csv}  FASTA: ${fa:-<none>}"
    done
    exit 0
fi

if [[ ${#TASKS[@]} -eq 0 ]]; then
    error "No valid tasks after filtering."
    exit 1
fi

# Run workers
if [[ ${JOBS} -le 1 ]]; then
    for task in "${TASKS[@]}"; do
        csv="${task%%|*}"
        fa="${task##*|}"
        process_sample "${csv}" "${fa}"
    done
else
    # Export all variables needed by worker subshells
    export SCRIPT_DIR SCRIPT_PY_DIR CENTROANNO_BIN PYTHON THREADS VIZ_MODE STRICT
    export MAX_HOR_LEN MIN_HOR_COPIES MIN_HOR_IDENTITY MIN_MAX_COPIES
    export DENSITY_THR GAP_FACTOR MIN_SIZE MAX_SIZE MIN_MONO_LEN MONO_GAP_BP
    export IDENTITY_BOUNDARY_THR NO_OFFSET SKIP_HORS SKIP_CENTRO SKIP_VIZ SKIP_REFINE

    printf '%s\n' "${TASKS[@]}" | xargs -P "${JOBS}" -I {} bash -c '
        task="$1"
        csv="${task%%|*}"
        fa="${task##*|}"
        _CP_WORKER_MODE=1 source "${SCRIPT_DIR}/continue_pipeline.sh"
        process_sample "${csv}" "${fa}"
    ' _ {}
fi

log "========================================================================"
log "PIPELINE COMPLETE"
log "  Processed      : ${#TASKS[@]} sample(s)"
log "  Skipped (regex): ${SKIPPED_NOTARGET}"
log "  Skipped (FASTA): ${SKIPPED_NOFASTA}"
log "========================================================================"
