#!/usr/bin/env python3
"""
centroFinder.py

CenSatArray candidate region detector based on CentroAnno outputs.
Inspired by:
  - CentroVision  : TR density + gap merge as hard constraint
  - CentIER       : k-mer entropy (inverted: low entropy = satellite)
  - RepeatOBserver: Shannon diversity extremes
  - CentriVision  : multi-signal voting

Algorithm (4-stage fusion):
  Stage 1: TR Density      — monomer coverage fraction per window
  Stage 2: Complexity      — inverted 4-mer Shannon entropy per window
  Stage 3: Conservation    — mean monomer identity per window
  Stage 4: Merge & Predict — density threshold + gap merge + entropy boundary

Coordinate handling:
  CentroAnno run on subregions produces relative coordinates in the CSV,
  with the subregion encoded in the sequence name (e.g. chrY:10500000-10950000).
  The script auto-detects this prefix and converts coordinates to absolute
  so they match the full-chromosome FASTA.

Inputs:
  --mono-csv   : CentroAnno *_decomposedResult.csv
  --fasta      : Chromosome/reference FASTA
  --out-prefix : Output path prefix

Outputs:
  {prefix}_censatarray_candidates.bed
  {prefix}_censatarray_signals.png
  {prefix}_censatarray_report.txt
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ============================================================================
# Stage 0: Parse CentroAnno output + coordinate offset
# ============================================================================

def parse_seq_prefix(seq_name):
    """Extract chromosome offset from 'chrY:10500000-10950000' style names.

    Returns (chrom, offset) where offset is added to all coordinates.
    If no range prefix, offset = 0.
    """
    seq_name = seq_name.strip()
    m = re.match(r'^([^:]+):(\d+)-(\d+)$', seq_name)
    if m:
        chrom = m.group(1)
        offset = int(m.group(2))
        return chrom, offset
    return seq_name, 0


_HEADER_KEYWORDS = {
    'monomer', 'mononer', 'sequence', 'start', 'end',
    'identity', 'position', 'length', 'name', 'estimated', 'strand',
}


def _looks_like_header(row):
    """Return True if a CSV row looks like a header line."""
    return any(
        any(kw in cell.lower() for kw in _HEADER_KEYWORDS)
        for cell in row
    )


def parse_mono_csv(csv_path):
    """Parse decomposedResult.csv (anno-sat with header OR anno-asm without header).

    Returns (records, auto_offset) where:
      - records: list of monomer dicts
      - auto_offset: offset extracted from sequence name (e.g. chrY:10500000-10950000)
                     For anno-asm (e.g. chr1) offset is 0 because coords are absolute.

    Formats supported:
      1. anno-sat (header present, 6 cols):
         sequence name, monomer name, start position, end position, estimated identity, length
      2. anno-asm (NO header, 7 cols):
         sequence name, monomer name, strand, start position, end position, estimated identity, length
      3. Old anno-sat (header present, 6 cols, 'mononer name'):
    """
    records = []
    auto_offset = 0
    with open(csv_path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        first_row = next(reader)

        if _looks_like_header(first_row):
            # -------------------------- anno-sat with header --------------------------
            fieldnames = [f.strip() for f in first_row]

            name_key = 'monomer name' if 'monomer name' in fieldnames else 'mononer name'
            start_key = 'start position' if 'start position' in fieldnames else 'start_pos'
            end_key = 'end position' if 'end position' in fieldnames else 'end_pos'
            seq_key = 'sequence name' if 'sequence name' in fieldnames else 'sequence_name'
            has_strand = 'strand' in fieldnames or 'Strand' in fieldnames

            fh.seek(0)
            dict_reader = csv.DictReader(fh)
            first_row_seen = True
            for row in dict_reader:
                raw_name = row[name_key].strip()
                strand = row.get('strand', row.get('Strand', '+')).strip() if has_strand else '+'
                if raw_name.startswith('Pos:') and '_' in raw_name:
                    idx = raw_name.index('_')
                    raw_name = raw_name[idx + 1:]

                # For old anno-sat without explicit strand column, infer from trailing '
                if not has_strand and raw_name.endswith("'"):
                    raw_name = raw_name[:-1]
                    strand = '-'

                mono_id = raw_name + "'" if strand == '-' else raw_name
                id_val = row.get('estimated identity', row.get('identity', '0'))
                try:
                    identity = float(id_val)
                except ValueError:
                    identity = 0.0

                seq_name = row.get(seq_key, '').strip()
                if first_row_seen and seq_name:
                    _, auto_offset = parse_seq_prefix(seq_name)
                    first_row_seen = False

                length = int(row[end_key]) - int(row[start_key]) + 1
                records.append({
                    'seq_name': seq_name,
                    'mono_id' : mono_id,
                    'strand'  : strand,
                    'start'   : int(row[start_key]),
                    'end'     : int(row[end_key]),
                    'identity': identity,
                    'length'  : length,
                })
        else:
            # -------------------------- anno-asm without header --------------------------
            # 7 fixed columns: seq_name, mono_name, strand, start, end, identity, length
            rows = [first_row]
            rows.extend(reader)

            first_row_seen = True
            for row in rows:
                if len(row) < 7:
                    continue
                seq_name = row[0].strip()
                raw_name = row[1].strip()
                strand = row[2].strip()
                if raw_name.startswith('Pos:') and '_' in raw_name:
                    idx = raw_name.index('_')
                    raw_name = raw_name[idx + 1:]

                mono_id = raw_name + "'" if strand == '-' else raw_name
                try:
                    identity = float(row[5])
                except ValueError:
                    identity = 0.0
                try:
                    length = int(row[6])
                except ValueError:
                    length = int(row[4]) - int(row[3]) + 1

                if first_row_seen:
                    _, auto_offset = parse_seq_prefix(seq_name)
                    first_row_seen = False

                records.append({
                    'seq_name': seq_name,
                    'mono_id' : mono_id,
                    'strand'  : strand,
                    'start'   : int(row[3]),
                    'end'     : int(row[4]),
                    'identity': identity,
                    'length'  : length,
                })
    return records, auto_offset


# ============================================================================
# Stage 1: TR Density (monomer coverage)
# ============================================================================

def compute_tr_density(seq_len, mono_records, window, step):
    """Return per-window monomer coverage fraction (0~1)."""
    n_win = (seq_len - window) // step + 1
    density = np.zeros(n_win, dtype=np.float32)
    cov = np.zeros(seq_len, dtype=np.int32)
    for m in mono_records:
        s = max(0, m['start'])
        e = min(seq_len, m['end'] + 1)
        if s < e:
            cov[s:e] = 1
    for i in range(n_win):
        w_s = i * step
        w_e = w_s + window
        density[i] = np.sum(cov[w_s:w_e]) / window
    return density


# ============================================================================
# Stage 2: Sequence Complexity (inverted k-mer Shannon entropy)
# ============================================================================

def kmer_entropy(seq_str, k=4):
    n = len(seq_str) - k + 1
    if n <= 0:
        return 0.0
    counts = defaultdict(int)
    for i in range(n):
        counts[seq_str[i:i+k]] += 1
    total = sum(counts.values())
    H = 0.0
    for cnt in counts.values():
        p = cnt / total
        H -= p * np.log2(p)
    return H


def compute_entropy(seq, window, step, k=4):
    """Sliding-window Shannon entropy."""
    seq = str(seq).upper()
    seq_len = len(seq)
    n_win = (seq_len - window) // step + 1
    entropy = np.zeros(n_win, dtype=np.float32)
    for i in range(n_win):
        entropy[i] = kmer_entropy(seq[i*step : i*step + window], k)
    return entropy


# ============================================================================
# Stage 3: Conservation (mean monomer identity)
# ============================================================================

def compute_conservation(seq_len, mono_records, window, step):
    """Per-window mean monomer identity."""
    n_win = (seq_len - window) // step + 1
    cons = np.zeros(n_win, dtype=np.float32)
    id_sum = np.zeros(seq_len, dtype=np.float32)
    id_cnt = np.zeros(seq_len, dtype=np.int32)
    for m in mono_records:
        s = max(0, m['start'])
        e = min(seq_len, m['end'] + 1)
        if s < e:
            id_sum[s:e] += m['identity']
            id_cnt[s:e] += 1
    for i in range(n_win):
        w_s = i * step
        w_e = w_s + window
        c = np.sum(id_cnt[w_s:w_e])
        if c > 0:
            cons[i] = np.sum(id_sum[w_s:w_e]) / c
        else:
            cons[i] = 0.0
    return cons


# ============================================================================
# Stage 3b: Monomer length (mean monomer length per window)
# ============================================================================

def compute_monomer_length(seq_len, mono_records, window, step):
    """Per-window mean monomer length in bp.

    Centromeric alpha-satellite monomers are typically ~150-200 bp,
    whereas microsatellites are <50 bp. This signal helps distinguish
    centromeric repeats from other tandem repeats.
    """
    n_win = (seq_len - window) // step + 1
    mlen = np.zeros(n_win, dtype=np.float32)
    len_sum = np.zeros(seq_len, dtype=np.float32)
    len_cnt = np.zeros(seq_len, dtype=np.int32)
    for m in mono_records:
        s = max(0, m['start'])
        e = min(seq_len, m['end'] + 1)
        if s < e:
            len_sum[s:e] += m['length']
            len_cnt[s:e] += 1
    for i in range(n_win):
        w_s = i * step
        w_e = w_s + window
        c = np.sum(len_cnt[w_s:w_e])
        if c > 0:
            mlen[i] = np.sum(len_sum[w_s:w_e]) / c
        else:
            mlen[i] = 0.0
    return mlen


def monomer_length_score(mean_length, min_mono_len=100):
    """Convert mean monomer length to a 0-1 score.

    Alpha-satellite monomers (~150-200 bp) score highest.
    Sub-100 bp repeats (e.g. microsatellites) are filtered out.
    Very large monomers (>500 bp) score moderately.
    """
    if mean_length <= 0 or mean_length < min_mono_len:
        return 0.0
    elif mean_length <= 250:
        return 1.0
    elif mean_length <= 500:
        return 1.0 - (mean_length - 250) / 250.0 * 0.5
    else:
        return 0.5


def size_score(region_bp, target=1_000_000):
    """Reward larger regions up to a target size.

    Regions smaller than the target get a linear penalty (floor 0.5),
    so that a 1 Mb centromere scores higher than a 500 Kb fragment.
    """
    if region_bp >= target:
        return 1.0
    else:
        return 0.5 + 0.5 * (region_bp / target)


# ============================================================================
# Stage 4: Merge & Predict
# ============================================================================

def gap_merge(high_mask, gap_factor=3):
    """Merge high-density windows separated by <= gap_factor windows."""
    regions = []
    in_region = False
    r_start = 0
    i = 0
    while i < len(high_mask):
        if high_mask[i] and not in_region:
            in_region = True
            r_start = i
        elif not high_mask[i] and in_region:
            # Look ahead to check gap size
            j = i + 1
            while j < len(high_mask) and not high_mask[j]:
                j += 1
            gap = j - i
            if gap > gap_factor:
                in_region = False
                regions.append((r_start, i - 1))
                i = j  # skip past the gap
                continue
            # else: gap is small, stay in region; i will be incremented
        i += 1
    if in_region:
        regions.append((r_start, len(high_mask) - 1))
    return regions


def refine_boundary(entropy, rs, re, step, window, side='left'):
    """Refine boundary using entropy gradient.

    For left boundary: find the window with steepest entropy rise
    (transition from low-entropy satellite to high-entropy non-satellite)
    within ±2 windows of the edge.

    For right boundary: find steepest entropy drop.
    """
    margin = 2
    if side == 'left':
        lo = max(0, rs - margin)
        hi = min(len(entropy) - 1, rs + margin)
        # Find max gradient (rise) in entropy = sharpest transition OUT of satellite
        best_grad = -1e9
        best_idx = rs
        for i in range(lo, hi):
            grad = entropy[i + 1] - entropy[i]
            if grad > best_grad:
                best_grad = grad
                best_idx = i
        return best_idx * step
    else:
        lo = max(0, re - margin)
        hi = min(len(entropy) - 1, re + margin)
        # Find max gradient (rise) going left-to-right = sharpest transition INTO satellite
        best_grad = -1e9
        best_idx = re
        for i in range(lo, hi):
            grad = entropy[i + 1] - entropy[i]
            if grad > best_grad:
                best_grad = grad
                best_idx = i
        # Use the right side of the steepest transition
        return min((best_idx + 1) * step + window, (len(entropy) - 1) * step + window)


def predict_centromere(density, entropy, conservation, mono_len, step, window,
                       density_thr=0.8, gap_factor=3,
                       min_size_bp=50000, max_size_bp=10_000_000,
                       min_mono_len=100, min_mono_gap_bp=100000):
    """Multi-signal fusion prediction.

    Primary filter   : TR density >= density_thr  (CentriVision hard constraint)
    Secondary filter : gap merge
    Tertiary filter  : min/max size constraints
    Refinement       : entropy gradient at boundaries
    Ranking          : score = density * conservation * length_score
                       (higher = more centromere-like)
    """
    n_win = len(density)
    high_mask = density >= density_thr
    raw_regions = gap_merge(high_mask, gap_factor)

    candidates = []
    for rs, re in raw_regions:
        if re < rs:
            continue
        # Convert window indices to base coordinates
        start_bp = rs * step
        end_bp = min(re * step + window, (n_win - 1) * step + window)

        # Filter by size
        span = end_bp - start_bp
        if span < min_size_bp:
            continue
        if span > max_size_bp:
            continue

        # Refine boundaries using entropy gradient
        left_refine = refine_boundary(entropy, rs, re, step, window, side='left')
        right_refine = refine_boundary(entropy, rs, re, step, window, side='right')
        start_bp = max(start_bp, left_refine)
        end_bp = min(end_bp, right_refine)
        span = end_bp - start_bp
        if span < min_size_bp or span > max_size_bp:
            continue

        # Split refined region by monomer length:
        # keep only sub-regions where mean monomer length >= min_mono_len.
        # Short gaps (< min_mono_gap_bp) of low monomer length are bridged to avoid
        # over-fragmentation; longer gaps split the region.
        refined_rs = start_bp // step
        refined_re = min((end_bp - 1) // step, n_win - 1)
        if refined_re < refined_rs:
            continue

        # Number of consecutive bad windows that span >= min_mono_gap_bp
        if min_mono_gap_bp > window:
            max_gap_win = max(1, (min_mono_gap_bp - window) // step + 1)
        else:
            max_gap_win = 1

        sub_rs = None
        last_good = None
        for i in range(refined_rs, refined_re + 1):
            if mono_len[i] >= min_mono_len:
                if sub_rs is None:
                    sub_rs = i
                last_good = i
            else:
                if sub_rs is not None and last_good is not None and (i - last_good) >= max_gap_win:
                    sub_re = last_good
                    sub_start = sub_rs * step
                    sub_end = min(sub_re * step + window, (n_win - 1) * step + window)
                    if sub_end - sub_start >= min_size_bp:
                        d_mean = float(np.mean(density[sub_rs:sub_re+1]))
                        c_mean = float(np.mean(conservation[sub_rs:sub_re+1]))
                        l_mean = float(np.mean(mono_len[sub_rs:sub_re+1]))
                        l_score = monomer_length_score(l_mean, min_mono_len)
                        s_score = size_score(sub_end - sub_start)
                        score = d_mean * c_mean * l_score * s_score
                        candidates.append({
                            'start': int(sub_start),
                            'end': int(sub_end),
                            'density_mean': d_mean,
                            'entropy_mean': float(np.mean(entropy[sub_rs:sub_re+1])),
                            'conservation_mean': c_mean,
                            'monomer_len_mean': round(l_mean, 1),
                            'score': round(score, 6)
                        })
                    sub_rs = None
                    last_good = None

        # trailing sub-region
        if sub_rs is not None and last_good is not None:
            sub_re = last_good
            sub_start = sub_rs * step
            sub_end = min(sub_re * step + window, (n_win - 1) * step + window)
            if sub_end - sub_start >= min_size_bp:
                d_mean = float(np.mean(density[sub_rs:sub_re+1]))
                c_mean = float(np.mean(conservation[sub_rs:sub_re+1]))
                l_mean = float(np.mean(mono_len[sub_rs:sub_re+1]))
                l_score = monomer_length_score(l_mean, min_mono_len)
                s_score = size_score(sub_end - sub_start)
                score = d_mean * c_mean * l_score * s_score
                candidates.append({
                    'start': int(sub_start),
                    'end': int(sub_end),
                    'density_mean': d_mean,
                    'entropy_mean': float(np.mean(entropy[sub_rs:sub_re+1])),
                    'conservation_mean': c_mean,
                    'monomer_len_mean': round(l_mean, 1),
                    'score': round(score, 6)
                    })

    # Sort by score descending (best centromere candidate first)
    candidates.sort(key=lambda x: x['score'], reverse=True)
    return candidates


# ============================================================================
# Visualization
# ============================================================================

def _fmt_bp(x):
    """Format base-pair coordinate with commas and Mb unit."""
    if x >= 1_000_000:
        return f"{x/1_000_000:.2f}M"
    elif x >= 1_000:
        return f"{x/1_000:.1f}k"
    else:
        return str(x)


def plot_signals(seq_len, window, step, density, entropy, conservation, mono_len,
                 candidates, out_path, chrom_name='chr'):
    """Plot all signals and predicted regions."""
    n_win = len(density)
    centers = np.arange(n_win) * step + window // 2

    fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=True,
                             gridspec_kw={'hspace': 0.06})

    # --- TR Density (top panel) ---
    ax = axes[0]
    ax.plot(centers, density, 'b-', lw=1.2)
    ax.axhline(0.8, color='r', ls='--', lw=1, label='θ=0.8')
    for c in candidates:
        ax.axvspan(c['start'], c['end'], color='green', alpha=0.2)

    # Annotate the best candidate (#1) with coordinates
    if candidates:
        best = candidates[0]
        mid = (best['start'] + best['end']) / 2
        span = best['end'] - best['start']
        label_text = (f"#1: {chrom_name}:{best['start']:,}-{best['end']:,} "
                      f"({_fmt_bp(span)}, score={best['score']:.3f})")
        # Place text near the top of the panel, centered on the candidate
        ax.text(mid, 0.92, label_text,
                ha='center', va='top', fontsize=9,
                color='darkgreen', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen',
                          edgecolor='green', alpha=0.7))

    ax.set_ylabel('TR Density')
    ax.set_title('CentroFinder: CenSatArray Candidate Detection')
    ax.legend(loc='upper right')
    ax.set_ylim(-0.05, 1.05)

    # --- Entropy ---
    ax = axes[1]
    ax.plot(centers, entropy, 'purple', lw=1)
    for c in candidates:
        ax.axvspan(c['start'], c['end'], color='green', alpha=0.2)
    ax.set_ylabel('Shannon Entropy (bits)')

    # --- Conservation ---
    ax = axes[2]
    ax.plot(centers, conservation, 'orange', lw=1)
    for c in candidates:
        ax.axvspan(c['start'], c['end'], color='green', alpha=0.2)
    ax.set_ylabel('Mean Identity')
    ax.set_ylim(-0.05, 1.05)

    # --- Monomer Length (bottom panel) ---
    ax = axes[3]
    ax.plot(centers, mono_len, 'darkgreen', lw=1)
    ax.axhline(100, color='r', ls='--', lw=0.8, label='α-sat threshold (100 bp)')
    ax.axhline(170, color='orange', ls='--', lw=0.8, label='α-sat optimal (170 bp)')
    for c in candidates:
        ax.axvspan(c['start'], c['end'], color='green', alpha=0.2)
    ax.set_ylabel('Mean Monomer Length (bp)')
    ax.set_xlabel('Genomic coordinate (bp)')
    ax.legend(loc='upper right')

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# Report
# ============================================================================

def write_report(report_path, candidates, seq_len, window, step, offset):
    with open(report_path, 'w') as fh:
        fh.write('CentroFinder Report\n')
        fh.write('=' * 60 + '\n\n')
        fh.write(f'Sequence length : {seq_len:,} bp\n')
        fh.write(f'Window size     : {window:,} bp\n')
        fh.write(f'Window step     : {step:,} bp\n')
        fh.write(f'Coordinate offset: {offset:,} bp\n')
        fh.write(f'Candidates found: {len(candidates)}\n\n')
        for i, c in enumerate(candidates, 1):
            span = c['end'] - c['start']
            marker = '  <-- BEST' if i == 1 else ''
            fh.write(f'Candidate {i}{marker}\n')
            fh.write(f'  Coordinates : {c["start"]:,} - {c["end"]:,} ({span:,} bp)\n')
            fh.write(f'  TR density  : {c["density_mean"]:.3f}\n')
            fh.write(f'  Entropy     : {c["entropy_mean"]:.3f}\n')
            fh.write(f'  Conservation: {c["conservation_mean"]:.3f}\n')
            fh.write(f'  Monomer len : {c["monomer_len_mean"]:.1f} bp\n')
            fh.write(f'  Score       : {c["score"]:.6f}\n\n')


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Predict CenSatArray regions from CentroAnno outputs')
    parser.add_argument('--mono-csv', required=True,
                        help='CentroAnno *_decomposedResult.csv')
    parser.add_argument('--fasta', required=True,
                        help='Reference FASTA (chromosome or region)')
    parser.add_argument('--out-prefix', required=True,
                        help='Output prefix (path + basename)')
    parser.add_argument('--window', type=int, default=50000,
                        help='Sliding window size in bp (default: 50000)')
    parser.add_argument('--step', type=int, default=10000,
                        help='Window step in bp (default: 10000)')
    parser.add_argument('--density-thr', type=float, default=0.8,
                        help='TR density threshold (default: 0.8)')
    parser.add_argument('--gap-factor', type=int, default=3,
                        help='Gap merge tolerance in window counts (default: 3)')
    parser.add_argument('--min-size', type=int, default=50000,
                        help='Minimum candidate size in bp (default: 50000)')
    parser.add_argument('--max-size', type=int, default=10_000_000,
                        help='Maximum candidate size in bp (default: 10000000)')
    parser.add_argument('--offset', type=int, default=None,
                        help='Manual coordinate offset in bp (overrides auto-detect)')
    parser.add_argument('--no-offset', action='store_true',
                        help='Disable auto-detected coordinate offset')
    parser.add_argument('--min-mono-len', type=int, default=80,
                        help='Minimum mean monomer length (bp) for centromere scoring '
                             '(sub-threshold windows get zero length score) [default: 100]')
    parser.add_argument('--mono-gap-bp', type=int, default=100000,
                        help='Maximum tolerated gap of low monomer-length (<min-mono-len) '
                             'inside a candidate region in bp. Gaps shorter than this are '
                             'bridged; longer gaps split the region. [default: 100000]')
    parser.add_argument('--identity-boundary-thr', type=float, default=0.8,
                        help='Trim candidate boundaries where mean monomer identity '
                             'is below this threshold. [default: 0.8]')
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load FASTA (multi-sequence support)
    # ------------------------------------------------------------------
    print(f'[centroFinder] Loading FASTA {args.fasta} ...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    print(f'[centroFinder] Loaded {len(seq_dict)} sequence(s)')

    print(f'[centroFinder] Parsing {args.mono_csv} ...')
    mono_records, auto_offset = parse_mono_csv(args.mono_csv)

    if not mono_records:
        print('[centroFinder] ERROR: No monomer records parsed from CSV.', file=sys.stderr)
        sys.exit(1)

    # Determine target sequence name from CSV records
    target_seq_name = mono_records[0].get('seq_name', '')
    if not target_seq_name:
        # Fallback: derive from CSV filename (basename before _decomposedResult.csv)
        target_seq_name = os.path.basename(args.mono_csv).replace('_decomposedResult.csv', '')

    # Handle coordinate offset
    if args.no_offset:
        offset = 0
        print('[centroFinder] Coordinate offset disabled by --no-offset')
    elif args.offset is not None:
        offset = args.offset
        print(f'[centroFinder] Manual coordinate offset: {offset:,} bp')
    else:
        offset = auto_offset
        if offset > 0:
            print(f'[centroFinder] Auto-detected coordinate offset: {offset:,} bp')

    # Apply offset to monomer coordinates
    if offset > 0:
        for m in mono_records:
            m['start'] += offset
            m['end'] += offset

    # Find matching FASTA record
    # Try exact match first, then fallback to chrom prefix extraction
    if target_seq_name in seq_dict:
        record = seq_dict[target_seq_name]
    else:
        chrom_from_prefix, _ = parse_seq_prefix(target_seq_name)
        if chrom_from_prefix in seq_dict:
            record = seq_dict[chrom_from_prefix]
        else:
            print(f"[centroFinder] ERROR: Sequence '{target_seq_name}' (or '{chrom_from_prefix}') "
                  f"not found in FASTA. Available: {list(seq_dict.keys())[:5]}...", file=sys.stderr)
            sys.exit(1)

    seq = str(record.seq).upper()
    seq_len = len(seq)
    print(f'[centroFinder] Target: {record.id} ({seq_len:,} bp)  |  '
          f'{len(mono_records)} monomer records')

    # Sanity check
    if mono_records:
        max_coord = max(m['end'] for m in mono_records)
        if max_coord > seq_len:
            print(f'[centroFinder] WARNING: max monomer coordinate ({max_coord:,}) '
                  f'exceeds sequence length ({seq_len:,}).')
            print(f'[centroFinder]          If you are using the same subregion FASTA '
                  f'that generated the CSV, use --no-offset.')
            sys.exit(1)

    # ------------------------------------------------------------------
    # Compute signals
    # ------------------------------------------------------------------
    print(f'[centroFinder] Stage 1: TR density (w={args.window}, s={args.step}) ...')
    density = compute_tr_density(seq_len, mono_records, args.window, args.step)

    print(f'[centroFinder] Stage 2: k-mer entropy ...')
    entropy = compute_entropy(seq, args.window, args.step, k=4)

    print(f'[centroFinder] Stage 3: conservation score ...')
    conservation = compute_conservation(seq_len, mono_records, args.window, args.step)

    print(f'[centroFinder] Stage 3b: monomer length ...')
    mono_len = compute_monomer_length(seq_len, mono_records, args.window, args.step)

    print(f'[centroFinder] Stage 4: fusion & prediction ...')
    candidates = predict_centromere(
        density, entropy, conservation, mono_len, args.step, args.window,
        density_thr=args.density_thr, gap_factor=args.gap_factor,
        min_size_bp=args.min_size, max_size_bp=args.max_size,
        min_mono_len=args.min_mono_len, min_mono_gap_bp=args.mono_gap_bp)
    print(f'[centroFinder] {len(candidates)} candidate region(s) found')

    # ------------------------------------------------------------------
    # Post-process: trim boundaries by block-level monomer identity + length
    # (signal-level window averages are too coarse for sharp boundaries;
    #  we evaluate small blocks of monomers from each edge inward.)
    # ------------------------------------------------------------------
    if candidates and args.identity_boundary_thr > 0:
        trimmed = []
        for c in candidates:
            region_monos = [m for m in mono_records
                            if m['start'] < c['end'] and m['end'] > c['start']]
            region_monos.sort(key=lambda m: m['start'])
            if not region_monos:
                trimmed.append(c)
                continue

            MIN_BLOCK_BP = 5000
            MIN_BLOCK_MONOS = 5

            # ---- Left trim ----
            new_start = c['start']
            i = 0
            while i < len(region_monos):
                block = []
                block_bp = 0
                j = i
                while j < len(region_monos):
                    block.append(region_monos[j])
                    block_bp += region_monos[j]['end'] - region_monos[j]['start']
                    j += 1
                    if block_bp >= MIN_BLOCK_BP and len(block) >= MIN_BLOCK_MONOS:
                        break
                if not block:
                    break
                avg_id = float(np.mean([m['identity'] for m in block]))
                avg_len = float(np.mean([m['length'] for m in block]))
                if avg_id < args.identity_boundary_thr or avg_len < args.min_mono_len:
                    # Bad block → skip and keep trimming
                    i = j
                    new_start = region_monos[min(j, len(region_monos)-1)]['start']
                else:
                    # Good block → stop left trim
                    new_start = region_monos[i]['start']
                    break

            # ---- Right trim ----
            new_end = c['end']
            i = len(region_monos) - 1
            while i >= 0:
                block = []
                block_bp = 0
                j = i
                while j >= 0:
                    block.append(region_monos[j])
                    block_bp += region_monos[j]['end'] - region_monos[j]['start']
                    j -= 1
                    if block_bp >= MIN_BLOCK_BP and len(block) >= MIN_BLOCK_MONOS:
                        break
                if not block:
                    break
                avg_id = float(np.mean([m['identity'] for m in block]))
                avg_len = float(np.mean([m['length'] for m in block]))
                if avg_id < args.identity_boundary_thr or avg_len < args.min_mono_len:
                    i = j
                    new_end = region_monos[max(j, 0)]['end']
                else:
                    new_end = region_monos[i]['end']
                    break

            if new_end - new_start < args.min_size:
                continue  # drop if too small after trimming

            c['start'] = int(new_start)
            c['end'] = int(new_end)

            # Re-calculate stats for trimmed window range
            n_win = len(density)
            trimmed_rs = c['start'] // args.step
            trimmed_re = min((c['end'] - 1) // args.step, n_win - 1)
            if trimmed_re >= trimmed_rs:
                c['density_mean'] = float(np.mean(density[trimmed_rs:trimmed_re+1]))
                c['entropy_mean'] = float(np.mean(entropy[trimmed_rs:trimmed_re+1]))
                c['conservation_mean'] = float(np.mean(conservation[trimmed_rs:trimmed_re+1]))
                c['monomer_len_mean'] = round(float(np.mean(mono_len[trimmed_rs:trimmed_re+1])), 1)
                l_score = monomer_length_score(c['monomer_len_mean'], args.min_mono_len)
                s_score = size_score(c['end'] - c['start'])
                c['score'] = round(c['density_mean'] * c['conservation_mean'] * l_score * s_score, 6)
            trimmed.append(c)

        candidates = trimmed
        candidates.sort(key=lambda x: x['score'], reverse=True)
        print(f'[centroFinder] {len(candidates)} candidate region(s) after block trimming')

    if candidates:
        best = candidates[0]
        print(f'[centroFinder] Best candidate: {best["start"]:,}-{best["end"]:,} '
              f'(score={best["score"]:.4f}, monomer_len={best["monomer_len_mean"]:.1f}bp)')

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    bed_path = args.out_prefix + '_censatarray_candidates.bed'
    with open(bed_path, 'w') as fh:
        for i, c in enumerate(candidates, 1):
            fh.write(f'{record.id}\t{c["start"]}\t{c["end"]}\t'
                     f'centroFinder_{i}\t{c["density_mean"]:.3f}\t+\n')
    print(f'[centroFinder] BED written: {bed_path}')

    plot_path = args.out_prefix + '_censatarray_signals.png'
    plot_signals(seq_len, args.window, args.step, density, entropy, conservation,
                 mono_len, candidates, plot_path, chrom_name=record.id)
    print(f'[centroFinder] Plot written: {plot_path}')

    # Safety sort: ensure report is always ordered by score descending
    candidates.sort(key=lambda x: x['score'], reverse=True)
    report_path = args.out_prefix + '_censatarray_report.txt'
    write_report(report_path, candidates, seq_len, args.window, args.step, offset)
    print(f'[centroFinder] Report written: {report_path}')

    print('[centroFinder] Done!')


if __name__ == '__main__':
    main()
