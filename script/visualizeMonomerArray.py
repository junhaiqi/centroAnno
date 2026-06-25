#!/usr/bin/env python3
"""
visualizeMonomerArray.py
========================
Comprehensive visualization for centroAnno outputs.

Generates:
  1. Tandem-repeat region map (cautils style)     (*_repeatRegions.{mode})
  2. HOR region map (cautils style)               (*_horRegions.{mode})
  3. HiCAT-style HOR tracks + diamond heatmap     (*_composite_*.png)
"""

import argparse
import csv
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from Bio import SeqIO


# =============================================================================
# HiCAT-style color map & helpers
# =============================================================================

HICAT_COLOR_MAP = {
    0: "#ffffff",
    1: "#4b3991", 2: "#2974af", 3: "#4a9da8", 4: "#57b894",
    5: "#7dd873", 6: "#c1f676", 7: "#ffff92", 8: "#fdda59",
    9: "#fb9e3f", 10: "#ee5624", 11: "#c9272e", 12: "#6a0023",
}

HICAT_THRESHOLDS = [
    (0.00, 0.50), (0.50, 0.65), (0.65, 0.75), (0.75, 0.80),
    (0.80, 0.85), (0.85, 0.88), (0.88, 0.90), (0.90, 0.92),
    (0.92, 0.94), (0.94, 0.96), (0.96, 0.98), (0.98, 1.01),
]


def jaccard_to_hicat_level(val):
    """Map Jaccard similarity (0–1) to HiCAT 12-level colour index."""
    for i, (lo, hi) in enumerate(HICAT_THRESHOLDS):
        if lo <= val < hi:
            return i + 1
    return 12


# =============================================================================
# Data parsers
# =============================================================================

def parse_mono_csv(path):
    records = []
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return records
        n_cols = len(header)
        for row in reader:
            if len(row) < n_cols:
                continue
            if n_cols == 6:
                raw_name = row[1]
                strand = "-" if raw_name.endswith("'") else "+"
                mono_name = raw_name.rstrip("'")
                records.append({
                    "seq_name": row[0], "mono_name": mono_name, "strand": strand,
                    "start": int(row[2]), "end": int(row[3]),
                    "identity": float(row[4]), "length": int(row[5]),
                })
            elif n_cols == 7:
                records.append({
                    "seq_name": row[0], "mono_name": row[1], "strand": row[2],
                    "start": int(row[3]), "end": int(row[4]),
                    "identity": float(row[5]), "length": int(row[6]),
                })
    return records


def parse_hors(path):
    hors = []
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return hors
        for row in reader:
            if len(row) < 12:
                continue
            hors.append({
                "seq_name": row[0], "region": row[1], "pattern": row[2],
                "compressed_pattern": row[3] if len(row) > 3 else row[2],
                "start_pos": int(row[8]), "end_pos": int(row[9]),
                "copies": int(row[10]),
                "mean_identity": float(row[11]) if len(row) > 11 else 0.0,
                "hor_len_monomers": int(row[5]),
            })
    return hors


def parse_centro_bed(path):
    regions = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            regions.append({
                "chrom": parts[0], "start": int(parts[1]), "end": int(parts[2]),
                "name": parts[3] if len(parts) > 3 else "",
                "score": parts[4] if len(parts) > 4 else "0",
            })
    return regions


# =============================================================================
# 1. Monomer array color-band map
# =============================================================================

def assign_mono_colors(records):
    mono_types = sorted(set(r["mono_name"] for r in records))
    n = len(mono_types)
    cmap = plt.get_cmap("tab20") if n <= 20 else plt.get_cmap("hsv")
    norm = n if n > 20 else 20
    colors = [cmap(i / norm) for i in range(n)]
    return {m: colors[i] for i, m in enumerate(mono_types)}


def darken(color, factor=0.55):
    r, g, b, a = color
    return (r * factor, g * factor, b * factor, a)


def draw_monomer_array(records, out_prefix, fmt="png",
                       hor_records=None, centro_regions=None,
                       fig_width_per_mb=12, dpi=300):
    if not records:
        print("[viz] No monomer records to visualize.", file=sys.stderr)
        return
    hor_records = hor_records or []
    centro_regions = centro_regions or []

    seq_groups = {}
    for r in records:
        seq_groups.setdefault(r["seq_name"], []).append(r)
    hor_by_seq = {}
    for h in hor_records:
        hor_by_seq.setdefault(h["seq_name"], []).append(h)
    centro_by_seq = {}
    for c in centro_regions:
        centro_by_seq.setdefault(c["chrom"], []).append(c)
    seq_names = sorted(seq_groups.keys())
    color_map = assign_mono_colors(records)

    n_seqs = len(seq_names)
    fig_height = max(4, n_seqs * 2.5)
    max_coord = max(r["end"] for r in records)
    fig_width = max(10, fig_width_per_mb * (max_coord / 1e6))
    fig, axes = plt.subplots(n_seqs, 1, figsize=(fig_width, fig_height),
                             squeeze=False, sharex=True)

    for idx, seq_name in enumerate(seq_names):
        ax = axes[idx, 0]
        recs = sorted(seq_groups[seq_name], key=lambda x: x["start"])
        hors = hor_by_seq.get(seq_name, [])
        centros = centro_by_seq.get(seq_name, [])

        for r in recs:
            base_color = color_map[r["mono_name"]]
            col = darken(base_color) if r["strand"] == "-" else base_color
            ax.barh(0.5, r["end"] - r["start"] + 1,
                    left=r["start"], height=0.8,
                    color=col, edgecolor="none", linewidth=0)

        if hors:
            for h in hors:
                ax.axvspan(h["start_pos"], h["end_pos"],
                           ymin=0.78, ymax=0.98,
                           color="crimson", alpha=0.18, zorder=2)

        if centros:
            for c in centros:
                ax.axvspan(c["start"], c["end"],
                           ymin=0.02, ymax=0.18,
                           color="limegreen", alpha=0.35, zorder=2)
                mid = (c["start"] + c["end"]) / 2
                ax.annotate("CEN", xy=(mid, 0.02), fontsize=8,
                            ha="center", va="bottom", color="darkgreen",
                            fontweight="bold", zorder=3)

        ax.set_ylim(0, 1)
        ax.set_yticks([0.5])
        ax.set_yticklabels([seq_name], fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Genomic position (bp)", fontsize=11)

        legend_patches = []
        if hors:
            legend_patches.append(mpatches.Patch(
                color="crimson", alpha=0.3,
                label=f"HOR region (n={len(hors)})"))
        if centros:
            legend_patches.append(mpatches.Patch(
                color="limegreen", alpha=0.4,
                label="Predicted CenSatArray"))
        if legend_patches:
            ax.legend(handles=legend_patches, loc="upper right",
                      fontsize=8, frameon=True, fancybox=True)

    all_starts = [r["start"] for r in records]
    all_ends = [r["end"] for r in records]
    pad = max(1, int((max(all_ends) - min(all_starts)) * 0.01))
    axes[0, 0].set_xlim(min(all_starts) - pad, max(all_ends) + pad)

    mono_legend = [mpatches.Patch(color=color_map[m], label=str(m))
                   for m in sorted(color_map.keys())]
    n_legend_cols = min(12, len(mono_legend))
    fig.legend(handles=mono_legend, loc="lower center",
               ncol=n_legend_cols, fontsize=8,
               title="Monomer types", title_fontsize=10,
               frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.suptitle("centroAnno Monomer Array Visualization", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    out_path = f"{out_prefix}_monomerArray.{fmt}"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[viz] Saved: {out_path}")


# =============================================================================
# 2. Tandem-repeat region map (cautils style)
# =============================================================================

def compute_repeat_regions(records, min_len=1000):
    """Merge consecutive monomer blocks into repeat regions per sequence,
    reporting the mode length for each region."""
    if not records:
        return []

    # Group by seq_name
    by_seq = {}
    for r in records:
        by_seq.setdefault(r["seq_name"], []).append(r)

    regions = []
    for seq_name, recs in by_seq.items():
        recs = sorted(recs, key=lambda r: r["start"])
        region_st = recs[0]["start"]
        region_ed = recs[0]["end"]
        len_list = [recs[0]["length"]]

        for r in recs[1:]:
            if abs(region_ed - r["start"]) < 100:
                region_ed = r["end"]
                len_list.append(r["length"])
            else:
                if region_ed - region_st >= min_len:
                    mode_len, _ = Counter(len_list).most_common(1)[0]
                    regions.append([seq_name, region_st, region_ed, mode_len])
                region_st = r["start"]
                region_ed = r["end"]
                len_list = [r["length"]]

        if region_ed - region_st >= min_len:
            mode_len, _ = Counter(len_list).most_common(1)[0]
            regions.append([seq_name, region_st, region_ed, mode_len])
    return regions


def draw_cautils_mono(regions, out_prefix, fmt="png"):
    """Draw cautils-style broken_barh figure colored by repeat-unit length."""
    if not regions:
        print("[viz] No repeat regions for mono map.", file=sys.stderr)
        return

    seq_names = sorted(set(r[0] for r in regions))
    seq_y = {s: i * 20 for i, s in enumerate(seq_names)}

    lengths = [r[3] for r in regions]
    min_len, max_len = min(lengths), max(lengths)
    norm = mcolors.Normalize(vmin=min_len, vmax=max(max_len, min_len + 1))
    cmap = plt.colormaps["viridis"]

    fig, ax = plt.subplots(figsize=(24, max(4, len(seq_names) * 1.5 + 1)))
    for r in regions:
        seq, st, ed, rep_len = r
        ypos = seq_y[seq]
        color = cmap(norm(rep_len))
        ax.broken_barh([(st, ed - st)], (ypos, 8), facecolors=color)

    for seq in seq_names:
        ypos = seq_y[seq]
        ax.text(0, ypos + 4, seq, va="center", ha="right",
                fontsize=11, fontweight="bold")

    ax.set_ylim(-10, max(seq_y.values()) + 20)
    ax.set_xlim(0, max(r[2] for r in regions) + 300)
    ax.set_xlabel("Genomic Coordinate (bp)", fontsize=13)
    ax.set_yticks([])
    ax.set_title("Tandem Repeat Unit Length Map", fontsize=15, fontweight="bold")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("Repeat Unit Length (bp)", fontsize=12)

    plt.tight_layout()
    out_path = f"{out_prefix}_repeatRegions.{fmt}"
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[viz] Saved: {out_path}")


# =============================================================================
# 3. HOR region map (cautils style)
# =============================================================================

def draw_cautils_hor(hors, out_prefix, fmt="png"):
    """Draw cautils-style broken_barh figure colored by HOR length (bp)."""
    if not hors:
        print("[viz] No HOR records for HOR map.", file=sys.stderr)
        return

    # Merge overlapping/adjacent HORs of the same pattern per sequence
    seq_hors = {}
    for h in hors:
        seq_hors.setdefault(h["seq_name"], []).append(h)

    merged = []
    for seq, hlist in seq_hors.items():
        hlist.sort(key=lambda x: x["start_pos"])
        for h in hlist:
            hor_len_bp = h["end_pos"] - h["start_pos"]
            if merged and merged[-1][0] == seq and merged[-1][4] == h["pattern"]:
                if abs(merged[-1][2] - h["start_pos"]) < 100:
                    merged[-1][2] = h["end_pos"]
                    merged[-1][3] = max(merged[-1][3], hor_len_bp)
                    continue
            merged.append([seq, h["start_pos"], h["end_pos"], hor_len_bp, h["pattern"]])

    seq_names = sorted(set(r[0] for r in merged))
    seq_y = {s: i * 20 for i, s in enumerate(seq_names)}

    lengths = [r[3] for r in merged]
    min_len, max_len = min(lengths), max(lengths)
    norm = mcolors.Normalize(vmin=min_len, vmax=max(max_len, min_len + 1))
    cmap = plt.colormaps["plasma"]

    fig, ax = plt.subplots(figsize=(24, max(4, len(seq_names) * 1.5 + 1)))
    for r in merged:
        seq, st, ed, hor_len, _ = r
        ypos = seq_y[seq]
        color = cmap(norm(hor_len))
        ax.broken_barh([(st, ed - st)], (ypos, 8), facecolors=color)

    for seq in seq_names:
        ypos = seq_y[seq]
        ax.text(0, ypos + 4, seq, va="center", ha="right",
                fontsize=11, fontweight="bold")

    ax.set_ylim(-10, max(seq_y.values()) + 20)
    ax.set_xlim(0, max(r[2] for r in merged) + 300)
    ax.set_xlabel("Genomic Coordinate (bp)", fontsize=13)
    ax.set_yticks([])
    ax.set_title("HOR Length Map", fontsize=15, fontweight="bold")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation="vertical")
    cb.set_label("HOR Span (bp)", fontsize=12)

    plt.tight_layout()
    out_path = f"{out_prefix}_horRegions.{fmt}"
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[viz] Saved: {out_path}")


# =============================================================================
# Shared helpers for sequence similarity
# =============================================================================

def kmer_set(seq, k=13):
    """Return set of k-mers for a sequence."""
    s = seq.upper()
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


# =============================================================================
# 4. HiCAT-style HOR tracks + diamond self-similarity heatmap
# =============================================================================

def draw_hicat_style(fasta_path, hors, out_prefix,
                     window_size=5000, max_windows=200,
                     show_n=5, diamond_scale=1.15, fmt="png", dpi=300,
                     centro_regions=None):
    """Draw HiCAT-style composite figure for predicted CenSatArray regions:
      • Upper-triangular diamond self-similarity heatmap (moddotplot style)
        forming an isosceles triangle with horizontal base at the bottom
      • Top-N HOR pattern tracks directly below the heatmap base
      • 12-level discrete colour scale (purple → dark red)

    If centro_regions is provided, only the best-scoring region per sequence
    is visualised; otherwise falls back to the full sequence.
    """
    if not hors:
        print("[viz] No HOR records for composite plot.", file=sys.stderr)
        return
    if not os.path.isfile(fasta_path):
        print(f"[viz] FASTA not found: {fasta_path}", file=sys.stderr)
        return

    seq_dict = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[rec.id] = str(rec.seq)

    seq_hors = {}
    for h in hors:
        seq_hors.setdefault(h["seq_name"], []).append(h)

    # Map: seq_name -> best CenSatArray region (by score descending)
    seq_best_region = {}
    if centro_regions:
        for reg in centro_regions:
            chrom = reg["chrom"]
            try:
                score = float(reg.get("score", 0))
            except ValueError:
                score = 0.0
            if chrom not in seq_best_region or score > seq_best_region[chrom]["score"]:
                seq_best_region[chrom] = {
                    "start": reg["start"], "end": reg["end"],
                    "name": reg.get("name", ""), "score": score
                }

    for seq_name, hlist in seq_hors.items():
        if seq_name not in seq_dict:
            print(f"[viz] Sequence '{seq_name}' not in FASTA, skipping composite plot.",
                  file=sys.stderr)
            continue

        full_seq = seq_dict[seq_name].upper()
        full_len = len(full_seq)

        # Determine region to plot
        if seq_name in seq_best_region:
            reg = seq_best_region[seq_name]
            reg_start = max(0, reg["start"])
            reg_end = min(full_len, reg["end"])
            if reg_end <= reg_start:
                print(f"[viz] Invalid CenSatArray region for {seq_name}, skipping.",
                      file=sys.stderr)
                continue
            seq = full_seq[reg_start:reg_end]
            seq_len = len(seq)
            region_name = reg["name"] or f"{reg_start}-{reg_end}"
            plot_label = f"{seq_name}_region_{region_name}"
        else:
            # Fallback: plot full sequence when no CenSatArray prediction
            seq = full_seq
            seq_len = full_len
            reg_start = 0
            plot_label = seq_name

        from collections import defaultdict
        pat_data = defaultdict(lambda: {"total_copies": 0, "spans": [], "compressed": ""})
        for h in hlist:
            # Only keep HORs that overlap the selected region
            h_st = h["start_pos"]
            h_ed = h["end_pos"]
            if h_ed <= reg_start or h_st >= reg_start + seq_len:
                continue
            p = h["pattern"]
            # Adjust coordinates relative to region start
            adj_st = max(0, h_st - reg_start)
            adj_ed = min(seq_len, h_ed - reg_start)
            pat_data[p]["total_copies"] += h["copies"]
            pat_data[p]["spans"].append((adj_st, adj_ed))
            pat_data[p]["compressed"] = h.get("compressed_pattern", p) or p
            pat_data[p]["total_span"] = sum(ed - st for st, ed in pat_data[p]["spans"])

        sorted_pats = sorted(
            pat_data.keys(),
            key=lambda p: pat_data[p]["total_span"],
            reverse=True
        )
        filter_pats = sorted_pats[:show_n]
        n_tracks = len(filter_pats)
        if not filter_pats:
            continue

        n_win = seq_len // window_size
        if n_win > max_windows:
            window_size = max(window_size, int(np.ceil(seq_len / max_windows)))
            n_win = seq_len // window_size
        if n_win < 10:
            print(f"[viz] {seq_name} too short ({n_win} windows) for HiCAT plot, skipping.",
                  file=sys.stderr)
            continue

        k = 15
        windows = [seq[i * window_size:(i + 1) * window_size] for i in range(n_win)]
        kmer_sets = [kmer_set(w, k) for w in windows]
        mat = np.zeros((n_win, n_win), dtype=np.float32)
        for i in range(n_win):
            mat[i, i] = 1.0
            for j in range(i + 1, n_win):
                jacc = jaccard(kmer_sets[i], kmer_sets[j])
                # Convert Jaccard to estimated sequence identity (ANI).
                # Mash formula: ANI ≈ (2*J/(1+J))^(1/k)
                if jacc > 0:
                    sim = min(1.0, (2.0 * jacc / (1.0 + jacc)) ** (1.0 / k))
                else:
                    sim = 0.0
                mat[i, j] = sim
                mat[j, i] = sim

        # ---- Plot ----
        fig, ax = plt.subplots(figsize=(12, 12))
        base_sequence_len = seq_len
        track_color = "#D14524"
        ws = window_size

        # -- Diamond heatmap: moddotplot coordinates --
        # cx = (col + row) * ws / 2,  cy = (row - col) * ws / 2
        # Upper-triangle (row >= col) becomes an isosceles triangle
        # with horizontal base at y = 0  (the diagonal)
        # Use PolyCollection for smooth, gap-free rendering with continuous colours.
        from matplotlib.collections import PolyCollection
        cmap = plt.cm.turbo
        # Use PowerNorm to stretch the high-identity band (0.90–1.0).
        # gamma=2 pulls apart the 95–100% region so subtle differences are visible.
        norm = plt.matplotlib.colors.PowerNorm(gamma=2, vmin=0.90, vmax=1.0)
        verts = []
        facecolors = []
        for col_idx in range(n_win):
            for row_idx in range(col_idx, n_win):
                sim = mat[col_idx, row_idx]
                cx = (col_idx + row_idx) * ws / 2
                cy = (row_idx - col_idx) * ws / 2
                ds = ws * diamond_scale / 2
                verts.append([
                    (cx, cy - ds),
                    (cx + ds, cy),
                    (cx, cy + ds),
                    (cx - ds, cy),
                ])
                facecolors.append(cmap(norm(sim)))
        poly = PolyCollection(
            verts, closed=True, facecolors=facecolors,
            edgecolors='none', antialiaseds=False, linewidths=0
        )
        ax.add_collection(poly)

        # -- HOR tracks (below heatmap, with extra gap to avoid overlap) --
        track_h = base_sequence_len / 80
        track_gap = base_sequence_len / 35   # wider gap so labels don't overlap
        heatmap_gap = track_gap * 1.5   # extra space between diamond base and first track
        for idx, pat in enumerate(filter_pats):
            track_y = -heatmap_gap - (idx + 1) * track_gap
            rect = mpatches.Rectangle(
                (0, track_y), base_sequence_len, track_h,
                color="#D0CECE"
            )
            ax.add_patch(rect)
            for st, ed in pat_data[pat]["spans"]:
                rect = mpatches.Rectangle(
                    (st, track_y), ed - st, track_h,
                    color=track_color, lw=0
                )
                ax.add_patch(rect)
            label = pat_data[pat]["compressed"] or pat
            plt.text(
                base_sequence_len + base_sequence_len / 100,
                track_y + track_h / 2,
                label, fontsize=9, va="center"
            )

        # -- Scale bar (below HOR tracks, MB units) --
        scale_y = -heatmap_gap - (n_tracks + 1) * track_gap - base_sequence_len / 200
        ax.add_patch(mpatches.Rectangle(
            (0, scale_y), base_sequence_len, base_sequence_len / 1000,
            color="black"
        ))
        # Dynamic tick count based on sequence length (avoid duplicate MB labels)
        if base_sequence_len < 1e6:
            n_ticks = 4
        elif base_sequence_len < 5e6:
            n_ticks = 6
        else:
            n_ticks = 11
        point_bar = base_sequence_len / (n_ticks - 1)
        for i in range(n_ticks):
            x = i * point_bar
            ax.add_patch(mpatches.Rectangle(
                (x, scale_y), base_sequence_len / 1000,
                -base_sequence_len / 100, color="black"
            ))
            mb_label = f"{x / 1e6:.1f}M"
            ax.text(
                x, scale_y - base_sequence_len / 35,
                mb_label, fontsize=7, ha="center"
            )

        # -- Colour legend (below scale bar): continuous turbo gradient 90-100% --
        legend_top = scale_y - base_sequence_len / 15
        legend_h = base_sequence_len / 80
        legend_w_total = base_sequence_len / 2.5
        n_grad = 50
        grad_w = legend_w_total / n_grad
        norm = plt.matplotlib.colors.PowerNorm(gamma=2, vmin=0.90, vmax=1.0)
        cmap = plt.cm.turbo
        for i in range(n_grad):
            lx = i * grad_w
            frac = i / (n_grad - 1)
            val = 0.90 + frac * 0.10
            c = cmap(norm(val))
            ax.add_patch(mpatches.Rectangle(
                (lx, legend_top - legend_h), grad_w, legend_h,
                color=c, ec="none", lw=0
            ))
        ax.add_patch(mpatches.Rectangle(
            (0, legend_top - legend_h), legend_w_total, legend_h,
            fill=False, ec="black", lw=0.5
        ))
        # Ticks at 90, 93, 96, 100 (spaced by PowerNorm gamma=2)
        for val, label in [(0.90, "90"), (0.93, "93"), (0.96, "96"), (1.0, "100")]:
            frac = ((val - 0.90) / 0.10) ** 2
            x = frac * legend_w_total
            ax.add_patch(mpatches.Rectangle(
                (x, legend_top - legend_h), base_sequence_len / 1000,
                -base_sequence_len / 200, color="black"
            ))
            ax.text(x, legend_top - legend_h * 2.5,
                    label, fontsize=6, ha="center")
        ax.text(legend_w_total / 2, legend_top - legend_h * 4.5,
                "Identity (%)", fontsize=8, ha="center", fontweight="bold")

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis("equal")

        safe_name = plot_label.replace(":", "_").replace("/", "_")
        out_path = f"{out_prefix}_composite_{safe_name}.{fmt}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight",
                    pad_inches=0.1, facecolor="white")
        plt.close(fig)
        print(f"[viz] Saved: {out_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive visualization for centroAnno outputs."
    )
    parser.add_argument("--mono-csv", required=True,
                        help="centroAnno *_decomposedResult.csv")
    parser.add_argument("--out-prefix", required=True,
                        help="Output file prefix (path + basename)")
    parser.add_argument("--mode", default="png",
                        choices=["png", "svg"],
                        help="Output format for all figures [default: png]")
    parser.add_argument("--hor-csv", default="",
                        help="detectHORs *_allHORs.csv")
    parser.add_argument("--centro-bed", default="",
                        help="centroFinder *_censatarray_candidates.bed")
    parser.add_argument("--fasta", default="",
                        help="Reference FASTA (required for composite plot)")
    parser.add_argument("--width-per-mb", type=float, default=12,
                        help="Figure width inches per Mb [default: 12]")
    parser.add_argument("--dpi", type=int, default=300,
                        help="DPI for raster formats [default: 300]")
    parser.add_argument("--composite", action="store_true",
                        help="Generate composite HOR + diamond heatmap")
    parser.add_argument("--composite-window", type=int, default=5000,
                        help="Composite diamond heatmap window size (bp) [default: 5000]")
    parser.add_argument("--composite-max-win", type=int, default=100,
                        help="Composite max windows [default: 100]")
    parser.add_argument("--composite-top-n", type=int, default=5,
                        help="Composite top-N HOR patterns to show [default: 5]")
    parser.add_argument("--composite-diamond-scale", type=float, default=1.0,
                        help="Diamond size scale factor (1.0=touching, <1=smaller, >1=overlap) [default: 1.0]")
    args = parser.parse_args()

    records = parse_mono_csv(args.mono_csv)
    if not records:
        print("[viz] ERROR: No records parsed from monomer CSV.", file=sys.stderr)
        sys.exit(1)

    hors = parse_hors(args.hor_csv) if args.hor_csv else []
    centros = parse_centro_bed(args.centro_bed) if args.centro_bed else []

    # 1. Repeat-region map (cautils style)
    repeat_regions = compute_repeat_regions(records)
    draw_cautils_mono(repeat_regions, args.out_prefix, fmt=args.mode)

    # 2. HOR map (cautils style)
    draw_cautils_hor(hors, args.out_prefix, fmt=args.mode)

    # 3. Composite plot (CenSatArray region only when --centro-bed given)
    if args.composite and args.fasta and hors:
        draw_hicat_style(args.fasta, hors, args.out_prefix,
                         window_size=args.composite_window,
                         max_windows=args.composite_max_win,
                         show_n=args.composite_top_n,
                         diamond_scale=args.composite_diamond_scale,
                         fmt=args.mode, dpi=args.dpi,
                         centro_regions=centros)
    elif args.composite and not args.fasta:
        print("[viz] Skipping composite plot: --fasta not provided.", file=sys.stderr)
    elif args.composite and not hors:
        print("[viz] Skipping composite plot: --hor-csv not provided.", file=sys.stderr)


if __name__ == "__main__":
    main()
