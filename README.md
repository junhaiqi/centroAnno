# centroAnno1.5

<p align="center">
  <img src="assets/centroAnno_logo.png" width="180" alt="centroAnno logo">
</p>

<p align="center">
  <b>Reference-free, end-to-end centromeric satellite annotation</b>
</p>

<p align="center">
  <a href="https://github.com/junhaiqi/centroAnno/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <img src="https://img.shields.io/badge/language-C%2B%2B11%20%7C%20Python3-blue" alt="C++11 / Python3">
  <img src="https://img.shields.io/badge/platform-Linux-green" alt="Linux">
</p>

---

**centroAnno** is a reference-free toolkit for the *de novo* identification, monomer-level decomposition, higher-order repeat (HOR) inference, and centromeric satellite array (CenSatArray) detection of centromere-associated satellite DNA directly from whole-genome assemblies, isolated centromeric subregions, or raw long-read sequencing data.

Unlike existing tools that require predefined monomer libraries, manually extracted centromeric sequences, or external post-processing, centroAnno performs **end-to-end hierarchical reconstruction** in a single command.

## ✨ Why centroAnno?

Centromeres are dominated by megabase-scale tandem repeats whose sequence organization is notoriously difficult to resolve. Existing methods typically address only *part* of the problem:

| Task | Typical limitation of existing tools | centroAnno |
|---|---|---|
| **Array localization** | Requires prior centromere coordinates or epigenetic maps | ✅ *De novo* from whole chromosomes |
| **Monomer inference** | Needs curated monomer templates | ✅ Reference-free inference |
| **HOR reconstruction** | Works only on pre-extracted regions | ✅ Runs on whole chromosomes |
| **Local-Nested HOR detection** | Not supported or requires manual analysis | ✅ local-nested HOR parsing |
| **Visualization** | Multiple disconnected tools | ✅ Publication-ready composite figures in one run |

In short, centroAnno unifies what previously required a patchwork of specialized tools and manual steps into one automated pipeline.

### Key Highlights

- 🔬 **Reference-free monomer inference** — no monomer library, repeat family, or HOR template required.
- 🧬 **Three analysis modes**:
  - `anno-asm`: whole chromosomes / assemblies (default)
  - `anno-sat-asm`: pre-extracted centromeric sequences (HiCAT/HORmon-style input)
  - `anno-read`: noisy long reads
- 🏗️ **Hierarchical HOR reconstruction** — canonical HORs, local-nested HORs (`A H^k B`) with compressed pattern output.
- 🎯 **CenSatArray detection** — literature-informed, structure-based scoring predicts the primary centromeric satellite array from monomer decomposition alone.
- 📊 **Publication-ready visualizations** — HiCAT-style composite plots combining diamond self-similarity heatmaps, HOR tracks, and coordinate scales.
- ⚡ **Fast** — *A. thaliana* in 18 minutes.
- 🌿 **Cross-species** — validated on human (T2T-CHM13) and plant (*Arabidopsis thaliana*) genomes.
- 🔧 **Modular & scriptable** — run the full pipeline with one command, or call each stage independently.

---

## 📦 Installation

```bash
git clone --recursive https://github.com/junhaiqi/centroAnno.git
cd centroAnno
./install.sh
```

> **One-command install**: `./install.sh` detects your OS, installs missing system/Python dependencies, initializes the `spoa` submodule, and compiles the binary. Use `--conda` to force conda for Python packages, or `--no-sudo` if you prefer manual system setup.

> If you cloned without `--recursive`, pull the submodule manually:
> ```bash
> git submodule update --init --recursive
> ```

### Requirements

| Component | Requirement |
|-----------|-------------|
| **centroAnno** (C++ core) | Linux, g++ ≥ 9.3.0, C++11, OpenMP, zlib |
| **detectHORs.py** | Python 3, numpy, pandas |
| **centroFinder.py** | Python 3, biopython, numpy, scipy, matplotlib, pandas |
| **visualizeMonomerArray.py** | Python 3, matplotlib, numpy |

```bash
pip install biopython numpy scipy matplotlib pandas
```

After compilation you will have:
- `./centroAnno` — the core monomer-decomposition binary
- `./run_centroAnno.sh` — the unified shell pipeline (recommended entry point)
- `./continue_pipeline.sh` — resume Stages 2–5 from existing CSVs
- `script/detectHORs.py` — HOR detection
- `script/centroFinder.py` — CenSatArray detection
- `script/visualizeMonomerArray.py` — monomer-array visualization
- `script/extractBedRegions.py` — BED region extraction

---

## 🚀 Quick Start

### 1. Whole chromosome or assembly (recommended)

```bash
./run_centroAnno.sh -i chr1.fa -o results/chr1 -t 16
```

### 2. Pre-extracted centromeric sequence (HiCAT/HORmon-like input)

```bash
./run_centroAnno.sh -i example/cen21.fa -o results/cen21 -x anno-sat-asm
```

### 3. Raw long reads

```bash
./run_centroAnno.sh -i reads.fa -o results/reads -x anno-read -t 8
```

When the pipeline finishes, check:
- `results/chr1/03_cenSatArray/*_censatarray_signals.png` — CenSatArray detection signals
- `results/chr1/04_visualization/*_composite_*.png` — HiCAT-style HOR tracks + diamond heatmap
- `results/chr1/05_cenSatArray_refined/` — high-resolution re-annotation of the top CenSatArray (optional, enabled by default in `anno-asm`)

---

## 📖 Usage

### `run_centroAnno.sh` — full end-to-end pipeline

```bash
./run_centroAnno.sh -i <in.fa> -o <out_dir> [options]
```

#### Main options

| Option | Description |
|--------|-------------|
| `-i FILE` | Input FASTA/FASTQ file **(required)** |
| `-o DIR` | Output directory **(required)** |
| `-x MODE` | Analysis mode: `anno-asm` (default), `anno-sat-asm`, or `anno-read` |
| `-m FILE` | Monomer template FASTA (**anno-asm only**; optional; skips *de novo* inference) |
| `-t INT` | Threads [default: 8] |
| `-k INT` | k-mer size [default: 10] |
| `-f FLOAT` | FPS cutoff [default: 0.6] |
| `-r FLOAT` | Repeat ratio cutoff [default: 0.3] |
| `-w INT` | Window size for template inference [default: 500000] |
| `-c INT` | Homopolymer compression (1=yes, 0=no) [default: 1] |
| `-e FLOAT` | DBSCAN identity cutoff [default: 0.95] |
| `-M INT` | Max monomers in a HOR [default: 50] |
| `-L INT` | Minimum sequence length to annotate [default: 5000] |
| `-A INT` | Max region length (**anno-asm only**) [default: 1000000] |
| `-N INT` | Min region length (genome mode) [default: 100] |
| `-F FLOAT` | Genome annotation identity cutoff [default: 0.8] |

#### Stage options

| Option | Default | Description |
|--------|---------|-------------|
| `--skip-hors` | — | Skip HOR detection (Stage 2) |
| `--skip-centrofinder` | — | Skip CenSatArray detection (Stage 3) |
| `--centrofinder` | — | Force CenSatArray detection in `anno-sat-asm` / `anno-read` |
| `--skip-visualize` | — | Skip visualization (Stage 4) |
| `--no-refine-centromere` | — | Skip Stage 5 high-resolution refinement (refinement is **enabled** by default in `anno-asm`) |
| `--refine-centromere` | — | Explicitly enable refinement |

#### Examples

```bash
# Skip HOR detection, only monomers + CenSatArray
./run_centroAnno.sh -i chr1.fa -o out/ --skip-hors

# Skip CenSatArray detection, only monomers + HORs
./run_centroAnno.sh -i chr1.fa -o out/ --skip-centrofinder

# Force CenSatArray detection even on a subregion
./run_centroAnno.sh -i cen.fa -o out/ -x anno-sat-asm --centrofinder

# Whole-genome run with high thread count and larger region limit
./run_centroAnno.sh -i genome.fa -o out/ -t 32 -A 2000000

# Subregion FASTA — coordinate offset is handled automatically
./run_centroAnno.sh -i cen21.fa -o out/ -x anno-sat-asm
```

> **Note on `--no-offset`**: For `anno-sat-asm` and `anno-read` modes the pipeline **automatically** passes `--no-offset` to `centroFinder.py` because the input FASTA is already the target subregion/reads. Additionally, if the FASTA sequence name matches the pattern `chr:start-end` (e.g. `CP068257.1:11699867-12031015`), `--no-offset` is also injected automatically in `anno-asm` mode, preventing coordinate overflow on subregion inputs. For `anno-asm` with whole-chromosome FASTA the offset is kept to allow coordinate conversion.

---

### `continue_pipeline.sh` — resume from existing Stage-1 CSVs

If you already ran Stage 1 (`centroAnno`) and have `*_decomposedResult.csv` files, use this script to run Stages 2–5 without re-running the C++ core.

```bash
./continue_pipeline.sh -d <bench_dir> [ -F <fasta_dir> | -f <fasta> ] [options]
```

#### FASTA requirement is stage-dependent

| Stage | Needs FASTA? | Notes |
|-------|--------------|-------|
| Stage 2 — HOR detection | ❌ No | Only needs the CSV |
| Stage 3 — CenSatArray detection | ✅ Yes | Uses the original input FASTA |
| Stage 4 — Visualization | ❌ No (basic) / ✅ Yes (composite) | Basic visualization needs only CSV + Stage 2 output |
| Stage 5 — Refinement | ✅ Yes | Extracts top CenSatArray from the original FASTA |

So **`-F`/`-f` is optional if you skip Stage 3 and Stage 5**.

#### Options

| Option | Description |
|--------|-------------|
| `-d DIR` | Directory containing `*_decomposedResult.csv` files **(required)** |
| `-F DIR` | Directory with source FASTA files |
| `-f FILE` | Single FASTA file to use for **all** CSVs (e.g. a multi-sequence FASTA) |
| `-j INT` | Parallel samples [default: 1] |
| `-t INT` | Threads for anno-sat-asm refinement [default: 16] |
| `--mode` | Visualization format: `png` or `svg` [default: `png`] |
| `--max-depth INT` | Max depth to search for CSVs [default: 2] |
| `--target-regex REGEX` | Only process CSVs whose basename matches the regex [default: `.*`] |
| `--skip-hors` | Skip Stage 2 |
| `--skip-centrofinder` | Skip Stage 3 |
| `--skip-visualize` | Skip Stage 4 |
| `--skip-refine` | Skip Stage 5 |
| `--strict` | Exit immediately if any stage command fails |
| `--dry-run` | Preview what would be run without executing |

#### Examples

```bash
# Run HOR detection + visualization only (no FASTA needed)
./continue_pipeline.sh -d results/ --skip-centrofinder --skip-refine

# Full Stages 2–5 with FASTA directory and 8 parallel jobs
./continue_pipeline.sh -d results/ -F Data/ -j 8 -t 16

# Single multi-sequence FASTA for all CSVs
./continue_pipeline.sh -d results/ -f genome.fa -j 4

# Only canonical human chromosomes
./continue_pipeline.sh -d results/ -F Data/ --target-regex '^chr([0-9]+|X|Y)$'

# Preview what would be processed
./continue_pipeline.sh -d results/ -F Data/ --dry-run
```

`continue_pipeline.sh` recursively searches for FASTA files under `-F`, matches them to CSVs by basename or parent-directory name, and reports how many samples were processed or skipped.

---

## 🧩 Standalone Usage

Each component can also be run independently for fine-grained control.

### 1. `centroAnno` — monomer inference / decomposition

```bash
# Whole chromosome / genome assembly (default mode)
./centroAnno -o test_output -x anno-asm chr1.fa

# Centromeric alpha-satellite sequence / assembly
./centroAnno -o test_output -x anno-sat-asm example/cen21.fa

# With monomer templates (Currently, it can only be used under anno-sat-asm!)
./centroAnno -o test_output -m templates.fa -x anno-sat-asm example/cen21.fa

# Sequencing reads
./centroAnno -o test_output -x anno-read reads.fa

# Stream input from another tool (e.g. AGC decompression)
./centroAnno -o test_output -x anno-asm -

# Pipe from AGC
agc getctg in.agc ctg1 ctg2 | ./centroAnno -o test_output -x anno-asm -
```

Full option list:
```
  -o STR     Output folder [required]
  -m STR     Monomer template FASTA [default: None]
  -k INT     k-mer size [default: 10]
  -f FLOAT   FPS cutoff [default: 0.6]
  -r FLOAT   Repeat ratio cutoff [default: 0.3]
  -w INT     Window size for template inference [default: 500000]
  -c BOOL    Homopolymer compression (1=yes, 0=no) [default: 1]
  -e FLOAT   DBSCAN identity cutoff [default: 0.95]
  -x STR     Mode: anno-asm | anno-sat-asm | anno-read [default: anno-asm]
  -t INT     Threads [default: 8]
  -M INT     Max monomers in HOR [default: 50]
  -L INT     Minimum sequence length [default: 5000]
  -A INT     Max region length (genome) [default: 1000000]
  -N INT     Min region length (genome) [default: 100]
  -F FLOAT   Genome annotation identity cutoff [default: 0.8]
  -S BOOL    Only scan repeats, no annotation (1=yes, 0=no) [default: 0]
```

### 2. `detectHORs.py` — HOR detection

```bash
python script/detectHORs.py \
    -i test_output/*_decomposedResult.csv \
    -o hor_results/sample \
    -m 50 -c 2 --min-identity 0.90 --min-max-copies 2
```

Outputs:
- `*_allHORs.csv` — every HOR occurrence
- `*_allHORPatterns.txt` — distinct canonical patterns per region
- `*_allHORStats.txt` — statistics per (region, pattern), sorted by total span
- `*_allHORPatterns_compressed.txt` — compressed nested patterns

### 3. `centroFinder.py` — CenSatArray detection

```bash
python script/centroFinder.py \
    --mono-csv test_output/*_decomposedResult.csv \
    --fasta chr1.fa \
    --out-prefix cf_results/sample \
    --window 50000 --step 10000 --density-thr 0.8 \
    --min-mono-len 80 --mono-gap-bp 100000
```

Outputs:
- `*_censatarray_candidates.bed` — predicted CenSatArray regions
- `*_censatarray_report.txt` — ranked candidate report
- `*_censatarray_signals.png` — 4-panel signal plot (TR density, entropy, conservation, monomer length)

### 4. `visualizeMonomerArray.py` — visualization

```bash
python script/visualizeMonomerArray.py \
    --mono-csv test_output/*_decomposedResult.csv \
    --out-prefix viz/mono \
    --hor-csv hor_results/*_allHORs.csv \
    --centro-bed cf_results/*_censatarray_candidates.bed \
    --mode png
```

---

## 📁 Output Structure

When using the pipeline, outputs are organized as follows:

```
out_dir/
├── 01_centroAnno/                    # Stage 1: Monomer decomposition
│   ├── *_decomposedResult.csv         # Monomer block list
│   └── *_monomerTemplates.fa         # Inferred monomer consensus sequences
│
├── 02_HORs/                          # Stage 2: HOR detection
│   ├── *_allHORs.csv                 # All HOR occurrences
│   ├── *_allHORPatterns.txt          # Distinct canonical patterns
│   ├── *_allHORStats.txt             # Statistics per (region, pattern)
│   └── *_allHORPatterns_compressed.txt # Nested compressed patterns
│
├── 03_cenSatArray/                    # Stage 3: CenSatArray detection
│   ├── *_censatarray_candidates.bed   # BED-format predictions
│   ├── *_censatarray_report.txt       # Ranked candidate report
│   └── *_censatarray_signals.png      # Signal plots
│
├── 04_visualization/                 # Stage 4: Visualization
│   ├── *_repeatRegions.{mode}        # Tandem repeat unit length map
│   ├── *_horRegions.{mode}           # HOR span map
│   └── *_composite_*.png              # HiCAT-style HOR tracks + diamond heatmap
│
└── 05_cenSatArray_refined/            # Stage 5: Optional high-resolution refinement
    └── */                             # Per-top-candidate refined annotation
```

### File format details

**`*_decomposedResult.csv`** (anno-sat format, header present)
```
sequence name,monomer name,start position,end position,estimated identity,length
```
- `monomer name`: monomer index; `12'` means reverse complement of monomer 12

**`*_decomposedResult.csv`** (anno-asm format, no header, 7 columns)
```
seq,mono,strand,start,end,identity,length
```
- Absolute coordinates; monomer names prefixed with `Pos:start:end_` for subregion tracking

**`*_allHORStats.txt`**
```
Region	HOR pattern	Occurrences	Total copies	Total span (bp)	Max copies	Mean identity	Median identity
```
Sorted by: region position ASC → total_span DESC → occ DESC → max_copies DESC

**`*_censatarray_candidates.bed`**
```
chrom	start	end	name	score	strand
```

---

## 📊 Visualization

### HiCAT-Style Composite Plot (`*_composite_*.png`)

A unified figure inspired by [HiCAT](https://github.com/zhangjizhong-chn/HICAT) that combines:

1. **Diamond self-similarity heatmap** — upper-triangular arrangement where each diamond represents pairwise window identity. The 12-level discrete colour scale ranges from purple (low identity) through green/yellow to dark red (high identity).
2. **Top-N HOR pattern tracks** — horizontal bar chart below the heatmap showing where each top HOR pattern occurs.
3. **Genomic scale bar** — coordinate ticks at the bottom.
4. **Colour legend** — 12-tier identity key.

**Standalone usage**:
```bash
python script/visualizeMonomerArray.py \
    --mono-csv test_output/*_decomposedResult.csv \
    --hor-csv hor_results/*_allHORs.csv \
    --fasta chr1.fa \
    --out-prefix viz/composite \
    --composite \
    --composite-window 3000 \
    --composite-top-n 5
```

| Option | Default | Description |
|--------|---------|-------------|
| `--composite` | — | Enable composite HOR + diamond heatmap |
| `--composite-window` | 5000 | Window size (bp) for pairwise comparison |
| `--composite-max-win` | 100 | Max number of windows |
| `--composite-top-n` | 5 | Number of top HOR patterns to display |
| `--composite-diamond-scale` | 1.0 | Diamond size scale factor |

### CenSatArray Signal Plot (`*_censatarray_signals.png`)

`centroFinder.py` produces a 4-panel plot showing:
1. **TR density** — tandem repeat density per window
2. **k-mer entropy** — sequence complexity
3. **Conservation score** — RepBlock-MonoTemp identity
4. **Mean monomer length** — biological length filter (α-satellite ~170 bp)

---

## ⚡ Performance Notes

- **Multi-sequence FASTA**: centroAnno natively supports multi-sequence FASTA input. Each sequence is processed independently, producing per-sequence `*_decomposedResult.csv` files.
- **Memory**: Sequence decomposition is the most memory-intensive step. Memory scales roughly linearly with thread count. For human chromosome 1 (~250 Mb), peak memory is ~3.9 GB with 8 threads.
- **Speed**: Human Chr1 annotation completes in ~2 hour on a modern server; 

- **Scaling**: For very large genomes, use `-A` to limit the maximum region length processed at once, it can reduce memory usage!.

---

## 🧪 Example Datasets

The `example/` directory contains small test files:
- `example/cen21.fa` — a human centromeric alpha-satellite sequence
- `example/AlphaSat.fa` — example monomer templates
- `example/simulated_genome.fasta` — a small simulated tandem-repeat genome

Run a quick test:
```bash
./run_centroAnno.sh -i example/cen21.fa -o results/cen21 -x anno-sat-asm
```

---

## 📚 Citation

If you use centroAnno in your research, please cite:

> Junhai Qi, Junchi Ma, Zheng Han, Renmin Han, Ting Yu, Guojun Li. De novo annotation of centromere with centroAnno. *bioRxiv* 2025.02.19.639205; doi: https://doi.org/10.1101/2025.02.19.639205

---

## 🤝 Acknowledgments

We thank [Meng Zhou](https://github.com/zhoudreames) for building the centroAnno v1.02 Singularity image.

## 📄 License

MIT License.
