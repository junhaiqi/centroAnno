## Getting Started

```bash
# Install centroAnno by conda
conda install qjh96::centroanno

# Install centroAnno by source codes
git clone https://github.com/junhaiqi/centroAnno.git
cd centroAnno && make -j8  # C++11 required to compile

# Run on test data (Default: Analyze centromeric alpha-satellite sequences/centromeric alpha-satellite assemblies (like HiCAT/HORmon/GRMhor) directly)
./centroAnno example/cen21.fa -o test

# Analyze the tandem repeats/HORs from given **chromosomes/assemblies/general sequences** (No prior information required):
./centroAnno $your_genome.fa -o $your_output -x anno-asm

# An example command for analyzing the human chromosome 1:
./centroAnno $your_chr1.fasta -o $your_chr1_out -x anno-asm

# Analyze the structure of centromere without template information from given **centromeric alpha-satellite sequences/centromeric alpha-satellite assemblies**:
./centroAnno $your_centromere.fa -o $your_output -x anno-sat-asm

# Analyze the structure of centromere with template information from given **centromeric alpha-satellite sequences/centromeric alpha-satellite assemblies**:
./centroAnno $your_centromere.fa -m $your_templates.fa -o $your_output -x anno-sat-asm

# Analyze the structure of centromere from given **sequencing reads**:
./centroAnno $your_sequencing_reads.fa -o $your_output -x anno-read

# Analyze and summarize the output of centroAnno, details in misc/misc/README.md:
python misc/misc/cautils.py $centroAnno_output_dir $your_analysis_dir
```


## Overview of centroAnno
centroAnno is a prior-independent tool for automatic and efficient centromere/tendem repeat structural analysis across multiple species. centroAnno supports the analysis of repeat units and higher-order tandem repeat units (HORs) in genome/assembly, centromere sequence, and single sequencing long read. In addition, we built a pipeline based on centroAnno to analyze repeats and HORs from noisy sequencing data, named [CentroRepeatAnalyzer](https://github.com/junhaiqi/CentroRepeatAnalyzer.git).
## Table of contents

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [Output](#output)
  * [Acknowledgments](#acknowledgments)
  * [License](#license)
  * [Cite](#cite)


## Requirements
centroAnno runs Linux and requires gcc 9.3.0+.


## Installation

```bash
git clone https://github.com/junhaiqi/centroAnno.git
cd centroAnno
make -j8  # C++11 required to compile
```
Then, there will be a binary file called centroAnno.

If necessary, you can recompile to get libspoa.a:

```bash
git clone --recursive https://github.com/rvaser/spoa.git
cd spoa
mkdir build
cd build
cmake -DCMAKE_CXX_FLAGS="-march=x86-64" .. or cmake -DCMAKE_CXX_FLAGS="-O2 -march=core2 -mtune=generic" .. (More conservative)
make -j8
cp lib/libspoa.a ../centroAnno/lib/
```

## Singularity Image

We also provide a [singularity image](https://doi.org/10.6084/m9.figshare.30138130.v1). A simple command to use it is as follows:

```bash
singularity exec /path/to/centroAnno.simg centroAnno $chr.cen.fa -o $chr.cen.denove -x anno-asm -t 20
```

## Usage

Basic command:
```bash
./centroAnno [Options:] <in.fa> -o <out_dir>
```
The specific parameters are as follows:
```bash
Version 1.0.2
Usage: ./centroAnno [Options:] <in.fa>
Options:
  -o STR     Specify the output folder [required parameters]
  -m STR     Specify the monomer template file with fasta type [default = None]
  -k INT     Specify the k-mer size [default = 11]
  -f FLOAT   Specify the fps cutoff [default = 0.6]
  -r FLOAT   Specify the repeat redio cutoff [default = 0.2]
  -w INT     Specify the window size for infering templates [default = 500000]
  -c BOOL    Specify the homopolymer compression (1: yes, 0: no) [default = 1]
  -e FLOAT   Specify the indentity cutoff for DBSCAN [default = 0.95]
  -x STR     Specify the annotated data type (anno-sat-asm: annotate centromeric alpha-satellite sequence (HiCAT/HORmon-like input), anno-asm: annotate chromosome/assembly, or anno-read: annotate sequencing reads)
[default is anno-sat-asm]
  -t INT     Specify the number of threads for template inference [default = 8]
  -M INT     Specify the maximum number of monomers that a HOR can contain [default = 50]
  -L INT     Specify the length cutoff that the annotated sequence needs to meet [default = 5000]
  -A INT     Specify the maxinum length cutoff that the annotated region in the genome needs to meet for speed [default = 1000000]
  -N INT     Specify the mininum length cutoff that the annotated region in the genome needs to meet for accuracy [default = 100]
  -F FLOAT   Specify the indentity cutoff for genome annotation [default = 0.8]
  -S BOOL    Specify the repeated sequences are scanned out without annotation (1: yes, 0: no) [default = 0]
  example command: ./centroAnno -o test test.fa
```

In fact, centroAnno's multi-threaded computing is mainly used for sequence decomposition problems (see https://academic.oup.com/bioinformatics/article/36/Supplement_1/i93/5870498). Sequence decomposition requires a lot of memory and is linearly related to the number of threads, the number of threads you set depends on the available memory. 

When we analyzed human CHM13 chromosome 1, annotation of chromosome 1 can be completed in ~1 hour, with a peak memory usage of ~3.9GB.

We have recently proposed a new low-memory algorithm for solving sequence decomposition problems, which is currently being tested, which will greatly reduce memory requirements for better parallel analysis, see [wsd](https://github.com/junhaiqi/wsd.git). 

## Example
Analyze the structure of centromere without template information from a given centromere sequence:
```bash
./centroAnno example/cen21.fa -o example/test
```

Analyze the structure of centromere using template information from a given centromere sequence:

```bash
./centroAnno example/cen21.fa -o example/test -m example/AlphaSat.fa
```

## Output
If the above example runs successfully, the folder 'test' will be the following files:

| File   | Description |
   |  :----:  | :----:  |
   | *_decomposedResult.csv  | Load the decompose results at monomer level |
   | *_horDecomposedResult.csv  | Load the decompose results at HOR level |
   | *_monomerTemplates.fa  | Load all monomers inferred by centroAnno |
   | *_HORs.fa  | Load all HORs inferred by centroAnno |

The * indicates the name of each sequence in the input fasta/fastq.gz file. When the mononer name is in the form of 12', it means the reverse complement of monomer 12. We provide a script (see `misc/misc/README.md`) to analyze and summarize the output of centroAnno.

## Acknowledgments
We thank [Meng Zhou](https://github.com/zhoudreames) for building the centroAnno v1.02 singularity image.

## License 
MIT License.

## Cite
Junhai Qi, Junchi Ma, Zheng Han, Renmin Han, Ting Yu, Guojun Li. De novo annotation of centromere with centroAnno. bioRxiv 2025.02.19.639205; doi: https://doi.org/10.1101/2025.02.19.639205
