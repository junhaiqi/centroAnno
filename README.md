## Getting Started

```bash
# Install centroAnno
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

# Analyze the structure of centromere with template information from given **sequencing reads**:
./centroAnno $your_sequencing_reads.fa -o $your_output -x anno-read
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
cmake -DCMAKE_CXX_FLAGS="-march=x86-64" ..
make -j8
cp lib/libspoa.a ../centroAnno/lib/
```


## Usage

Basic command:
```bash
./centroAnno [Options:] <in.fa>
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
  -x STR   Specify the annotated data type (anno-sat-asm: annotate centromeric alpha-satellite sequence (HiCAT/HORmon-like input), anno-asm: annotate chromosome/assembly, or anno-read: annotate sequencing reads)
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

When we analyzed human CHM13 chromosome 1, we used the following command:

```bash
./centroAnno human_chr/chr1.fasta -o CHM13_centroAnno_out/chr1 -k 10 -r 0.3 -L 10000 -G true -A 1000000
```

The log is:

```bash
        Command being timed: "./centroAnno human_chr/chr1.fasta -o CHM13_centroAnno_out/chr1 -k 10 -r 0.3 -L 10000 -G true -A 1000000"
        User time (seconds): 6015.51
        System time (seconds): 206.47
        Percent of CPU this job got: 153%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:07:46
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 4157412
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 111704253
        Voluntary context switches: 6051
        Involuntary context switches: 697541
        Swaps: 0
        File system inputs: 0
        File system outputs: 47856
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

This implies that using the default thread values, annotation of chromosome 1 can be completed in ~1 hour, with a peak memory usage of ~3.9GB.

In fact, we have recently proposed a new low-memory algorithm for solving sequence decomposition problems, which is currently being tested, which will greatly reduce memory requirements for better parallel analysis. Welcome your follow-up attention.

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

The * indicates the name of each sequence in the input fasta/fastq.gz file. When the mononer name is in the form of 12', it means the reverse complement of monomer 12.

## Acknowledgments
None.

## License 
MIT License.

## Cite
Junhai Qi, Junchi Ma, Zheng Han, Renmin Han, Ting Yu, Guojun Li. De novo annotation of centromere with centroAnno. bioRxiv 2025.02.19.639205; doi: https://doi.org/10.1101/2025.02.19.639205
