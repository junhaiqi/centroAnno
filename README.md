## Getting Started

```bash
# Install centroAnno
git clone https://github.com/junhaiqi/centroAnno.git
cd centroAnno && make -j8  # C++11 required to compile

# Run on test data
./centroAnno example/cen21.fa -o test

# Analyze the tandem repeats/HORs from a given **genome/assembly/general sequence**:
./centroAnno $your_genome.fa -o $your_output -G true

# An example command for analyzing the human chromosome 1:
./centroAnno $your_chr1.fasta -o $your_chr1_out -k 10 -r 0.3 -L 10000 -G true -A 1000000
```

# Analyze the structure of centromere without template information from a given **centromeric alpha-satellite sequence/centromeric alpha-satellite assembly**:
./centroAnno $your_centromere.fa -o $your_output

# Analyze the structure of centromere with template information from a given **centromeric alpha-satellite sequence/centromeric alpha-satellite assembly**:
./centroAnno $your_centromere.fa -m $your_templates.fa -o $your_output


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
   | Parameters   | Description |
   |  :----:  | :----:  |
   | -o  | Specify the output folder [required parameters] |
   | -m  | Specify the monomer template file with fasta type [default = None] |
   | -k  | Specify the k-mer size [default = 13] |
   | -f  | Specify the fps cutoff [default = 0.6] |
   | -r  | Specify the repeat redio cutoff [default = 0.2] |
   | -w  | Specify the window size for infering templates [default = 500000] |
   | -c  | Specify closing the homopolymer compression [default = false] |
   | -e  | Specify the indentity cutoff for DBSCAN [default = 0.95] |
   | -t  | Specify the number of threads for template inference [default = 8] |
   | -M  | Specify the maximum number of monomers that a HOR can contain [default = 50] |
   | -L  | Specify the length cutoff that the annotated sequence needs to meet [default = 5000] |
   | -A  | Specify the maxinum length cutoff that the annotated region in the genome needs to meet for speed [default = 100000] |
   | -N  | Specify the mininum length cutoff that the annotated region in the genome needs to meet for accuracy [default = 100] |
   | -F  | Specify the indentity cutoff for genome annotation [default = 0.8] |
   | -S  | Specify the repeated sequences are scanned out without annotation [default = false] |
   | -G  | Specify the tendem repeat of genome/assembly/read are annotated [default = false] |

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

Analyze the tendem repeat unit and HOR from a given genome/assembly/long sequencing read:

```bash
./centroAnno example/simulated_genome.fasta -o example/test -G true
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
