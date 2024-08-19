
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

The * indicates the name of each sequence in the input fasta/fastq.gz file.

## Acknowledgments
None.

## License 
MIT License.

## Cite
None.
