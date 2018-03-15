EMA: An aligner for barcoded short-read sequencing data
=======================================================
[![Build Status](https://travis-ci.org/arshajii/ema.svg?branch=master)](https://travis-ci.org/arshajii/ema) [![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/arshajii/ema/master/LICENSE)

EMA uses a latent variable model to align barcoded short-reads (such as those produced by [10x Genomics](https://www.10xgenomics.com)' sequencing platform). More information is available in [our paper](https://www.biorxiv.org/content/early/2017/11/16/220236). The full experimental setup is available [here](https://github.com/arshajii/ema-paper-data/blob/master/experiments.ipynb).

### Download and Compile
In a nutshell:

```
git clone --recursive https://github.com/arshajii/ema
cd ema
make
```

The `--recursive` flag is needed because EMA uses BWA's C API.

### Usage
Input FASTQs must first be preprocessed with `ema preproc` then sorted with `ema sort` (see below for more details).

```
usage: ./ema <preproc|count|align|help> [options]

preproc: preprocess barcoded FASTQ files
  -w <whitelist path>: specify whitelist [required]
  -n <num buckets>: number of barcode buckets to make [500]
  -h: apply Hamming-2 correction [off]
  -o: <output directory> specify output directory [required]
  -t <threads>: set number of threads [1]

count: performs preliminary barcode count (takes FASTQ via stdin)
  -w <whitelist path>: specify barcode whitelist [required]
  -o <output prefix>: specify output prefix [required]

align: choose best alignments based on barcodes
  -1 <FASTQ1 path>: first (preprocessed and sorted) FASTQ file [required]
  -2 <FASTQ2 path>: second (preprocessed and sorted) FASTQ file [required]
  -s <EMA-FASTQ path>: specify special FASTQ path [none]
  -x: multi-input mode; takes input files after flags and spawns a thread for each
  -r <FASTA path>: indexed reference [required]
  -o <SAM file>: output SAM file [stdout]
  -R <RG string>: full read group string (e.g. '@RG\tID:foo\tSM:bar') [none]
  -d: apply fragment read density optimization
  -p <platform>: sequencing platform (one of '10x', 'tru', 'cpt') [10x]
  -t <threads>: set number of threads [1]

help: print this help message
```

### Preprocessing (10x only)

(TODO)

### Parallelism

EMA supports several modes of parallelism:

1. Single-file parallelism: multiple worker threads work on a single FASTQ (specified with `-t <threads>`).
2. Multi-file parallelism: multiple FASTQs can be specified after the flags, and each will be assigned a thread (enabled with `-x`).
3. Multi-process parallelism: multiple EMA instances can be spawned to work on individual FASTQs.

Options 1 and 2 are advantageous because they use only a single copy of the data structures needed in alignment (i.e. those of BWA), but have a slight overhead due to multithreading.

The different options can be combined: for example, specifying `-t T` _and_ `-x` will spawn `T` threads to work on each input FASTQ.

Option 3 can be enabled by running a separate instance of EMA for each barcode bucket produced during preprocessing (this can also be used to sort the FASTQs with `ema sort`). A script using [GNU Parallel](https://www.gnu.org/software/parallel/) is provided in [`utils`](util/) and can be run as follows (in the same directory as the `bucket` folders):

```
EMAPATH=/path/to/ema/executable PICARDPATH=/path/to/picard/jar ./ema_wrapper.sh
```

### Output

EMA outputs a standard SAM file with several additional tags:

- `XG`: alignment probability
- `XC`: cloud identifier
- `XA`: alternate high-probability alignments
