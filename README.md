EMA: Expectation Maximization Aligner
-------------------------------------

EMA leverages an existing sensitive all-mapper to align barcoded short reads (such as those produced by 10X Genomics' sequencing platform) through expectation maximization.

### Build
A simple `make` should suffice. Requires GNU C99 or later.

### Usage
```
usage: ./ema <preproc|align|help> [options]

preproc: preprocess barcoded FASTQ files
  -1 <fastq1 path>: specify first FASTQ file [required]
  -2 <fastq2 path>: specify second FASTQ file [none]
  -w <whitelist path>: specify whitelist [required]
  -n <num buckets>: number of barcode buckets to make [20]
  -c <counts file>: specify preexisting barcode counts [none]

count: performs preliminary barcode count
  -1 <fastq1 path>: specify first FASTQ file [required]
  -w <whitelist path>: specify whitelist [required]
  -i: indicates FASTQ is interleaved
  -o <output file>: specify output file [stdout]

align: choose best alignments based on barcodes
  -s <SAM file>: multi-mappings in SAM format [required]
  -i <fai file>: fai file for reference used in mapping [required]
  -o <SAM file>: output SAM file [default: stdout]

help: print this help message
```

