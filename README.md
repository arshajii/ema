EMA: Expectation Maximization Aligner
-------------------------------------

EMA leverages an existing sensitive all-mapper to align barcoded short reads (such as those produced by 10X Genomics' sequencing platform) through expectation maximization.

### Build
Requires BWA library and GNU C99 or later. Can be compiled with `make`.

### Usage
Input FASTQs must first be preprocessed with asdasda 

```
usage: ./ema <preproc|sort|count|align|help> [options]

preproc: preprocess barcoded FASTQ files
  -1 <FASTQ1 path>: first FASTQ file [required]
  -2 <FASTQ2 path>: second FASTQ file [none]
  -w <whitelist path>: barcode whitelist [required]
  -n <num buckets>: number of barcode buckets to make [20]
  -c <counts file>: preexisting barcode counts [none]

sort: sort preprocessed FASTQs by barcode
  -1 <FASTQ1 path>: first FASTQ file [required]
  -2 <FASTQ2 path>: second FASTQ file [required]

count: performs preliminary barcode count
  -1 <FASTQ1 path>: first FASTQ file [required]
  -w <whitelist path>: barcode whitelist [required]
  -i: indicates FASTQ is interleaved
  -o <output file>: output file [stdout]

align: choose best alignments based on barcodes
  -1 <FASTQ1 path>: first (preprocessed and sorted) FASTQ file [required]
  -2 <FASTQ2 path>: second (preprocessed and sorted) FASTQ file [required]
  -r <FASTA path>: indexed reference [required]
  -o <SAM file>: output SAM file [default: stdout]
  -R <RG string>: full read group string (e.g. $'@RG\tID:foo\tSM:bar') [default: none]

help: print this help message
```

