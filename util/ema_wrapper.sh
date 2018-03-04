#!/bin/bash
set -e
set -o pipefail

OUTSAM="ema_out.sam"
OUTBAM="ema_out.sorted.bam"
OUTBAM_NO_DUPS="ema_out.sorted.dupsmarked.bam"
FINALBAM="ema_out.full.sorted.bam"

SHORT=r:R:t:
LONG=ref:,rg:,threads:

if [[ -z "$EMAPATH" ]]; then
    echo "error: must specify EMAPATH in environment"
    exit 1
fi

if [[ -z "$PICARDPATH" ]]; then
    echo "error: must specify PICARDPATH in environment"
    exit 2
fi

PARSED=$(getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    exit 3
fi
eval set -- "$PARSED"

t=1

while true; do
    case "$1" in
        -r|--ref)
            r="$2"
            shift 2
            ;;
        -R|--rg)
            R="$2"
            shift 2
            ;;
        -t|--threads)
            t="$2"
            shift 2
            ;;
        -h|--help)
            echo "required options:"
            echo "-r  indexed reference FASTA path"
            echo "-R  full read group string (eg. $'@RG\tID:foo\tSM:bar' -- notice the $)"
            echo "-t  number of jobs"
            shift
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            exit 4
            ;;
    esac
done

if [[ -z "$r" ]]; then
    echo "error: specify indexed reference FASTA with -r"
    exit 5
fi

if [[ -z "$R" ]]; then
    echo "error: specify full read group string with -R (eg. $'@RG\tID:foo\tSM:bar' -- notice the $)"
    exit 6
fi

echo "Aligning..."
parallel -j "$t" --xapply "$EMAPATH align -1 {1} -2 {2} -r $r -o {1//}/$OUTSAM -R '$R'" \
                           ::: bucket0*/*1.preproc.fastq \
                           ::: bucket0*/*2.preproc.fastq
bwa mem -t "$t" -M -R "${R//$'\t'/'\t'}" "$r" bucket_no_bc/*1.no_bc.fastq bucket_no_bc/*2.no_bc.fastq > "bucket_no_bc/$OUTSAM"

echo "SAM -> sorted BAM..."
parallel -j "$t" "samtools sort -m 5G -o {//}/$OUTBAM {}" ::: bucket*/$OUTSAM
rm bucket*/*.sam

echo "Removing duplicates..."
parallel -j "$t" "java -Xmx10g -jar $PICARDPATH MarkDuplicates I={} O={//}/$OUTBAM_NO_DUPS M={//}/marked_dup_metrics.txt READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX" ::: bucket*/$OUTBAM

echo "BAM merge..."
samtools merge -@ "$t" -f -c "$FINALBAM" bucket*/$OUTBAM_NO_DUPS
rm bucket*/*.bam

echo "Indexing..."
samtools index "$FINALBAM"

echo "Done!"

