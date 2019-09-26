#!/bin/bash
FILES=$1
OUTNAME=$2
GENOMESIZE=$3
DIR=$4
GENOMESIZE2=$( perl -e'my ($n, $x) = split("e", $ARGV[0]); print $n*10**$x;' -- $GENOMESIZE )

cd $DIR
samtools merge $OUTNAME.bam $FILES
samtools sort $OUTNAME.bam -o $OUTNAME.sorted.bam
samtools index $OUTNAME.sorted.bam
rm $OUTNAME.bam
bamCoverage --bam $OUTNAME.sorted.bam -o $OUTNAME.bw --normalizeUsing BPM --effectiveGenomeSize $GENOMESIZE2 --binSize 20 --smoothLength 60 -p 4
rm $OUTNAME.sorted.bam
rm $OUTNAME.sorted.bam.bai
