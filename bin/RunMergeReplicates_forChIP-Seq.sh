#!/bin/bash
FILES=$1
OUTNAME=$2
READTYPE=$3
FILENAME=$4
PEAKCALLDIR=$5
GENOMESIZE=$6
DIR=$7
GENOMESIZE2=$( perl -e'my ($n, $x) = split("e", $ARGV[0]); print $n*10**$x;' -- $GENOMESIZE )

cd $DIR
samtools merge $OUTNAME.bam $FILES
samtools sort $OUTNAME.bam -o $OUTNAME.sorted.bam
samtools index $OUTNAME.sorted.bam
samtools markdup -r $OUTNAME.sorted.bam $OUTNAME.final.bam
samtools sort $OUTNAME.final.bam -o $OUTNAME.final.sorted.bam
samtools index $OUTNAME.final.sorted.bam
rm $OUTNAME.final.bam
rm $OUTNAME.sorted.bam
rm $OUTNAME.sorted.bam.bai
rm $OUTNAME.bam
bamCoverage --bam $OUTNAME.final.sorted.bam -o $OUTNAME.bw --binSize 20 --normalizeUsing RPGC --effectiveGenomeSize $GENOMESIZE2 --extendReads 200 --smoothLength 60 -p 4
macs2 callpeak -t $OUTNAME.final.sorted.bam --buffer-size 2000 -f $READTYPE -g $GENOMESIZE -n $FILENAME.narrowpeaks --call-summits --outdir $PEAKCALLDIR/$FILENAME.narrowpeaks
macs2 callpeak -t $OUTNAME.final.sorted.bam --buffer-size 2000 -f $READTYPE -g $GENOMESIZE -n $FILENAME.broadpeaks --broad --broad-cutoff 0.1 --outdir $PEAKCALLDIR/$FILENAME.broadpeaks
rm $OUTNAME.final.sorted.bam
rm $OUTNAME.final.sorted.bam.bai