#!/bin/bash
FASTQ1=$1
QUAL=$2
BOWTIE_PATH=$3
RSEM_PATH=$4
RSEM_INDEX=$5
OUTPUT=$6
FASTQC_OUT=$7
BAMTOBIGWIG=$8
CHROMSIZE=$9
GENOMESIZE=${10}
GENOMESIZE2=$( perl -e'my ($n, $x) = split("e", $ARGV[0]); print $n*10**$x;' -- $GENOMESIZE )

echo Started RSEM Run !!
$RSEM_PATH --output-genome-bam --num-threads 4 $QUAL --bowtie2 --bowtie2-path $BOWTIE_PATH $FASTQ1 $RSEM_INDEX $OUTPUT
samtools sort $OUTPUT.genome.bam -o $OUTPUT.genome.sorted.bam
samtools index $OUTPUT.genome.sorted.bam
echo Converting BAM to BigWig!!
bamCoverage --bam $OUTPUT.genome.sorted.bam -o $OUTPUT.bw --effectiveGenomeSize $GENOMESIZE2 --normalizeUsing BPM --binSize 20 --smoothLength 60 -p 4
echo Removing Temp Files !!
rm $OUTPUT.bedGraph
rm $OUTPUT.sorted.bedGraph
rm $OUTPUT.genome.bam
rm $OUTPUT.transcript.bam
echo Started Fastqc Run !!
fastqc -o $FASTQC_OUT $FASTQ1