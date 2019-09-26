#!/bin/bash
FASTQ1=$1
FASTQ2=$2
QUAL=$3
BOWTIE_PATH=$4
RSEM_PATH=$5
RSEM_INDEX=$6
OUTPUT=$7
FASTQC_OUT=$8
BAMTOBIGWIG=$9
CHROMSIZE=${10}
GENOMESIZE=${11}
GENOMESIZE2=$( perl -e'my ($n, $x) = split("e", $ARGV[0]); print $n*10**$x;' -- $GENOMESIZE )

echo Started RSEM Run !!
$RSEM_PATH --output-genome-bam --paired-end --num-threads 4 $QUAL --bowtie2 --bowtie2-path $BOWTIE_PATH $FASTQ1 $FASTQ2 $RSEM_INDEX $OUTPUT
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
fastqc -o $FASTQC_OUT $FASTQ2