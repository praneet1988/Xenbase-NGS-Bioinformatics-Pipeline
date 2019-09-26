#!/bin/bash


SAMPLE=$1
DATADIR=$2
PREFETCH=$3
FASTQDUMP=$4
SRAPATH=$5

$PREFETCH $SAMPLE
mv $SRAPATH/$SAMPLE.sra $DATADIR
$FASTQDUMP --split-3 -O $DATADIR $DATADIR/$SAMPLE.sra
echo "Downloaded Technical Rep $SAMPLE from SRA and converted to fastq"