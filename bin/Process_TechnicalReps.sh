#!/bin/bash
echo "Merge Technical reps Fastq and process"


SEQ=$1
TYPE=$2
DATADIR=$3
FILENAME=$4
FILENAME1=$5
SPECIES=$6
GSEFOLDER=$7
SCRIPTLOCATION=$8



if [ $SEQ == 'RNA-SEQ' ] && [ $TYPE == 'SINGLE-END' ]; then
	cd $DATADIR
	cat *.fastq >> $FILENAME.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-RNASeq_SingleEnd $DATADIR/$FILENAME.fastq $SPECIES $DATADIR phred33
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.isoforms.results $GSEFOLDER
	mv *.bai $GSEFOLDER
fi

if [ $SEQ == 'RNA-SEQ' ] && [ $TYPE == 'PAIRED-END' ]; then
	cd $DATADIR
	cat *_1.fastq >> $FILENAME.fastq
	cat *_2.fastq >> $FILENAME1.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-RNASeq_PairedEnd $DATADIR/$FILENAME.fastq $DATADIR/$FILENAME1.fastq $SPECIES $DATADIR phred33
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.isoforms.results $GSEFOLDER
	mv *.bai $GSEFOLDER
fi

if [ $SEQ == 'ChIP-TF' ] && [ $TYPE == 'SINGLE-END' ]; then
    cd $DATADIR
	cat *.fastq >> $FILENAME.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_SingleEnd $DATADIR/$FILENAME.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi

if [ $SEQ == 'ATAC' ] && [ $TYPE == 'SINGLE-END' ]; then
    cd $DATADIR
	cat *.fastq >> $FILENAME.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_SingleEnd $DATADIR/$FILENAME.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi

if [ $SEQ == 'ChIP-Epigenetic' ] && [ $TYPE == 'SINGLE-END' ]; then
    cd $DATADIR
	cat *.fastq >> $FILENAME.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_SingleEnd $DATADIR/$FILENAME.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi

if [ $SEQ == 'ChIP-TF' ] && [ $TYPE == 'PAIRED-END' ]]; then
	cd $DATADIR
	cat *_1.fastq >> $FILENAME.fastq
	cat *_2.fastq >> $FILENAME1.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_PairedEnd $DATADIR/$FILENAME.fastq $DATADIR/$FILENAME1.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi

if [ $SEQ == 'ATAC' ] && [ $TYPE == 'PAIRED-END' ]]; then
	cd $DATADIR
	cat *_1.fastq >> $FILENAME.fastq
	cat *_2.fastq >> $FILENAME1.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_PairedEnd $DATADIR/$FILENAME.fastq $DATADIR/$FILENAME1.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi

if [ $SEQ == 'ChIP-Epigenetic' ] && [ $TYPE == 'PAIRED-END' ]]; then
	cd $DATADIR
	cat *_1.fastq >> $FILENAME.fastq
	cat *_2.fastq >> $FILENAME1.fastq
	cd $SCRIPTLOCATION
	perl CSBB-v3.0_MacOS.pl Process-ChIP-ATAC_PairedEnd $DATADIR/$FILENAME.fastq $DATADIR/$FILENAME1.fastq $SPECIES $DATADIR phred33 $SEQ
	cd $DATADIR
	mv *.bam $GSEFOLDER
	mv *.bw $GSEFOLDER
	mv *.bai $GSEFOLDER
	mkdir MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS MACS2_RESULTS.$FILENAME
	mv MACS2_RESULTS.$FILENAME $GSEFOLDER
fi