########################## Merge Reps and Generate Visualization Files #################
########################## version 1.0 #################################################
my $runinfo=$ARGV[0]; chomp $runinfo; ####### Run File
my $dirtoread=$ARGV[1]; chomp $dirtoread; ######### DATADIR
my $ScriptLocation=$ARGV[2]; chomp $ScriptLocation; ##### Script Location
my $maindir=$ARGV[3]; chomp $maindir; ################### Main Directory used for output


my $DIR=$dirtoread;

my $counter=0; my %hash_storeGSE;
open(F,"$runinfo");
while(my $data = <F>)
{
 $data =~ s/^\s+|\s+$//g;
 chomp $data;
 $counter++;
 if($counter==1)
 {}
 else
 {
  my @arr=split(/\t/,$data);
  if($arr[9] eq "NULL")
  {}
  else
  {
   my $GSE=$arr[2]."\_".$arr[1]."\_".$arr[9];
   my $store_value=$arr[1]."\t".$arr[6]."\t".$arr[7]."\t".$arr[8]."\t".$arr[9];
   $hash_storeGSE{$GSE}{$store_value}=0;
  }
 }
}
close F;
foreach my $m(sort keys %hash_storeGSE)
{
 my $createdir="Replicates_Merged";
 my $createdir1="Replicates_Merged/MACS2_RESULTS";
 if(-d $createdir)
 {
  mkdir($createdir1);
 }
 else
 {
  mkdir($createdir);
  mkdir($createdir1);
 }
}
my $RNAseqmergereps=$ScriptLocation."\/"."./RunMergeReplicates_forRNA-Seq.sh";
my $DNAseqmergereps=$ScriptLocation."\/"."./RunMergeReplicates_forChIP-Seq.sh";
my $genomesize_laevis="3.1e9";
my $genomesize_trop="1.7e9";
my $genomesize_human="2.7e9";
my $genomesize_mouse="1.87e9";
my $genomesize_zebrafish="1.2e8";
foreach my $m(sort keys %hash_storeGSE)
{
 my $dircheck=$maindir."\/".$m;
 if($dircheck eq $dirtoread)
 {
  my $countPeakCall=0;
  my @arr_split_m=split(/\_/,$m);
  my $species=$arr_split_m[1];
  print "Processesing $m for merging Replicates\n";
  foreach my $n(keys %{$hash_storeGSE{$m}})
  {
   my %hash_SRR_treatment; my %hash_SRR_control; my $header_name1; my $header_name2; my $filestomerge_treatment; my $name_of_merged_treatment; my $filestomerge_control; my $name_of_merged_control;
   my @arr=split(/\t/,$n);
   if(($arr[3] eq "NULL")||($arr[4] eq "NULL"))
   {
    print LOG "\nRead Type or Assay Type is set to NULL \.\.\. can't proceed for $n record for $m\n";
    break;
   }
   else
   {
    if(($arr[4] eq "ChIP-TF")||($arr[4] eq "ATAC")||($arr[4] eq "ChIP-Epigenetic"))
    {
     if($arr[3] eq "SINGLE-END")
     {
      print "yes\n";
      my @arr_treatment=split(/\,/,$arr[1]);
      my $numberoftreatmentReps=scalar(@arr_treatment);
      if($numberoftreatmentReps > 1)
      {
       foreach my $s1(@arr_treatment)
       {
        $hash_SRR_treatment{$s1}=0;
       }
       foreach my $x(sort keys %hash_SRR_treatment)
       {
        $filestomerge_treatment=$filestomerge_treatment." ".$x.".final.sorted.bam";
        $name_of_merged_treatment=$name_of_merged_treatment."\_".$x;
       }
       print "$filestomerge_treatment\n";
       my $revheader_name1=reverse($name_of_merged_treatment); chop $revheader_name1; my $final_mergedname=reverse($revheader_name1);
       my $revheader_name1=reverse($filestomerge_treatment); chop $revheader_name1; my $final_FilesMerge=reverse($revheader_name1);
       my $MergedBam="Replicates_Merged/".$final_mergedname;
       my $peakcalldir="Replicates_Merged/MACS2_RESULTS";
       if($species eq "xenopus-laevis")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_laevis", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_laevis", "$DIR");
        }
       }
       if($species eq "xenopus-trop")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_trop", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_trop", "$DIR");
        }
       }
       if($species eq "human")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_human", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_human", "$DIR");
        }
       }
       if($species eq "mouse")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_mouse", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_mouse", "$DIR");
        }
       }
       if(($arr[4] eq "ATAC")||($arr[4] eq "ChIP-TF"))
       {
        my $dir_to_delete=$maindir.$peakcalldir."\/".$final_mergedname."\."."broadpeaks";
        system('rm', '-rf', $dir_to_delete);
       }
       if($arr[4] eq "ChIP-Epigenetic")
       {
        my $dir_to_delete=$maindir.$peakcalldir."\/".$final_mergedname."\."."narrowpeaks";
        system('rm', '-rf', $dir_to_delete);
       }
      }
     }
     if($arr[3] eq "PAIRED-END")
     {
      my @arr_treatment=split(/\,/,$arr[1]);
      my $numberoftreatmentReps=scalar(@arr_treatment);
      if($numberoftreatmentReps > 1)
      {
       foreach my $s1(@arr_treatment)
       {
        $hash_SRR_treatment{$s1}=0;
       }
       foreach my $x(sort keys %hash_SRR_treatment)
       {
        $filestomerge_treatment=$filestomerge_treatment." ".$x."_1.final.sorted.bam";
        $name_of_merged_treatment=$name_of_merged_treatment."\_".$x;
       }
       my $revheader_name1=reverse($name_of_merged_treatment); chop $revheader_name1; my $final_mergedname=reverse($revheader_name1);
       my $revheader_name1=reverse($filestomerge_treatment); chop $revheader_name1; my $final_FilesMerge=reverse($revheader_name1);
       my $MergedBam="Replicates_Merged/".$final_mergedname;
       my $peakcalldir="Replicates_Merged/MACS2_RESULTS";
       if($species eq "xenopus-laevis")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_laevis", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_laevis", "$DIR");
        }
       }
       if($species eq "xenopus-trop")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_trop", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_trop", "$DIR");
        }
       }
       if($species eq "human")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_human", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_human", "$DIR");
        }
       }
       if($species eq "mouse")
       {
        if($arr[3] eq "SINGLE-END")
        {
         my $readtype="BAM";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_mouse", "$DIR");
        }
        if($arr[3] eq "PAIRED-END")
        {
         my $readtype="BAMPE";
         system("$DNAseqmergereps", "$final_FilesMerge", "$MergedBam", "$readtype", "$final_mergedname", "$peakcalldir", "$genomesize_mouse", "$DIR");
        }
       }
       if(($arr[4] eq "ATAC")||($arr[4] eq "ChIP-TF"))
       {
        my $dir_to_delete=$maindir.$peakcalldir."\/".$final_mergedname."\."."broadpeaks";
        system('rm', '-rf', $dir_to_delete);
       }
       if($arr[4] eq "ChIP-Epigenetic")
       {
        my $dir_to_delete=$maindir.$peakcalldir."\/".$final_mergedname."\."."narrowpeaks";
        system('rm', '-rf', $dir_to_delete);
       }
      }
     }
    }
   }
  }
 }
}