########### Script is part of Xenbase Public data processing pipeline ########
########### Version 1.5 (Changes to Paired End File Format and now handles technical reps and bug fixes) ######################################################
print "Create CSBB Run\-Info files per GSE\n";
print "Developed by Praneet Chaturvedi\n";
sleep(1);
########## Contact Praneet for Bugs and Issues ##############################
my $runinfo=$ARGV[0]; chomp $runinfo; ####### Run Info from Xenbase Database
my $ScriptLocation=$ARGV[1]; chomp $ScriptLocation; ########## ScriptsLocation
my $datadir=$ARGV[2]; chomp $datadir; ######################## Main Folder for writing results
my $sralocation=$ARGV[3]; chomp $sralocation; ################ PATH TO DOWNLOAD SRA files
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
open(OUT,">List_Of_GSEs_to_process.txt");
foreach my $m(sort keys %hash_storeGSE)
{
 if(-d $m)
 {
  print "\n$m folder already present and thus will use the same to write the results\n";
 }
 else
 {
  print "\nCreating $m folder to write CSBB pipeline results\n";
  mkdir($m);
  my $DEdir=$m."\/"."DE";
  my $PeakCallingDir=$m."\/"."PeakCall";
  mkdir($DEdir);
  mkdir($PeakCallingDir);
 }
 print OUT "$m\n";
}
close OUT;
open(LOG,">>LogFile.txt");

my $ScriptDownloadTechreps=$ScriptLocation."\/"."./Download_TechnicalReps.sh";
my $ScriptProcessTechreps=$ScriptLocation."\/"."./Process_TechnicalReps.sh";
my $path_to_srafile=$sralocation."\/"."sra";
foreach my $m(sort keys %hash_storeGSE)
{
 my $countPeakCall=0;
 print "Processesing $m for creating subsquent sample based RunInfo files\n";
 foreach my $n(keys %{$hash_storeGSE{$m}})
 {
  my %hash_SRR_treatment; my %hash_SRR_control; my $numberOfcontrols; my $numberOftreatments; my $header_name1; my $header_name2;
  my @arr=split(/\t/,$n);
  if(($arr[3] eq "NULL")||($arr[4] eq "NULL"))
  {
   print LOG "\nRead Type or Assay Type is set to NULL \.\.\. can't proceed for $n record for $m\n";
   break;
  }
  else
  {
   if($arr[4] eq "RNA-SEQ")
   {
    if(($arr[2] eq "NULL")||($arr[2] eq ""))
    {
     print LOG "No Control replicates/sample found for $n record of $m\n";
     my @arr_treatment=split(/\,/,$arr[1]);
     foreach my $s1(@arr_treatment)
     {
      $hash_SRR_treatment{$s1}=0;
     }
     foreach my $x(sort keys %hash_SRR_treatment)
     {
      if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
      }
      else
      {
       my $outfile_SRR=$x."\_"."RunInfo_RNA-Seq";
       open(OUT,">$m/$outfile_SRR.txt");
       print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
       print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
       close OUT;
      }
     }
    }
    else
    {
     my @arr_treatment=split(/\,/,$arr[1]);
     my @arr_control=split(/\,/,$arr[2]);
     $numberOfcontrols=scalar(@arr_control); $numberOftreatments=scalar(@arr_treatment);
     foreach my $s1(@arr_treatment)
     {
      $hash_SRR_treatment{$s1}=0;
     }
     foreach my $s1(@arr_control)
     {
      $hash_SRR_control{$s1}=0;
     }
     foreach my $s1(sort keys %hash_SRR_control)
     {
      $header_name1=$header_name1."\-".$s1;
     }
     foreach my $s1(sort keys %hash_SRR_treatment)
     {
      $header_name2=$header_name2."\-".$s1;
     }
     my $revheader_name1=reverse($header_name1); chop $revheader_name1; my $finalheader_name1=reverse($revheader_name1);
     my $revheader_name2=reverse($header_name2); chop $revheader_name2; my $finalheader_name2=reverse($revheader_name2);
     my $SampleFile_forDE=$finalheader_name1."\_"."vs"."\_".$finalheader_name2;
     if($arr[3] eq "SINGLE-END")
     {
      open(OUTDE,">$m/DE/$SampleFile_forDE.txt");
      print OUTDE "\#$SampleFile_forDE\t$numberOfcontrols\t$numberOftreatments\n";
      foreach my $x(sort keys %hash_SRR_control)
      {
       if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       print OUTDE "$x\n";
      }
       else
       {
        my $outfile_SRR=$x."\_"."RunInfo_RNA-Seq";
        open(OUT,">$m/$outfile_SRR.txt");
        print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
        print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
        close OUT;
        print OUTDE "$x\n";
       }
      }
      foreach my $x(sort keys %hash_SRR_treatment)
      {
       if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       print OUTDE "$x\n";
      }
       else
       {
        my $outfile_SRR=$x."\_"."RunInfo_RNA-Seq";
        open(OUT,">$m/$outfile_SRR.txt");
        print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
        print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
        close OUT;
        print OUTDE "$x\n";
       }
      }
      close OUTDE;
     }
     if($arr[3] eq "PAIRED-END")
     {
      open(OUTDE,">$m/DE/$SampleFile_forDE.txt");
      print OUTDE "\#$SampleFile_forDE\t$numberOfcontrols\t$numberOftreatments\n";
      foreach my $x(sort keys %hash_SRR_control)
      {
       if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
        my $pairedcontrol=$x."_1";
        print OUTDE "$pairedcontrol\n";
      }
       else
       {
        my $outfile_SRR=$x."\_"."RunInfo_RNA-Seq";
        open(OUT,">$m/$outfile_SRR.txt");
        print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
        print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
        close OUT;
        my $pairedcontrol=$x."_1";
        print OUTDE "$pairedcontrol\n";
       }
      }
      foreach my $x(sort keys %hash_SRR_treatment)
      {
       if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       my $pairedtreatment=$x."_1";
        print OUTDE "$pairedtreatment\n";
      }
       else
       {
        my $outfile_SRR=$x."\_"."RunInfo_RNA-Seq";
        open(OUT,">$m/$outfile_SRR.txt");
        print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
        print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
        close OUT;
        my $pairedtreatment=$x."_1";
        print OUTDE "$pairedtreatment\n";
       }
      }
      close OUTDE;
     }
    }
   }
   if(($arr[4] eq "ATAC")||($arr[4] eq "ChIP-TF")||($arr[4] eq "ChIP-Epigenetic"))
   {
    if(($arr[2] eq "NULL")||($arr[2] eq ""))
    {
     print LOG "No Input replicates/sample found for $n record of $m\n";
     my @arr_treatment=split(/\,/,$arr[1]);
     foreach my $s1(@arr_treatment)
     {
      $hash_SRR_treatment{$s1}=0;
     }
     foreach my $x(sort keys %hash_SRR_treatment)
     {
      if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
      }
      else
      {
       my $outfile_SRR=$x."\_"."RunInfo_DNA-Seq";
       open(OUT,">$m/$outfile_SRR.txt");
       print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
       print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
       close OUT;
      }
     }
    }
    else
    {
     my @arr_treatment=split(/\,/,$arr[1]);
     my @arr_control=split(/\,/,$arr[2]);
     $numberOfcontrols=scalar(@arr_control); $numberOftreatments=scalar(@arr_treatment);
     foreach my $s1(@arr_treatment)
     {
      $hash_SRR_treatment{$s1}=0;
     }
     foreach my $s1(@arr_control)
     {
      $hash_SRR_control{$s1}=0;
     }
     $countPeakCall++;
     my $SampleFile_forPeakCall=$m."\_"."PeakCall_Set".$countPeakCall;
     open(OUTDE,">$m/PeakCall/$SampleFile_forPeakCall.txt");
     print OUTDE "\#$SampleFile_forPeakCall\t$numberOfcontrols\t$numberOftreatments\n";
     foreach my $x(sort keys %hash_SRR_control)
     {
      if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       print OUTDE "$x\n";
      }
      else
      {
       my $outfile_SRR=$x."\_"."RunInfo_DNA-Seq";
       open(OUT,">$m/$outfile_SRR.txt");
       print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
       print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
       close OUT;
       print OUTDE "$x\n";
      }
     }
     foreach my $x(sort keys %hash_SRR_treatment)
     {
      if($x =~ /_/)
      {
       print "Technical Replicates found \.\.\. Please see run time will increase\n";
       sleep(1);
       print "\nProcess Update \: Download, Convert to fastq, combine and process\n";
       my @arr_treatment_technical=split(/\_/,$x);
       print "Setting up a new folder for technical replicates set\: $x\n";
       my $technicalRepFolder=$datadir."\/".$m."\/".$x;
       if(-d $technicalRepFolder)
       {}
       else
       {
        mkdir($technicalRepFolder);
       }
       foreach my $techrep(@arr_treatment_technical)
       {
        my $pathtostore=$technicalRepFolder;
        my $path_to_sratoolkit_prefetch=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./prefetch";
        my $path_to_sratoolkit_fastqdump=$ScriptLocation."\/"."Modules/sratoolkit.2.9.0-mac64/bin/./fastq-dump";   
        system("$ScriptDownloadTechreps", "$techrep", "$pathtostore", "$path_to_sratoolkit_prefetch", "$path_to_sratoolkit_fastqdump", "$path_to_srafile");
       }
       if($arr[3] eq "SINGLE-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x;
        my $combinedtechfile2="NA";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       if($arr[3] eq "PAIRED-END")
       {
        my $pathtostore_new=$technicalRepFolder;
        my $path_to_gsefolder=$datadir."\/".$m;
        my $combinedtechfile1=$x."_1";
        my $combinedtechfile2=$x."_2";
        system("$ScriptProcessTechreps", "$arr[4]", "$arr[3]", "$pathtostore_new", "$combinedtechfile1", "$combinedtechfile2", "$arr[0]", "$path_to_gsefolder", "$ScriptLocation");
       }
       print OUTDE "$x\n";
      }
      else
      {
       my $outfile_SRR=$x."\_"."RunInfo_DNA-Seq";
       open(OUT,">$m/$outfile_SRR.txt");
       print OUT "SRRid\tSpecies\tAssay\tReadType\tSampleName\n";
       print OUT "$x\t$arr[0]\t$arr[4]\t$arr[3]\t$x\n";
       close OUT;
       print OUTDE "$x\n";
      }
     }
     close OUTDE;
    }
   }
  }
 }  
}
close LOG;