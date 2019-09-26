use File::chdir;
my $list=$ARGV[0]; chomp $list; ########### List of GSE's to process
my $dirtoread=$ARGV[1]; chomp $dirtoread; ##### Directory to read
my $ScriptLocation=$ARGV[2]; chomp $ScriptLocation; #### Location of scripts
my $RunFile=$ARGV[3]; chomp $RunFile; ################## RunFile
my $sralocation=$ARGV[4]; chomp $sralocation; ########## sralocation
my %hash_GSE;
open(F,"$list");
while(my $data = <F>)
{
 $data =~ s/^\s+|\s+$//g;
 chomp $data;
 my $path=$dirtoread."\/".$data;
 print "$path\n";
 $hash_GSE{$path}=$data;
}
close F;

my $pathtoCSBB=$ScriptLocation."\/"."CSBB-v3.0_Linux.pl";
my $pathtoExpMatandDEscript=$ScriptLocation."\/"."GenerateExpMat_and_PerformDE.pl";
my $pathtoMergeRepsRNA=$ScriptLocation."\/"."MergeReps_and_CreateBigwigs-RNA-SEQ.pl";
my $pathtoMergeRepsCHIP=$ScriptLocation."\/"."MergeReps_and_CreateBigwigs-CHIP-SEQ.pl";
foreach my $m(sort keys %hash_GSE)
{
 opendir(DIR,"$m");
 foreach my $file(readdir(DIR))
 {
  if($file =~ /.txt/)
  {
   my $filepath=$m."\/".$file;
   chdir $ScriptLocation;
   system("perl", "$pathtoCSBB", "ProcessPublicData", "$filepath", "$m", "$sralocation");
   print "Processed $file\n";
  }
 }
 closedir DIR;
 my @arr=split(/\_/,$hash_GSE{$m});
 if($arr[2] eq "RNA-SEQ")
 {
  print "CSBB has processed all the SRA Data \.\.\. Casting off steps to generate Expression matrix and perform DE analysis if Needed\n";
  chdir $ScriptLocation;
  system("perl", "$pathtoExpMatandDEscript", "$m", "$ScriptLocation");
  print "Running Merge Replicates step \.\.\. Time taking step \.\.\. please be patient\n";
  chdir $m;
  system("perl", "$pathtoMergeRepsRNA", "$RunFile", "$m", "$ScriptLocation", "$dirtoread");
 }
 if(($arr[2] eq "ChIP-TF")||($arr[2] eq "ChIP-Epigenetic")||($arr[2] eq "ATAC"))
 {
  print "CSBB has processed all the SRA Data \.\.\.\n";
  print "Running Merge Replicates step \.\.\. Time taking step \.\.\. please be patient\n";
  chdir $m;
  system("perl", "$pathtoMergeRepsCHIP", "$RunFile", "$m", "$ScriptLocation", "$dirtoread");
 }
}