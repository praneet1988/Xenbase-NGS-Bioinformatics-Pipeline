################ Generate Expression Matrices and Perform DE analysis ##############
################ version 1.1 (Changed Filtering and added capability to perform DE with no reps) #######################################################
print "\nGenerate Expression Matrices and Perform DE analysis\n";
print "Developed by Praneet Chaturvedi\n";
sleep(1);

my $dir_to_read=$ARGV[0]; chomp $dir_to_read;    ######### DATADIR
my $ScriptLocation=$ARGV[1]; chomp $ScriptLocation; ###### Script Locations
my $pathtoCSBB=$ScriptLocation."\/"."CSBB-v3.0_MacOS.pl";

my @species_temp=split(/\//,$dir_to_read);
my $species_temp1=pop(@species_temp);
my @species_temp2=split(/\_/,$species_temp1);
my $species=$species_temp2[1];
print "\nObtained all the RNA\-SEQ based processed GSE dir\'s\n";
print "\nDeploying CSBB to generate TPM and Counts Matrix per GSE\n";
sleep(1.5);
system("perl", "$pathtoCSBB","Generate-TPM-Counts-Matrix","$dir_to_read","$species","$dir_to_read");
print "\nGenerating Genome Wide Matrices completed \.\. proceeding to DE analysis\n";

my $dir=$dir_to_read."\/"."FilesForDE";
if(-d $dir)
{}
else
{
 print "Creating directory $dir\n";
 mkdir($dir);
}
my $dir_to_open=$dir_to_read."\/"."DE";
opendir(F,"$dir_to_open");
foreach my $f(readdir(F))
{
 if($f=~/.txt/)
 {
  open(R,"$dir_to_open/$f");
  print "\nReading $f to create specific set based Counts Matrix for DE analysis\n";
  my $counter=0; my @samples; my $numberofcontrol; my $numberoftreatment;
  while(my $data = <R>)
  {
   $data =~ s/^\s+|\s+$//g;
   chomp $data;
   $counter++;
   if($counter == 1)
   {
    my @arr_split=split(/\t/,$data);
    $numberofcontrol=$arr_split[1]; $numberoftreatment=$arr_split[2];
   }
   else
   {
    push(@samples,$data);
   }
  }
  close R;
  print "Reading $dir_to_read\/Genes_Counts_Matrix.txt\n";
  open(R,"$dir_to_read/Genes_Counts_Matrix.txt");
  my $counter1=0; my @header_file; my %hash_counts;
  while(my $data = <R>)
  {
   chomp $data;
   $counter1++;
   if($counter1==1)
   {
    my @arr_split=split(/\t/,$data,2);
    @header_file=split(/\t/,$arr_split[1]);
   }
   else
   {
    my @arr_split=split(/\t/,$data,2);
    my @exp=split(/\t/,$arr_split[1]);
    for(my $i=0;$i<@exp;$i++)
    {
     my $value=$header_file[$i]."\t".$exp[$i];
     my $key=$arr_split[0];
     $hash_counts{$key}{$value}=0;
    }
   }
  }
  close R;
  my $header_final;
  foreach my $x(@samples)
  {
   $header_final=$header_final."\t".$x;
  }
  open(OUTFILE,">$dir_to_read/FilesForDE/$f");
  print OUTFILE "Gene$header_final\n";
  foreach my $x(sort keys %hash_counts)
  {
   my $finalcounts="";
   foreach my $a(@samples)
   {
    my $value_to_store="";
    foreach my $y(keys %{$hash_counts{$x}})
    {
     my @arr=split(/\t/,$y);
     if($arr[0] eq $a)
     {
      $value_to_store=$arr[1];
     }
    }
    if($value_to_store eq "")
    {
     $finalcounts=$finalcounts."\t"."NA";
    }
    else
    {
     $finalcounts=$finalcounts."\t".$value_to_store;
    }
   }
   print OUTFILE "$x$finalcounts\n";
  }
  close OUTFILE;
  print "Generated Specific Counts File in $dir_to_read\/FilesForDE\/ for $f\n";
  my $path_to_counts=$dir_to_read."\/"."FilesForDE"."\/".$f;
  print "\nDeploying CSBB to run Differential Expression on $dir_to_read\/FilesForDE\/$f\n";
  if(($numberofcontrol == 1)&&($numberoftreatment == 1))
  {
   my $dirtoRnew=$ScriptLocation."\/"."RUVseq_NoReps.r"; my $pathforR=$dir_to_read."\/"."FilesForDE"; my $filenameforR=$f; my $countstofilter_new=0; my $samplestofilter_new=0;
   system(Rscript,$dirtoRnew,$pathforR,$filenameforR,$numberofcontrol,$numberoftreatment,$countstofilter_new,$samplestofilter_new);
   print "\nR Run complete \.\.\. Now writing header with perl\n";
   my $Rfile = $pathforR."\/"."temporaryfile"."\."."txt";
   my $outfile_new = "DE_Using_UQ_Normalization"."\_".$filenameforR;
   open(OUTDE,">$pathforR/$outfile_new");
   print OUTDE "Gene\tlogFC\tlogCPM\tLR\tPvalue\tFDR\n";
   my $counter1=0;
   open(V,"$Rfile");
   while(my $data = <V>)
   {
    $data =~ s/^\s+|\s+$//g;
    chomp $data;
    print "$data\n";
    $counter1++;
    if($counter1==1)
    {}
    else
    {
     print OUTDE "$data\n";
    }
   }
   close V;
   close OUTDE;
   unlink $Rfile;
  }
  else
  {
   system("perl", "$pathtoCSBB", "DifferentialExpression", "$path_to_counts", "$numberofcontrol", "$numberoftreatment", "0", "0", "UpperQuantile");
  }
 }
}
closedir F;