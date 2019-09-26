######Copyright Xenbase###########################
######General Public License v3###################
##############Script version:v1###################


############### Please Contact Praneet Chaturvedi or Joshua Fortriede on Github for bugs or issues ##############


############### Parameters ##############
my $RunFile = $ARGV[0]; chomp $RunFile;           ############### RunFile specifying samples to process
my $DirectoryOut = $ARGV[1]; chomp $DirectoryOut; ############### Output Directory
my $OS = $ARGV[2]; chomp $OS;                     ############### OS type
my $srapath=$ARGV[3]; chomp $srapath;             ############### path to SRA Folder to download
#########################################


############## Information ###############
print "\nXenbase GEO NGS Analysis Pipeline \-\-version 1.0\n";
print "\nPlease cite Xenbase and NGS pipline when using for your research\n";
print "\nPlease check back every month for new releases and bug fixes\n";
###########################################

############### Creating Output Directory #################
if(-d $DirectoryOut)
{
 print "\n$DirectoryOut is already present \.\. Pipeline will use the same to write the results\n";
 sleep(1);
}
else
{
 mkdir($DirectoryOut);
 print "\nPipeline has created $DirectoryOut for writing analysis results\n";
}
###########################################################


############## Get Processing Scripts Location ############
use Cwd;
my $cwd = getcwd;
my $ScriptsLocation=$cwd."\/"."bin";
###########################################################


############# NGS Processing #############################

if($OS eq "UNIX")
{
 my $NGSshell=$ScriptsLocation."\/"."./NGS_Processing_UNIX.sh";
 print "$NGSshell\n";
 system("$NGSshell", "$DirectoryOut", "$ScriptsLocation", "$RunFile", "$srapath");
}

if($OS eq "LINUX")
{
 my $NGSshell=$ScriptsLocation."\/"."./NGS_Processing_LINUX.sh";
 system("$NGSshell", "$DirectoryOut", "$ScriptsLocation", "$RunFile", "$srapath");
}
###########################################################