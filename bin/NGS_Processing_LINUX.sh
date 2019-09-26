#!/bin/bash


############This script is part of NGS processing########
############Script version: v1###########################


###########Parameters#######
DIR=$1
ScriptLocation=$2
RunFile=$3
SraLocation=$4
############################


########Set Output Dir######
echo Setting Output Directory
cd $DIR
############################

######Creating RunFiles required by CSBB##
perl $ScriptLocation/Create_CSBB_RunInfo_Files_perGSE_LINUX.pl $RunFile $ScriptLocation $DIR $SraLocation
echo RunFiles generated
##########################################

##########Deploy CSBB for NGS############
perl $ScriptLocation/DeployCSBB_ForNGSProcessing_LINUX.pl $DIR/List_Of_GSEs_to_process.txt $DIR $ScriptLocation $RunFile $SraLocation
##########################################