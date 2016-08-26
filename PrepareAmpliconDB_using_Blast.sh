#!/bin/sh
# program:     prepareAmpliconDB_using_Blast.sh
# author:      Gloria Rendon
# date:        August, 2016
# description: a program that takes five inputs  and generates a database of amplicons. 
# The OUTPUT fasta sequences include taxonomic information
#param1: input folder, where fasta files of all amplicons are located 
# param2: output folder, where the database will be written to 
# param3: file extension to look for inside the folder, for example .fasta .faa .fa .fna
# param4: output name, string to use as prefix for all output files
# param5: scriptdir, location of  all scripts needed to run this program          
#########################################################################################
set -x
echo `date`

if [ $# != 5 ]
then
     echo -e "\n\nprogram $0 stopped at line $LINENO\nREASON=Parameters mismatch\n"
     echo -e "Rerun like this: $0 <inputdir> <outputdir> <FastafileExtension> <outputName> <scriptdir>\n"
     echo -e "Example: $0 ./ ./ fasta AmpliconsDB /home/groups/hpcbio_shared/brianallan_group/src\n\n"
     echo -e "where:\n<inputdir> is the folder with the amplicon files \n"
     echo -e "<outputdir> is the folder where the database witll be created \n"
     echo -e "<FastafileExtension> is the file extension of the amplicon files, for example fasta\n"
     echo -e "<outputName> is the string to use as prefix on all output files\n"
     echo -e "<scriptdir> is the folder where all the scripts of the pipeline are located\n\n\n"
     exit 1;
fi


inputdir=$1
outputdir=$2
inExt=$3
prefix=$4
scriptdir=$5

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining variables and filenames               ################"
echo -e "\n##############################################################################">&2; set -x;  



outBlastfile=${prefix}_vs_NT_Blast.txt
newfastadbname=${prefix}_tmp.fa
newtablename=${prefix}_tmp.tsv
outfastadbname=${prefix}.fa
outtablename=${prefix}.tsv
redotablename=${prefix}_without_taxids_REDO.txt


set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining the tools                             ################"
echo -e "\n##############################################################################">&2; set -x;  

blastdb=/home/mirrors/NCBI/BLAST_DBS/20160717/nt
module load blast+/2.2.31

set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   Main Program starts here                      #################"
echo -e "\n##############################################################################">&2; set -x;  

if [ ! -d $scriptdir ]
then
    echo $scriptdir directory not found
    exit 1;
fi

if [ ! -d $outputdir ]
then
    echo $outputdir directory not found
    exit 1;
fi

if [ ! -d $inputdir ]
then
    echo $inputdir directory not found
    exit 1;
fi



set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   run   generateDB.pl                           #################"
echo -e "\n##############################################################################">&2; set -x;  

# this script generates two files $newfastadbname and $newtablename
# the fasta file has all the fasta sequences minus duplicates, a unique sequence identifier is assigned to each sequence
# the table file has details for each sequence such as the list of files that contain this sequence
# we do not assume that the input fasta sequences have headers with complete information such as accession and taxid

cd $inputdir

perl $scriptdir/generateDB.pl $inExt $newfastadbname  $newtablename
echo `date`

echo check that files were created

if [ ! -s $newfastadbname ]
then
     echo $newfastadbname file not created. generateDB.pl failed.
     exit 1
fi
if [ ! -s $newtablename ]
then
     echo $newtablename file not created. generateDB.pl failed.
     exit 1
fi

#sleep 5s

set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   run   blastn                                  #################"
echo -e "\n##############################################################################">&2; set -x; 

# this blast command generates a blast report in table format named $outBlastfile
# this step is run to validate the information in the header of the sequences AND to obtain taxonomic information

blastn -task megablast -db $blastdb -query $newfastadbname -num_alignments 100 -num_threads 10 -max_hsps 1 -out $outBlastfile -outfmt "6 std qlen slen  sacc staxids sscinames"

echo `date`

echo check that files were created

if [ ! -s $outBlastfile ]
then
     echo $outBlastfile file not created. blastn failed.
     exit 1
fi

#sleep 5s

set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   run   ExpandDB.pl                             #################"
echo -e "\n##############################################################################">&2; set -x;  
# this script takes the results of the previous two commands, merges the information and generates two files $xtrafastadbname $xtratablename
# the fasta file has sequences with the right taxonomic information in the header
# the table file also has the taxonomic information added 

perl $scriptdir/ExpandDB.pl $outBlastfile $newfastadbname $newtablename $outfastadbname $outtablename $redotablename

echo `date`

echo check that files were created

if [ ! -s $outfastadbname ]
then
     echo $xoutfastadbname file not created. ExpandDB.pl failed.
     exit 1
fi

if [ ! -s $outtablename ]
then
     echo $outtablename file not created. ExpandDB.pl failed.
     exit 1
fi

set +x;
echo -e "\n##############################################################################">&2;      
echo -e "\n############   move results to outputdir and cleanup         #################">&2;
echo -e "\n##############################################################################">&2; 
set -x;  

rm $outBlastfile 
rm $newfastadbname 
rm $newtablename 
mv $outfastadbname $outputdir/$outfastadbname
mv $outtablename   $outputdir/$outtablename
mv $redotablename   $outputdir/$redotablename
echo `date`

set +x;
echo -e "\n##############################################################################">&2;      
echo -e "\n############   done. exiting now                             #################">&2;
echo -e "\n##############################################################################">&2; 
set -x;  

