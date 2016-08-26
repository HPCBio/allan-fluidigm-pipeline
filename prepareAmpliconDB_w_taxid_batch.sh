#!/bin/sh
# program:     prepare_ampliconDB_w_taxid.sh
# author:      Gloria Rendon
# date:        August, 2016
# description: a program that takes five inputs  and generates a database of amplicons. 
# The OUTPUT fasta sequences include taxonomic information
# param1: input folder, where fasta files of all amplicons are located 
# param2: output folder, where the database will be written to 
# param3: file extension to look for inside the folder, for example .fasta .faa .fa .fna
# param4: output name, string to use as prefix for all output files
# param5: email, is the email address of the recipient of the biocluster execution updates            
#########################################################################################
set -x
echo `date`

if [ $# != 5 ]
then
     echo -e "\n\nprogram $0 stopped at line $LINENO\nREASON=Parameters mismatch\n"
     echo -e "Rerun like this: $0 <inputdir> <outputdir> <FastafileExtension> <outputName> <email>\n"
     echo -e "Example: $0 ./ ./ fasta AmpliconsDB grendon@illinois.edu\n\n"
     echo -e "where:\n<inputdir> is the folder with the amplicon files \n"
     echo -e "<outputdir> is the folder where the database witll be created \n"
     echo -e "<FastafileExtension> is the file extension of the amplicon files, for example fasta\n"
     echo -e "<outputName> is the string to use as prefix on all output files\n"
     echo -e "<email> is the email address of the recipient of the biocluster execution updates\n\n\n"
     exit 1;
fi

inputdir=$1
outputdir=$2
inputExt=$3
outputPrex=$4
email=$5

set +x;
echo -e "\n##############################################################################">&2;  
echo -e "\n############   sanity check                                         ##########">&2;
echo -e "\n##############################################################################" >&2; 
set -x;
if [ ! -d $outputdir ]
then
    echo $outputdir directory not found. creating it
    mkdir -p $outputdir
else
    echo $outputdir directory exists. resetting it
    cd $outputdir
    rm ${outputPrex}*
fi

if [ ! -d $inputdir ]
then
    echo $inputdir directory not found
    exit 1;
fi

cd $inputdir
numfiles=`ls -1 *$inputExt| wc -l`

if [ $numfiles -lt 1 ]
then
    echo $inputdir directory is empty
    exit 1;
    
fi




set +x;
echo -e "\n##############################################################################" >&2;  
echo -e "\n############   defining variables for inputs and outputs     #################" >&2;
echo -e "\n##############################################################################" >&2; 
set -x; 

rootdir=/home/groups/hpcbio_shared/brianallan_group
scriptdir=$rootdir/src/  
tmpdir=$rootdir/src/tmp_logs


set +x;
echo -e "\n##############################################################################" >&2;  
echo -e "\n############   defining variables for running jobs in the cluster     ########" >&2;
echo -e "\n##############################################################################">&2; 
set -x;  

nodes=1
threads=12
queue=default
mem=4gb

set +x;
echo -e "\n##############################################################################" >&2;  
echo -e "\n############   DO NOT EDIT THE CODE BEYOND THIS POINT         ################" >&2;
echo -e "\n##############################################################################">&2; 
set -x;  
            
qsub1=$tmpdir/qsub.prepareAmpliconDB_with_taxid

echo "#PBS -S /bin/bash" > $qsub1
echo "#PBS -N prepareAmpliconDB_with_taxid" >> $qsub1
echo "#PBS -M $email" >> $qsub1
echo "#PBS -m ae" >> $qsub1
echo "#PBS -e $tmpdir/qsub.prepareAmpliconDB_with_taxid.er" >> $qsub1
echo "#PBS -o $tmpdir/qsub.prepareAmpliconDB_with_taxid.ou" >> $qsub1
echo "#PBS -l nodes=$nodes:ppn=$threads" >> $qsub1
echo "#PBS -q $queue" >> $qsub1
echo "set -x" >> $qsub1
echo -e "\n\n" >> $qsub1
echo -e "$scriptdir/PrepareAmpliconDB_using_Blast.sh $inputdir $outputdir $inputExt $outputPrex $scriptdir" >> $qsub1
`chmod g+r $qsub1 `
alnjobid=`qsub $qsub1`
echo `date`

set +x;
echo -e "\n##############################################################################"  >&2;   
echo -e "\n############   Done. Exiting now                                     #########"  >&2;
echo -e "\n##############################################################################"  >&2;
set -x;

