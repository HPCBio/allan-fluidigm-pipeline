#!/bin/sh
# program:     preprocess_reads.sh
# author:      Gloria Rendon
# date:        April, 2017
# description: a program that generates a configuration file specified in the line INPUTFILES of the runfile
#              INPUTFILES is a text file with fluidigm-demultiplexed paired raw reads. If the file exists, it will be reset.
#              The line GENERATE_INPUTS of the configuration file must be set of YES in order to run this script 
#              
# IMPORTANT:   run this script from a COMPUTE NODE
#########################################################################################

runfile=$1        #full path to configuration file

set -x
echo `date`

if [ $# != 1 ]
then
	echo -e "\n\nprogram $0 stopped at line $LINENO\nREASON=Parameters mismatch\n\n\nRerun like this: $0 <runfile>\n\n"
	echo -e "where <runfile> is the configuration file for this run\n"
	exit 1;
fi
if [ `expr ${#runfile}` -lt 1 ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$runfile file not specified. Exiting now."
	exit 1;
fi
if [ ! -s $runfile ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$runfile empty file. Exiting now."
	exit 1;
fi

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining variables for inputs and outputs     #################"
echo -e "\n##############################################################################" >&2; set -x; 
##################################################################################################
############   INPUTFILES, is a text file. One sample per row, columns are separated by <SPACE>
############   col1=sample col2=barcode col3-read1 col4=read2
############   regexpr for read names: $sample-$barcode_<sequencingindex>_R{1,2}.fastq
##################################################################################################

rootdir=$( cat $runfile | grep "PROJECTDIR" | cut -d "=" -f2 )
rawreadsdir=$( cat $runfile | grep "RAW_READS_DIR" | cut -d "=" -f2 )
infiles=$( cat $runfile | grep "INPUTFILES" | cut -d "=" -f2 )                 #raw reads to process
samplefile=$( cat $runfile | grep "SAMPLEFILE" | cut -d "=" -f2 )              #one samplename per line
barcodefile=$( cat $runfile | grep "BARCODEFILE" | cut -d "=" -f2 )            #one barcodename per line
generateInput=$( cat $runfile | grep "GENERATE_INPUTS" | cut -d "=" -f2 )

set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   sanity  check                                 #################"
echo -e "\n##############################################################################">&2; set -x;  

if [ `expr ${#generateInput}` -lt 1 ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=GENERATE_INPUTS value not specified. Exiting now."
	exit 1;
fi

if [ $generateInput != "YES" -a  $generateInput != "Y" ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=GENERATE_INPUTS invalid value. Exiting now."
	exit 1;
fi

if [ ! -d $rootdir ]
then
	echo -e "$0 stopped at line $LINENO\sREASON=$rootdir directory not found"
	exit 1;
fi

if [ ! -d $rawreadsdir ]
then
	echo -e "$0 stopped at line $LINENO\sREASON=$rawreadsdir directory not found"
	exit 1;
fi



if [ `expr ${#barcodefile}` -lt 1 ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$barcodefile file not specified. Exiting now."
	exit 1;
fi
if [ ! -s $barcodefile ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$barcodefile empty file. Exiting now."
	exit 1;
fi

if [ `expr ${#samplefile}` -lt 1 ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$samplefile file not specified. Exiting now."
	exit 1;
fi

if [ ! -s $samplefile ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$samplefile empty file. Exiting now."
	exit 1;
fi

if [ `expr ${#infiles}` -lt 1 ]
then
	echo -e "$0 stopped at line $LINENO\nREASON=$infiles file not specified. Exiting now."
	exit 1;
fi

if [ -s $infiles ]
then
	echo -e "$infiles file exists. Resetting now."
	truncate -s 0 $infiles
fi



set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   parameter checked. Everything ok              #################"
echo -e "\n##############################################################################">&2; set -x;


echo -e "\n##############################################################################"    
echo -e "\n############   analysis loop starts  here.                           #########" 
echo -e "\n##############################################################################"  

listBarcodes=$( cat $barcodefile ) 

cd $rawreadsdir


while read sample
do

	echo -e "\n############  start processing $sample.                      #########" 

	if [ `expr ${#sample}` -lt 1 ]
	then
		echo -e "########## skip empty line           ###########"
		continue
	fi

	if [ ! -d $rawreadsdir/$sample ]
	then
		echo -e "########## $rawreadsdir/$sample directory not found. skip sample=$sample"
		continue
	fi
	
	for barcode in $listBarcodes
	do
		echo -e "########## generate line for sample=$sample barcode=$barcode ###########"
		
		R1=$( ls $rawreadsdir/${sample}/${sample}-${barcode}_*_R1.fastq )
		R2=$( ls $rawreadsdir/${sample}/${sample}-${barcode}_*_R2.fastq )
		
		if [ `expr ${#R1}` -lt 1 ]
		then
			echo -e "########## no R1 found for sample=$sample barcode=$barcode ###########"
			continue
		fi
		if [ `expr ${#R2}` -lt 1 ]
		then
			echo -e "########## no R2 found for sample=$sample barcode=$barcode ###########"
			continue
		fi

		#echo -e "$sample $barcode $R1 $R2\n" >> $infiles
		echo -e "$sample $barcode $R1 $R2" >> $infiles
	done
	
	echo -e "\n############  end processing $sample.                      #########" 

done < $samplefile


echo -e "\n##############################################################################"    
echo -e "\n############   analysis loop ends  here.                             #########" 
echo -e "\n############   Done. Exiting now                                     #########"
echo -e "\n##############################################################################"  

