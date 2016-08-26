#!/bin/sh
# program:     prepare_fluidigmInput.sh
# author:      Gloria Rendon
# date:        April, 2016
# description: a program that takes as input files of fluidigm-demultiplexed paired reads
#              and runs: trim, stitch, fastq_to_fasta to get them ready for vsearch
#########################################################################################
infiles=$1
set -x
echo `date`

if [ $# != 1 ]
then
     echo -e "\n\nprogram $0 stopped at line $LINENO\nREASON=Parameters mismatch\n\n\nRerun like this: $0 <infiles>\n\n"
     echo -e "where <infiles> is a file with four columns separated by spaces: \n"
     echo -e "col1=amplicon col2=demultiplexedSample col3=read1 col4=read2\n\n\n"
     exit 1;
fi

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   infiles is a text file with four columns separated by space   #"
echo -e "\n############   col1=amplicon col2=demultiplexedSample col3=read1 col4=read2  #"
echo -e "\n##############################################################################" >&2; set -x;



set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining variables for inputs and outputs     #################"
echo -e "\n##############################################################################" >&2; set -x; 

rootdir=/home/groups/hpcbio_shared/brianallan_group
outdir=$rootdir/results/2016-08-09-prepare-reads  
tmpdir=$rootdir/src/tmp_logs
email=lpfreder@illinois.edu


set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining the tools                             ################"
echo -e "\n##############################################################################">&2; set -x;  

trimMod=trimmomatic/0.33
fastqcMod=fastqc/0.11.4
pearMod=pear/0.9.5
fastxMod=fastx_toolkit/0.0.14
adapters=/home/apps/trimmomatic/trimmomatic-0.33/adapters/TruSeq3-PE-2.fa

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   defining variables for running jobs in the cluster     ########"
echo -e "\n##############################################################################">&2; set -x;  

nodes=1
threads=8
queue=default
mem=4gb

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   DO NOT EDIT THE CODE BEYOND THIS POINT         ################"
echo -e "\n##############################################################################">&2; set -x;  


set +x;
echo -e "\n##############################################################################"      
echo -e "\n############   sanity  check                                 #################"
echo -e "\n##############################################################################">&2; set -x;  

if [ `expr ${#infiles}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=infiles not specified. Exiting now."
    exit 1;
fi
if [ ! -e $infiles ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=input file not found $infiles"
    exit 1;
fi
if [ ! -d $rootdir ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$rootdir directory not found"
    exit 1;
fi
if [ ! -d $outdir ]
then
    echo -e "creating output folder $outdir";
    `mkdir -p $outdir`
fi
if [ ! -d $tmpdir ]
then
    echo -e "creating output folder $tmpdir";
    `mkdir $tmpdir`
fi

set +x;
echo -e "\n##############################################################################"  
echo -e "\n############   parameter checked. Everything ok              #################"
echo -e "\n##############################################################################"  

echo -e "\n##############################################################################"    
echo -e "\n############   analysis loop start here.                             #########" 
echo -e "\n############   Reading infiles, one line at a time                   #########"
echo -e "\n############   col1=amplicon col2=demultiplexedSample col3-read1 col4=read2  #"
echo -e "\n##############################################################################" >&2; set -x; 

# these files will contain the samples that will not be analyzed further after this point
skipped=$outdir/input_files_with_no_reads
emptyOut=$outdir/output_files_with_no_reads
truncate -s 0 $skipped
truncate -s 0 $emptyOut

while read sampledetail
do
    echo -e "\n############# processing next line in file...                #################\n"

    if [ `expr ${#sampledetail}` -lt 1 ]
    then
        set +x;
        echo -e "\n##############################################################################"  
        echo -e "\n############ skipping empty line                               ###############"
        echo -e "\n##############################################################################" >&2; set -x;  
    else
        set +x;
        echo -e "\n##############################################################################"  
        echo -e "\n############ parsing line: $sampledetail                       ###############"
        echo -e "\n##############################################################################"  
        echo -e "\n###    col1=amplicon col2=demultiplexedSample col3-read1 col4=read2        ###\n" >&2; set -x;        


        amplicon=$( echo $sampledetail | cut -d ' ' -f1 )
        demultiplexed=$( echo $sampledetail | cut -d ' ' -f2 )
        R1=$( echo $sampledetail | cut -d ' ' -f3 )
        R2=$( echo $sampledetail | cut -d ' ' -f4 )
        empty="NO"

        set +x;
        echo -e "\n##############################################################################" 
        echo -e "\n############ lined parsed. now validating info                 ###############" 
        echo -e "\n##############################################################################" >&2; set -x;
        
	if [ `expr ${#amplicon}` -lt 1 ]
 	then
	    echo -e "$0 stopped at line $LINENO\nREASON=empty string for ampliclon"
	    exit 1;
        fi
           
	if [ `expr ${#demultiplexed}` -lt 1 ]
 	then
	    echo -e "$0 stopped at line $LINENO\nREASON=empty string for sample name"
	    exit 1;
        fi
           
	if [ `expr ${#R1}` -lt 1 ]
 	then
	    echo -e "$0 stopped at line $LINENO\nREASON=empty string for R1"
	    exit 1;
        fi    
 	if [ `expr ${#R2}` -lt 1 ]
 	then
	    echo -e "$0 stopped at line $LINENO\nREASON=empty string for R2"
	    exit 1;
        fi 
        
	if [ ! -s $R1 ]
	then
	    echo -e "empty file read1 $R1"
	    empty="Y";
        fi

	if [ ! -s $R2 ]
	then
	    echo -e "empty file read2 $R2"
	    empty="Y";
        fi
           
        if [ $empty != "NO" ]
        then
            set +x;
            echo -e "\n##############################################################################" 
            echo -e "\n############ no reads. skipping this line   $sampledetails     ###############"
            echo -e "\n##############################################################################" >&2; set -x;             
            echo -e "$R1\t$R2\n" >> $skipped
        else
            set +x;
            echo -e "\n##############################################################################" 
            echo -e "\n############ info validated. proceed with analysis             ###############"
            echo -e "\n##############################################################################"    
            echo -e "\n############ define output folders and files  for the nalaysis ###############"  
            echo -e "\n##############################################################################" >&2; set -x; 
            
            outputdir=$outdir/$amplicon/$demultiplexed
            fqdir2=$outputdir/FastQC-trimmed
            fqdir1=$outputdir/FastQC-raw               
            mkdir -p $fqdir1
            mkdir -p $fqdir2              

            b1=${demultiplexed}_R1
            b2=${demultiplexed}_R2
            pairedR1=${b1}.paired.fq
            pairedR2=${b2}.paired.fq
            merged=${demultiplexed}_stitched
            fq2fa=${demultiplexed}_SearchReady.fasta

            set +x;
            echo -e "\n##############################################################################" 
            echo -e "\n############ output folder: $outputdir                          ##############"
            echo -e "\n############ output file: $fq2fa                                ##############" 
            echo -e "\n############ creating the qsub script and launching it          ##############"
            echo -e "\n##############################################################################">&2; set -x; 
            
	    qsub1=$tmpdir/qsub.prepareReads.$demultiplexed
	       
	    echo "#PBS -S /bin/bash" > $qsub1
	    echo "#PBS -N prepareReads_${demultiplexed}" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -e $tmpdir/qsub.prepareReads.${demultiplexed}.er" >> $qsub1
	    echo "#PBS -o $tmpdir/qsub.prepareReads.${demultiplexed}.ou" >> $qsub1
	    echo "#PBS -l nodes=$nodes:ppn=$threads" >> $qsub1
	    echo "#PBS -q $queue" >> $qsub1
            echo "set -x" >> $qsub1
            echo -e "\n\n" >> $qsub1
            echo "module load $fastqcMod"  >> $qsub1
            echo "cd $outputdir"  >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step1 fastqc on raw reads" >> $qsub1
            echo "fastqc -o $fqdir1 -t $threads $R1" >> $qsub1
            echo "fastqc -o $fqdir1 -t $threads $R2" >> $qsub1
            echo -e "\n\n" >> $qsub1
            echo "module unload $fastqcMod"  >> $qsub1
            echo "module unload java/1.8.0_65"  >> $qsub1
            echo "module load $trimMod"  >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step2 trim raw reads" >> $qsub1
	    echo "java -classpath /home/apps/trimmomatic/trimmomatic-0.33/trimmomatic-0.33.jar \
   org.usadellab.trimmomatic.TrimmomaticPE \
   -threads $threads \
   -phred33 \
   -trimlog $outputdir/${demultiplexed}.trim.log \
   $R1 $R2 \
   $outputdir/${b1}.paired.fq $outputdir/${b1}.unpaired.fq \
   $outputdir/${b2}.paired.fq $outputdir/${b2}.unpaired.fq \
   HEADCROP:20 SLIDINGWINDOW:3:15 " >> $qsub1
            echo "\$exitcode=\$?"  >> $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "if [ \$exitcode -ne 0 ]"  >> $qsub1
            echo "then"  >> $qsub1
            echo "   echo program stopped at line=$LINENO reason=trimmomatic failed on $sampledetail"  >> $qsub1
            echo "   exit \$exitcode"  >> $qsub1
            echo "fi"  >> $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "module unload $trimMod"  >> $qsub1
            echo "module unload java"  >> $qsub1
            echo "module load $fastqcMod"  >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step3 fastqc on trimmed reads" >> $qsub1
            echo "fastqc -o $fqdir2 -t $threads $pairedR1" >> $qsub1
            echo "fastqc -o $fqdir2 -t $threads $pairedR2" >> $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "module unload $fastqcMod"  >> $qsub1
            echo "module load $pearMod"  >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step4 stitch trimmed reads" >> $qsub1
            echo "pear -f $pairedR1 -r $pairedR2 -o $merged"   >> $qsub1
            echo "\$exitcode=\$?"  >> $qsub1            
            echo "if [ \$exitcode -ne 0 ]"  >> $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "then"  >> $qsub1
            echo "   echo program stopped at line=$LINENO reason=pear failed on $sampledetail"  >> $qsub1
            echo "   exit \$exitcode"  >> $qsub1
            echo "fi"  >> $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "module unload $pearMod"  >> $qsub1
            echo "module unload perl/5.16.1_hpcbio" >> $qsub1
            echo "module load $fastxMod"  >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step5  convert fq2fa and concat files" >> $qsub1
            echo "fastq_to_fasta -v -i ${merged}.assembled.fastq -o ${merged}.assembled.fasta"   >> $qsub1
            echo "fastq_to_fasta -v -i ${merged}.unassembled.forward.fastq -o ${merged}.unassembled.forward.fasta"   >> $qsub1            
            echo "cat ${merged}.assembled.fasta ${merged}.unassembled.forward.fasta > $fq2fa" >> $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo step6 report empty fasta files" >> $qsub1
            echo "if [ ! -s $fq2fa ]" >>  $qsub1
            echo "then" >>  $qsub1
            echo "   echo $fq2fa is empty" >> $qsub1
            echo "   cat $fq2fa >> $emptyOut" >>  $qsub1
            echo "fi" >>  $qsub1
            echo "echo `date`" >>  $qsub1
            echo -e "\n\n" >> $qsub1            
            echo "echo exiting now" >> $qsub1
	    `chmod g+r $qsub1 `
	    alnjobid=`qsub $qsub1`
	    echo `date`
	fi  # file with read
    fi     # nonempty line      
done < $infiles


echo -e "\n##############################################################################"    
echo -e "\n############   analysis loop ends  here.                             #########" 
echo -e "\n############   Done. Exiting now                                     #########"
echo -e "\n##############################################################################"  

