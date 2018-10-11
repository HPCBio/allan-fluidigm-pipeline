#!/bin/bash
# program:     prepare_fluidigmInput_forV2.sh
# author:      Gloria Rendon
# date:        August, 2016
# description: to launch vsearch for samples and amplicons specified in the configuration file and to produce summaries of the results by SAMPLE
#########################################################################################

set -x
echo `date`

set +x;
echo -e "\n\n##############################################################################" >&2;
echo -e "\n############   STARTS VSEARCH_BY_SAMPLE_FORV2.SH                       #########" >&2;
echo -e "\n##############################################################################\n\n" >&2; 
set -x;

if [ $# -lt 1 ]
then
     echo -e "Program $0 to launch vsearch for samples in SAMPLEFILE and amplicons in AMPLICONFILE \n";
     echo -e "SAMPLEFILE and AMPLICONFILE are specified in the configuration file\n\n"
     echo -e "Program $0 stopped at line $LINENO\nREASON=Parameters mismatch. \n\n"
     echo -e "Rerun this program like this:\n$0 <configuration file>\n"
     echo -e "Where <configuration file> is a text file with information needed to run this program\n\n";
     exit 1;
fi

runfile=$1

if [ ! -s $runfile ]
then
     echo "$runfile configuration file not found. Exiting now\n\n"
     exit 1;
fi


set +x;
echo -e "\n\n##############################################################################">&2;  
echo -e "\n\n############   PARSING CONFIGURATION FILE AND SANITY CHECK     ###############">&2;
echo -e "\n\n##############################################################################" >&2; 
set -x; 

ampliconfile=$( cat $runfile | grep -w AMPLICONFILE | cut -d '=' -f2 )
samplefile=$( cat $runfile | grep -w SAMPLEFILE | cut -d '=' -f2 )
database=$( cat $runfile | grep -w DATABASE_FASTA | cut -d '=' -f2 )
xrefdb=$( cat $runfile | grep -w DATABASE_TSV | cut -d '=' -f2 )
maxhits=$( cat $runfile | grep -w MAXHITS | cut -d '=' -f2 )
perc_id=$( cat $runfile | grep -w MIN_PERCENT_IDENTITY | cut -d '=' -f2 )
percentile=$( cat $runfile | grep -w PERCENTILE_VALUE | cut -d '=' -f2 )
min_coverage=$( cat $runfile | grep -w MIN_COVERAGE | cut -d '=' -f2 )
email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w PROJECTDIR | cut -d '=' -f2 )
inputdir=$( cat $runfile | grep -w PREPARED_READS_DIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TEMPDIR | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w VSEARCH_OUTPUTDIR | cut -d '=' -f2 )
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
vsearchMod=$( cat $runfile | grep -w VSEARCHMOD | cut -d '=' -f2 ) 
nodes=$( cat $runfile | grep -w NODES | cut -d '=' -f2 )
threads=$( cat $runfile | grep -w THREADS | cut -d '=' -f2 )
queue=$( cat $runfile | grep -w QUEUE | cut -d '=' -f2 )
mem=$( cat $runfile | grep -w MEMORY | cut -d '=' -f2 )
cleanup=$( cat $runfile | grep -w REMOVE_TMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

set +x; echo -e "\n\n\n############ default values for certain search parameters\n" >&2; set -x;

default_perc_id=0.90
dafault_maxhits_noV=50
dafault_maxhits_V=50
default_mincov=50
default_percentile=5


set +x; echo -e "\n\n\n############ checking vsearch module\n" >&2; set -x;

`module load $vsearchMod`
exitcode=$?
if [ $exitcode -ne 0 ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$vsearchMod module not loaded"
    exit $exitcode;
fi

set +x; echo -e "\n\n\n############ checking that the file with amplicon names  exists\n" >&2; set -x;

if [ `expr ${#ampliconfile}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=file with amplicon names was not specified"
    exit 1;
elif [ ! -s $ampliconfile ]
then 
    echo -e "$0 stopped at line $LINENO\sREASON=$ampliconfile file is empty "
    exit 1;
fi


set +x; echo -e "\n\n\n############  checking that the file with sample names  exists\n" >&2; set -x;

if [ `expr ${#samplefile}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=file with sample names was not specified"
    exit 1;
elif [ ! -s $samplefile ]
then 
    echo -e "$0 stopped at line $LINENO\sREASON=$samplefile file is empty"
    exit 1;
fi


set +x; echo -e "\n\n\n############  checking that the database  exists \n" >&2; set -x;

if [ `expr ${#database}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=database in fasta format was not specified"
    exit 1;
elif [ ! -s $database ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$database database in fasta format was not found"
    exit 1;
fi

if [ `expr ${#xrefdb}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=database in tsv format was not specified"
    exit 1;
elif [ ! -s $xrefdb ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$xrefdb database in tsv format was not found"
    exit 1;
fi

set +x; echo -e "\n\n\n############   checking that the maxhits parameter was specified \n" >&2; set -x;

if [ $maxhits -lt 1 ]
then
    echo -e "maxhits parameter must be a positive integer. A default value will be used"
    if [[ $ampliconfile =~  /^V/ ]]
    then
         maxhits=$dafault_maxhits_V
    else    
         maxhits=$dafault_maxhits_noV
    fi
fi

set +x; echo -e "\n\n\n############   checking that the maxhits parameter was specified as a real number \n" >&2; set -x;

if [ $perc_id -gt 1.0 ]
then
    echo -e "percent identity parameter must be real value between 0.0 and 1.0. A default value will be used"
    perc_id=$default_perc_id
fi

set +x; echo -e "\n\n\n############  checking that the min_coverage was specified  \n" >&2; set -x;

if [ $min_coverage -lt 1 ]
then
    echo -e "min_coverage parameter must be a positive integer. A default value will be used"
    min_coverage=$default_mincov
fi

set +x; echo -e "\n\n\n############  checking that the percentile was specified  \n" >&2; set -x;

if [ $percentile -lt 1 ]
then
    echo -e "percentile parameter must be a positive integer. A default value will be used"
    percentile=$default_percentile
fi

set +x; echo -e "\n\n\n############  checking that project folder exists  \n" >&2; set -x;

if [ ! -d $rootdir ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$rootdir directory not found"
    exit 1;
fi

set +x; echo -e "\n\n\n############  checking that input folder exists  \n" >&2; set -x;

if [ ! -d $inputdir ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$inputdir directory not found"
    exit 1;
fi

set +x; echo -e "\n\n\n############  checking that tmp folder exists  \n" >&2; set -x;
if [ ! -d $tmpdir ]
then
    echo -e "creating tmp folder $tmpdir"
    `mkdir $tmpdir`
fi

set +x; echo -e "\n\n\n############  checking that output folder exists  \n" >&2; set -x;

if [ ! -d $outputdir ]
then
    echo -e "creating output folder $outputdir";
    `mkdir -p $outputdir`
else
    echo -e "resetting output folder $outputdir";
    `rm -rf $outputdir`
    `mkdir -p $outputdir`
fi

set +x;
echo -e "\n\n##############################################################################">&2;  
echo -e "\n\n############   PREP WORK: COPY FILES TO OUTPUT FOLDER          ###############">&2;
echo -e "\n\n##############################################################################" >&2; 
set -x; 

cp $ampliconfile $outputdir/AMPLICON_names.txt
cp $samplefile   $outputdir/SAMPLE_names.txt
cp $xrefdb       $outputdir/database.tsv
cp $database     $outputdir/database.fa

set +x; echo -e "\n\n\n############  checking that the copy commands worked  \n" >&2; set -x;

if [ ! -s $outputdir/AMPLICON_names.txt ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$outputdir/AMPLICON_names.txt file not copied"
    exit 1;
fi
if [ ! -s $outputdir/SAMPLE_names.txt ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$outputdir/SAMPLE_names.txt file not copied"
    exit 1;
fi
if [ ! -s $outputdir/database.tsv ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$outputdir/database.tsv file not copied"
    exit 1;
fi
if [ ! -s $outputdir/database.fa ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$outputdir/database.tsv file not copied"
    exit 1;
fi

set +x;
echo -e "\n\n##############################################################################">&2;  
echo -e "\n\n############   PREP WORK: DEFINE VARIABLES AND GENERATE HEADER OF OUTPUT FILE ">&2;
echo -e "\n\n##############################################################################" >&2; 
set -x; 

numsamples=`wc -l $outputdir/SAMPLE_names.txt`           # to count the number of unique samples processed
numamplicons=`wc -l $outputdir/AMPLICON_names.txt`       # to count the number of unique amplicons processed

ampliconfile=$outputdir/AMPLICON_names.txt               # local copy of config file
samplefile=$outputdir/SAMPLE_names.txt                   # local copy of config file
xrefdb=$outputdir/database.tsv                           # local copy of config file
database=$outputdir/database.fa                          # local copy of config file
prefixVsearchName=VsearchFiles
cmdline=$outputdir/vsearch_cmd
truncate -s 0 $cmdline

set +x; echo -e "\n\n\n############  generate header for output file  \n" >&2; set -x;


prefixSummaryName=vsearchSummary_percIdent_${perc_id}_maxhits_${maxhits}_mincov_$min_coverage
echo "command1 run: vsearch --usearch_global <inputfile> --db <database.fa>  --threads <threads>  --notrunclabels --id $perc_id --maxhits $maxhits  --blast6out <outputfile>" > $cmdline
echo "command2 run: filter_vsearch_output $min_coverage <outputfile.tsv> <filtered.outputfile.tsv> <outputfile.aln.txt> <iltered.outputfile.aln.txt>" >> $cmdline
echo -e "\n\nSample\tAmplicon\tEXTENDED_SEQID\tTOTAL_READS_IN_DEMULTIPLEXED_SAMPLE\tTOTAL_HITS\tMIN_percIdent\tMAX_percIdent\tAVG_percIdent\tMIN_alnLen\tMAX_alnLen\tAVG_alnLen\tMIN_Coverage\tMAX_Coverage\tAVG_Coverage\t${percentile}_PERCENTILE" >> $cmdline

set +x;
echo -e "\n\n##############################################################################">&2;  
echo -e "\n\n############   ANALYSIS BLOCK STARTS HERE                      ###############">&2;
echo -e "\n\n##############################################################################" >&2; 
set -x; 


set +x;
echo -e "\n\n##############################################################################" >&2;   
echo -e "\n\n############   START LOOP1 OVER SAMPLES                              #########" >&2; 
echo -e "\n\n##############################################################################" >&2; 
set -x; 


while read sample
do
	if [ `expr ${#sample}` -gt 0 ]
	then
		set +x;  echo -e "\n############  Processing sample=$sample           #################\n\n\n" >&2; set -x;

		set +x;  echo -e "\n############  Preparatory steps: folders, files and variables\n\n\n">&2; set -x;


		
		outdir=$outputdir/$sample

		if [ ! -d $outdir ]
		then
			echo -e "creating output folder $outdir"
			`mkdir -p $outdir`
		else
			echo -e "output folder $outdir exists. RESETTING IT"
			rm -rf $outdir/*.tsv
			rm -rf $outdir/*.txt		
		fi
		
		cd $outdir
		
		set +x;  echo -e "\n############  file to hold all vsearch commands for this sample\n\n">&2; set -x;
		
		searchLines=SearchLines_$sample		
		truncate -s 0 $searchLines
		
		set +x;  echo -e "\n############  files to complete summary report\n\n">&2; set -x;
		
 		ln -s $ampliconfile AMPLICON_names.txt
 		ln -s $xrefdb database.tsv
 		
		set +x;
		echo -e "\n\n#####################################################################\n" >&2;
		echo -e "############    START LOOP2 OVER AMPLICONS                          #####\n" >&2;
		echo -e "############    TO PACKAGE ALL VSEARCH COMMANDS FOR THIS SAMPLE     #####\n" >&2;		
		echo -e "#########################################################################\n" >&2; 
		set -x; 

        
		while read amplicon
		do
			if [ `expr ${#amplicon}` -gt 1 ]
			then

				set +x; echo -e "\n############ processing sample=${sample} and amplicon=$amplicon  ###############\n\n" >&2; set -x;
				
				prefix=${amplicon}-${sample}
				query_file=${prefix}_SearchReady.fasta
				inputfile=$inputdir/$amplicon/$prefix/$query_file	            
				outputfile1=$outdir/${prefix}_vsearch_temp_prefiltered.tsv
				outputfile2=$outdir/${prefix}_vsearch.tsv
				outputaln1=$outdir/${prefix}_vsearch_temp_prefiltered.aln.txt
				outputaln2=$outdir/${prefix}_vsearch.aln.txt

				if [ ! -s $inputfile ]
				then
					set +x;  echo -e "\n############ $inputfile empty file. skipping this entry #######\n" >&2; set -x;
				else
					set +x;  echo -e "\n############ creating the search/filter commands for $inputfile ########\n" >&2; set -x;
					
					echo "vsearch --usearch_global $inputfile --db $database  --threads $threads  --notrunclabels  --id $perc_id --maxhits $maxhits  --alnout $outputaln1 --blast6out $outputfile1 " >> $searchLines
					echo "perl $scriptdir/filter_vsearch_output.pl $min_coverage $outputfile1 $outputfile2 $outputaln1 $outputaln2" >> $searchLines
					if [ $cleanup == "YES" ]
					then

						echo "rm $outputfile1" >> $searchLines
						echo "rm $outputaln1" >> $searchLines
					fi #end cleanup
					
				fi #end if inputfile
			fi # end if empty line
			
			set +x; echo -e "\n############ done processing sample=${sample} and amplicon=$amplicon Next one please ###############\n\n" >&2; set -x;
			
		done  < $ampliconfile

		set +x;
		echo -e "\n\n#####################################################################\n"
		echo -e "############   END LOOP2 OVER    AMPLICONS                        #######\n"
		echo -e "############   STILL INSIDE LOOP1 OVER SAMPLES                    #######\n"
		echo -e "#########################################################################\n"


		echo -e "\n\n#####################################################################\n"
		echo -e "############   Launching the qsub job with the searches           #######\n"
		echo -e "#########################################################################\n\n" >&2; 
		set -x; 

		outputfile=${prefixSummaryName}_$sample
		truncate -s 0 $outputfile
		cat $cmdline >>  $outputfile

		
		qsub1=$tmpdir/qsub.vsearch.$sample
		echo -e "#PBS -S /bin/bash" > $qsub1
		echo -e "#PBS -N vsearch.$sample" >> $qsub1
		echo -e "#PBS -M $email" >> $qsub1
		echo -e "#PBS -m ae" >> $qsub1
		echo -e "#PBS -e $tmpdir/qsub.vsearch.$sample.er" >> $qsub1
		echo -e "#PBS -o $tmpdir/qsub.vsearch.$sample.ou" >> $qsub1
		echo -e "#PBS -l nodes=$nodes:ppn=$threads" >> $qsub1
		echo -e "#PBS -q $queue" >> $qsub1
		echo -e "set -x" >> $qsub1
		echo -e "module load $vsearchMod"  >> $qsub1
		echo -e "cd $outdir"  >> $qsub1
		echo -e "echo \`date\`" >>  $qsub1
		echo -e "echo ############ running vsearches for sample=$sample" >>  $qsub1
		cat $searchLines >>  $qsub1
		echo -e "echo \`date\`" >>  $qsub1
		echo -e "echo ############ running vsearch summary for sample=$sample" >>  $qsub1
		echo -e "perl $scriptdir/vsearch_calcSummary_by_sample_forV2.pl  $outdir $sample $outputfile $maxhits $percentile" >> $qsub1
		
		if [ $cleanup == "YES" ]
		then
		        echo -e "if [ -s $outputfile ]"  >> $qsub1
		        echo -e "then"  >> $qsub1		        
			echo -e "    unlink AMPLICON_names.txt" >> $qsub1
			echo -e "    unlink database.tsv" >> $qsub1	
			echo -e "    rm $searchLines"  >> $qsub1
		        echo -e "fi"  >> $qsub1			
		fi

		`chmod g+r $qsub1`
		thisjob=`qsub $qsub1`
		
     fi # end if empty sample
     
     set +x; echo -e "\n############ end processing sample=$sample Next one please ###############\n\n" >&2; set -x;
     
done < $samplefile

set +x;
echo -e "\n\n###############################################################################\n">&2;
echo -e "############    END LOOP1 OVER    SAMPLES                              ############\n">&2;
echo -e "###################################################################################\n">&2; 
set -x; 

set +x;
echo -e "\n##############################################################################" >&2;    
echo -e "\n############   analysis ends  here.                                  #########" >&2;
echo -e "\n############   $numsamples total samples processed                   #########" >&2;
echo -e "\n############   $numamplicons total amplicons processed               #########" >&2;
echo -e "\n############   $numsamples total qsub jobs submitted                 #########" >&2;
echo -e "\n############   $numamplicons total files processed per qsub jobs     #########" >&2;
echo -e "\n##############################################################################" >&2;
echo -e "\n##############################################################################" >&2;

echo -e "\n\n##############################################################################" >&2;
echo -e "\n############   ENDS VSEARCH_BY_SAMPLE_FORV2.SH                          #########" >&2;
echo -e "\n##############################################################################\n\n" >&2; 
set -x;
