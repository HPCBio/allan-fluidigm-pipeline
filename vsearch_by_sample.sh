#!/bin/bash
set -x
echo `date`

if [ $# -lt 6 ]
then
     echo -e "Program $0 to launch vsearch for samples and amplicons specified in the input\nAnd to produce summaries of the results by SAMPLE\n\n";
     echo -e "Program $0 stopped at line $LINENO\nREASON=Parameters mismatch. \n\n"
     echo -e "Rerun this program like this:\n$0 <amplicon_names.txt> <sample_names.txt> <database.fasta> <database.tsv> <max_hits> <min_percent_identity> <min_coverage>\n"
     echo -e "Example:\n$0 amplicons.txt samples.txt database.fa database.tsv 3 0.90 50\n\n";
     exit 1;
fi

set +x;
echo -e "\n\n##############################################################################">&2;  
echo -e "\n\n############   defining variables for inputs and outputs     #################">&2;
echo -e "\n\n##############################################################################" >&2; 
set -x; 

ampliconfile=$1   #parameter for specifying file with amplicon names, one per line
samplefile=$2     #parameter for specifying file with sample names, one per line
database=$3       #parameter for specifying amplicon database in fasta format to be used for the search 
xrefdb=$4         #parameter for specifying amplicon database in tsv format
maxhits=$5        #parameter for specifying maxinum number of hits
perc_id=$6        #parameter for specifying minumum percent identity, a real number between 0.0 and 1.0
min_coverage=$7   #parameter for specifying minimum alignment coverage, an integer number between  0 and 100

echo -e "\n############  global variables\n" 

email=lpfreder@illinois.edu
rootdir=/home/groups/hpcbio_shared/brianallan_group
inputdir=$rootdir/results/2016-08-09-prepare-reads
tmpdir=$rootdir/src/tmp_logs
outputdir=$rootdir/results/2016-08-18_vsearch_summary_by_sample_with_taxid
scriptdir=$rootdir/src

echo -e "\n############  default values for parameters\n" 

default_perc_id=0.90
dafault_maxhits_noV=10
dafault_maxhits_V=10
default_mincov=50


echo -e "\n############  load tools\n" 

vsearchMod=vsearch/1.0.7

echo -e "\n############  cluster variables\n"   

nodes=1
threads=2
queue=default
mem=8gb

set +x;
echo -e "\n\n##############################################################################">&2;      
echo -e "\n\n############   DO NOT EDIT THE SCRIPT BEYOND THIS POINT      #################">&2;
echo -e "\n\n##############################################################################">&2; 
set -x; 

set +x;
echo -e "\n\n##############################################################################">&2;      
echo -e "\n\n############   Main Program starts here                      #################"
echo -e "\n\n##############################################################################">&2; 
set -x;  

echo -e checking that the program exists
`module load $vsearchMod`
exitcode=$?
if [ $exitcode -ne 0 ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$vsearchMod module not loaded"
    exit $exitcode;
fi

echo -e checking that the file with amplicon names was specified at the command line and that it exists
if [ `expr ${#ampliconfile}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=file with amplicon names was not specified"
    exit 1;
elif [ ! -s $ampliconfile ]
then 
    echo -e "$0 stopped at line $LINENO\sREASON=$ampliconfile file not found"
    exit 1;
fi

echo -e checking that the file with sample names was specified at the command line and that it exists
if [ `expr ${#samplefile}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=file with sample names was not specified"
    exit 1;
elif [ ! -s $samplefile ]
then 
    echo -e "$0 stopped at line $LINENO\sREASON=$samplefile file not found"
    exit 1;
fi


echo -e checking that the database was specified at the command line and that it exists
if [ `expr ${#database}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=database was not specified"
    exit 1;
elif [ ! -s $database ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$database database file not found"
    exit 1;
fi

echo -e checking that the database was specified at the command line and that it exists
if [ `expr ${#xrefdb}` -lt 1 ]
then
    echo -e "$0 stopped at line $LINENO\nREASON=database in tsv format was not specified"
    exit 1;
elif [ ! -s $xrefdb ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$xrefdb database file not found"
    exit 1;
fi
echo -e checking that the maxhits parameter was specified at the command line 
if [ `expr ${#maxhits}` -lt 1 ]
then
    echo -e "maxhits parameter was not specified. A default value will be used"
    if [[ $ampliconfile =~  /^V/ ]]
    then
         maxhits=$dafault_maxhits_V
    else    
         maxhits=$dafault_maxhits_noV
    fi
fi

echo -e checking that the maxhits parameter was specified at the command line as a real number 
if [ `expr ${#perc_id}` -gt 1 ]
then
    echo -e "percent identity parameter was not specified as a real value between 0.0 and 1.0. A default value will be used"
    perc_id=$default_perc_id
fi

echo -e checking that the min_coverage was specified at the command line 
if [ `expr ${#min_coverage}` -lt 1 ]
then
    echo -e "min_coverage parameter was not specified. A default value will be used"
    min_coverage=$default_mincov
fi

echo -e checking that input files exist
if [ ! -d $rootdir ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$rootdir directory not found"
    exit 1;
fi
if [ ! -d $inputdir ]
then
    echo -e "$0 stopped at line $LINENO\sREASON=$inputdir directory not found"
    exit 1;
fi

if [ ! -d $tmpdir ]
then
    echo -e "creating tmp folder $tmpdir"
    `mkdir $tmpdir`
fi

if [ ! -d $outputdir ]
then
    echo -e "creating output folder $outputdir";
    `mkdir $outputdir`
fi

echo -e "\n############  making a copy to the output folder\n"
echo -e "\n############  it doesn't take a lot of room and it is good for self-documentation\n"
 

cp $ampliconfile $outputdir/AMPLICON_names.txt
cp $samplefile   $outputdir/SAMPLE_names.txt
cp $xrefdb       $outputdir/database.tsv
cp $database     $outputdir/database.fa


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

echo -e "\n############  define more variables\n"

ampliconfile=$outputdir/AMPLICON_names.txt
samplefile=$outputdir/SAMPLE_names.txt
xrefdb=$outputdir/database.tsv
database=$outputdir/database.fa
prefixVsearchName=VsearchFiles
cmdline=$outputdir/vsearch_cmd
truncate -s 0 $cmdline
prefixSummaryName=vsearchSummary_percIdent_${perc_id}_maxhits_${maxhits}_mincov_$min_coverage
echo "command1 run: vsearch --usearch_global <inputfile> --db <database.fa>  --threads <threads>  --notrunclabels --id $perc_id --maxhits $maxhits  --blast6out <outputfile>" > $cmdline
echo "command2 run: filter_vsearch_output $min_coverage <outputfile.tsv> <filtered.outputfile.tsv> <outputfile.aln.txt> <iltered.outputfile.aln.txt>" >> $cmdline
echo -e "\n\nSample\tAmplicon\tEXTENDED_SEQID\tTOTAL_READS_IN_DEMULTIPLEXED_SAMPLE\tTOTAL_HITS\tMIN_percIdent\tMAX_percIdent\tAVG_percIdent\tMIN_alnLen\tMAX_alnLen\tAVG_alnLen\tMIN_Coverage\tMAX_Coverage\tAVG_Coverage" >> $cmdline


set +x;
echo -e "\n\n##############################################################################" >&2;   
echo -e "\n\n############   START loop one over SAMPLES                           #########" >&2; 
echo -e "\n\n##############################################################################" >&2; 
set -x; 


while read sample
do
	if [ `expr ${#sample}` -gt 0 ]
	then
		echo -e "\n############        Processing sample=$sample           #################\n" 

		echo -e "\n############        Preparatory steps: folders, files and variables\n"		
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
		
		# file to hold filenames of all vsearch files for this sample
		VsearchFiles=$outputdir/${prefixVsearchName}_$sample
		listFile=""
		truncate -s 0 $VsearchFiles
		
		cd $outdir
		
		# file to hold all vsearch commands for this sample
		searchLines=SearchLines_$sample		
		truncate -s 0 $searchLines
		
		# files to complete summary report
 		ln -s $ampliconfile AMPLICON_names.txt
 		ln -s $xrefdb database.tsv
 		
		set +x;
		echo -e "\n\n#####################################################################\n" >&2;
		echo -e "############    START loop two over amplicons           #################\n" >&2;
		echo -e "#########################################################################\n" >&2; 
		set -x; 

        
		while read amplicon
		do
			if [ `expr ${#amplicon}` -gt 1 ]
			then

				echo -e "\n############ processing $amplicon  ###############\n"
				
				prefix=${amplicon}_${sample}
				query_file=${prefix}_SearchReady.fasta
				inputfile=$inputdir/$amplicon/$prefix/$query_file	            
				outputfile1=$outdir/${prefix}_vsearch_temp_prefiltered.tsv
				outputfile2=$outdir/${prefix}_vsearch.tsv
				outputaln1=$outdir/${prefix}_vsearch_temp_prefiltered.aln.txt
				outputaln2=$outdir/${prefix}_vsearch.aln.txt

				if [ ! -s $inputfile ]
				then
					echo -e "\n############ $inputfile empty file. skipping this entry #######\n" 
				else
					echo -e "\n############ creating the search/filter commands for $inputfile ########\n"
					echo "vsearch --usearch_global $inputfile --db $database  --threads $threads  --notrunclabels  --id $perc_id --maxhits $maxhits  --alnout $outputaln1 --blast6out $outputfile1 " >> $searchLines
					echo "perl $scriptdir/filter_vsearch_output.pl $min_coverage $outputfile1 $outputfile2 $outputaln1 $outputaln2" >> $searchLines
					echo "rm $outputfile1" >> $searchLines
					echo "rm $outputaln1" >> $searchLines
					listfiles=${listfiles}:${outputfile2}
				fi #inputfile
			fi # empty line
		done  < $ampliconfile

		set +x;
		echo -e "\n\n#####################################################################\n"
		echo -e "############   END loop two over AMPLICONS                        #######\n"
		echo -e "############   STILL inside the loop one over SAMPLES             #######\n"
		echo -e "#########################################################################\n"


		echo -e "\n\n#####################################################################\n"
		echo -e "############   Launching the qsub job with the searches           #######\n"
		echo -e "#########################################################################\n\n" >&2; 
		set -x; 

		echo $listfiles > $VsearchFiles
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
		echo -e "echo running vsearches for sample=$sample" >>  $qsub1
		cat $searchLines >>  $qsub1
		echo -e "echo \`date\`" >>  $qsub1
		echo -e "echo running vsearch summary for sample=$sample" >>  $qsub1
		echo -e "perl $scriptdir/vsearch_calcSummary_by_sample.pl  $outdir $sample $VsearchFiles $outputfile $maxhits " >> $qsub1
		echo -e "unlink AMPLICON_names.txt" >> $qsub1
		echo -e "unlink database.tsv" >> $qsub1		
		`chmod g+r $qsub1`
		thisjob=`qsub $qsub1`
		rm $searchLines
     fi # empty sample
done < $samplefile

set +x;
echo -e "\n\n###############################################################################\n">&2;
echo -e "############    END loop one over SAMPLES                              ############\n">&2;
echo -e "###################################################################################\n">&2; 
set -x; 


set +x;
echo -e "\n\n###############################################################################\n"
echo -e "############    Exiting now                                   #################\n"
echo -e "###############################################################################\n\n" >&2; set -x; 