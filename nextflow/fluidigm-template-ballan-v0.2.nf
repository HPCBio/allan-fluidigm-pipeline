#!/usr/bin/env nextflow
/*
B_ALLAN FLUIDIGM PIPELINE
This is the generic template to process B-Allan Fluidigm data. It performs read preparation and clustering with VSEARCH
*/


/*Some hacks to deal with dates in output paths*/

import java.text.SimpleDateFormat

// version
version = 0.2


// Credit to Phil Ewels for this segment to generate the help message at the command-line

def helpMessage() {
    log.info"""
    =========================================
     B_ALLAN FLUIDIGM PIPELINE v${version} 
    =========================================

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run -c <nxf-amplicon.conf>  <fluidigm-template-amplicon.nf>
    
    where nxf-amplicon.conf is the configuration file for a particular amplicon
    
    To override existing values from the command line, please use these parameters:

	--Amplicon                    The name of the Amplicon as it appears in the sequencing report and in the filenames of the reads
	--SeqLot                      Sequencing lot or batch number that indicates the provenance of the Amplicon. Valid value: String without spaces
	--email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	--reads                       Input data( this string must be surrounded with quotes and should contain a full path and regular expression)
	--projectdir                  Path to project directory
	--scriptdir                   Path to directory with scripts
	--dbdir                       Path to directory with database of amplicons
	--resultsdir                  Path where results will be written to
	--out			      Filename of output report. Please specify its full path

    Read preparation options. Parameters for R2 will be ignored if skipStitch is set to YES:
	--skipStitch                  Skip stitching of reads. This value depends on the amplicon. Valid values: 'YES' or 'NO'
	--IsVtarget                   This value depends on the amplicon. Set it to 'YES' for all-V targets such as V1 or V4. Set it to 'NO' for all other targets.
	--readLen                     Read length. Valid value: a positive integer
	--R1minlen                    Minimum sequence length after qctrim of R1. Valid value: a positive integer
	--R1maxlen                    Maximum sequence length of R1 after qctrim. Valid value: a positive integer
	--R1adaptorLen                Length of adaptor for R1. This value comes from the sequencing report. Valid value: a positive integer
	--R2minlen                    Minimum sequence length of R2 after qctrim. Valid value: a positive integer 
	--R2maxlen                    Maximum sequence length of R2 after qctrim. Valid value: a positive integer   
	--R2adaptorLen                Length of adaptor for R2. This value comes from the sequencing report. Valid value: a positive integer
	--maxSeqLen                   Maximum sequence length of STITCHED read. Valid value: a positive integer
	--minSeqLen                   Minimum sequence length of STITCHED read. Valid value: a positive integer
	--minqual                     Minimum base quality. Valid value: a positive integer <= 40                        
	--minoverlap                  Minimum length of overlap. Valid value: a positive integer. Default=20

   Clustering options:
	--percidentity  = '0.90'      percent identity cutoff       
	--maxhitsV      = '50'        maximum number of hits to report per query; applicable to all V amplicons              
	--maxhitsNOV    = '10'        maximum number of hits to report per query; applicable to all NON-V amplicons 
	--mincov        = '50'        minimum coverage cutoff to use when post-processing hits   
	--percentile    = '5'         percentile cutoff
      
    Other optional parameters:
        -with-dag                     Graphic of the execution report in PDF format. Valid value:a file name
        -with-timeline                Timeline of all processes executed durin the execution. Report in HTML format. Valid value:a file name
	-qs                           Maximum number of concurrent processes.  Valid value: a positive integer. Default value 1. Recommended value 20.  
        --help                        To generate this message
     
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Configurable variables

// To convert this script into a template, simply set the value to false for all variables that start with params. 

/* Main variables*/
params.Amplicon      = false
params.SeqLot        = false
params.email         = false
params.projectdir    = false
params.reads         = false
params.srcdir        = false
params.dbdir         = false
params.resultsdir    = false

/*Variables for read preparation. parameters for R2 will be ignored if stitching is not performed*/
params.skipStitch    = false
params.IsVtarget     = false
params.readLen       = false
params.R1minlen      = false
params.R1maxlen      = false
params.R1adaptorLen  = false
params.R2minlen      = false
params.R2maxlen      = false
params.R2adaptorLen  = false
params.maxSeqLen     = false
params.minSeqLen     = false
params.minqual       = false
params.minoverlap    = false
params.maxhits       = '10'

/*Variables for clustering with VSEARCH*/
params.percidentity  = false
params.maxhitsV      = false
params.maxhitsNOV    = false
params.mincov        = false
params.percentile    = false

/* Path for folders with temp and final results  */
timestamp            = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())
readprepPath         = "${params.resultsdir}/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/readprep"
vsearchPath          = "${params.resultsdir}/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/vsearchResults"
params.out           = "${vsearchPath}/${params.Amplicon}_seqLot${params.SeqLot}_vsearchSummaryReport.txt"

/* Path for internal DB, programs and scripts*/
params.adaptors             = "${params.srcdir}/adaptors.fa"
params.ScriptFilterVsearch  = "${params.srcdir}/filter_vsearch_output.pl" 
params.ScriptSummaryVsearch = "${params.srcdir}/vsearch_calcSummary.pl" 
params.DBfasta              = "${params.dbdir}/DB_version2_combined_and_filtered_2017.fasta"
params.DBtsv                = "${params.dbdir}/DB_version2_combined_and_filtered_2017.tsv"



/*Biocluster2 options. List memory in gigabytes like in suggestions below*/
myQueue       = 'budget'
trimMemory    = '5'
trimCPU       = '4'           // this line becomes the slurm directive for number of CPUs
trimThreads   = '4'           // this line becomes the param number-of-threads used by the triming tool
fastqcMemory  = '2'
fastqcCPU     = '1'
fastxMemory   = '2'
fastxCPU      = '1'
vsearchMemory = '5'
vsearchCPU    = '4'           // this line becomes the slurm directive for number of CPUs
vsearchThread = '4'           // this line becomes the param number-of-threads used by the vsearch tool 
pearMemory    = '5'
pearCPU       = '4'

/*Biocluster2 module names and versions*/
pearMod       = 'PEAR/0.9.8-IGB-gcc-4.9.4'
fastxMod      = 'FASTX-Toolkit/0.0.14-IGB-gcc-4.9.4'
vsearchMod    = 'VSEARCH/2.4.3-IGB-gcc-4.9.4'
trimVersion   = '0.36'
trimMod       = "Trimmomatic/0.36-Java-1.8.0_121"
trimHome      = '\$EBROOTTRIMMOMATIC'
trimJar       = "${trimHome}/trimmomatic-${trimVersion}.jar"
trimJava      = 'Java/1.8.0_121'
nxfJava       = 'Java/1.8.0_152'
params.Java   = " -Xmx2g -XX:ParallelGCThreads=${trimThreads} "  // this line to force java to behave correctly with resources

/* trim options for stitching case */
R1trimOptions   = " HEADCROP:${params.R1adaptorLen} "
R2trimOptions   = " HEADCROP:${params.R2adaptorLen} "
stitchOptions   = " -v ${params.minoverlap} -t ${params.minSeqLen} -m ${params.maxSeqLen} -q ${params.minqual} "

/* trim options for non-stitching case */
trimOptions     = " HEADCROP:${params.R1adaptorLen} ILLUMINACLIP:${params.adaptors}:1:10:7 LEADING:${params.minqual} TRAILING:${params.minqual} SLIDINGWINDOW:4:20 CROP:${params.R1maxlen} "



// Validate inputs
if(  params.skipStitch != 'YES'  &&  params.skipStitch != 'NO' ){
    exit 1, "Invalid value for ${params.skipStitch}.  Valid options: 'YES', 'NO'\n"
}

if( params.skipStitch == 'YES' ){ 
	params.singleEnd     = '1'
	if( !params.R1maxlen || !params.R1adaptorLen || !params.minqual ){ exit 1, "Missing values for R1maxlen R1adaptorLen minqual"	}
	params.readprep      = " trimmomatic: ${trimOptions} "
}

if( params.skipStitch == 'NO' ){ 
	params.singleEnd     = '2' 
	if( !params.R1adaptorLen || !params.R2adaptorLen || !params.minoverlap ){ exit 1, "Missing values for R1adaptorLen R2adaptorLen minoverlap" }	
	if( !params.minSeqLen || !params.maxSeqLen || !params.minqual ){ exit 1, "Missing values for minSeqLen maxSeqLen minqual" }	
	params.readprep      = " trimmomatic: ${R1trimOptions} + ${R1trimOptions} PEAR: ${stitchOptions} "
}


// Create a channel for input read files. It should work with single reads and with paired reads

Channel
    .fromFilePairs( params.reads, size: -1, flat: true )
    .set { read_files_trimming  }


// Redmine Ticket must be set for project mgmt tracking
params.ticket = 7000
if (params.ticket == 0) exit 1, "Must set Redmine ticket for pipeline summary to be sent"


// Credit to Phil Ewels for this little bit of self-documentation
run_info = """
          
====================================
 B_Allan FLUIDIGM PIPELINE: Fluidigm read preparation and clustering for amplicon: ${params.Amplicon} batch: ${params.SeqLot}
====================================

Pipeline ver                : ${version}
Amplicon                    : ${params.Amplicon}
SequencingLot               : ${params.SeqLot}
Reads                       : ${params.reads}
Output file                 : ${params.out}


Read Processing options
SKIP stitching reads        : ${params.skipStitch}
read length                 : ${params.readLen}
read preparation cmds       : ${params.readprep}

Clustering options
percent_identity cutoff     : ${params.percidentity}
min_coverage cutoff         : ${params.mincov}
is V target                 : ${params.IsVtarget}
maxhits to keep for non-V targets : ${params.maxhitsNOV}
maxhits to keep for all-V targets : ${params.maxhitsV}


Tools and versions
VSEARCH                     : ${vsearchMod}
Fastx-Toolkit               : ${fastxMod}
Trimmomatic                 : ${trimMod}

Paths and other options
Current user                : $USER
Current path                : $PWD
Script dir                  : ${params.srcdir}
Database dir                : ${params.dbdir}
Output dir                  : ${params.resultsdir} 
Project dir                 : ${params.projectdir}
Ticket                      : ${params.ticket}

====================================

"""

log.info run_info


/*
* Step 1. ReadPreparation
* In this process we perform these analyses:
* --Trimming done with Trimmomatic
* --Stitching done with PEAR. It may be skipped for certains targets.
* --Format conversion with fastx-toolkit
*
* Warnings: 
* there are some tweaks that were added to make Nextflow/0.26 happy
* --The java cmd now includes params for managing heap size, -Xmx, and number of threads -XX:ParallelGCThreads
* --We perform swapping of Java modules to resolve conflicts with different versions used by nextflow and other java-based apps 
*/

process readPreparation {
    executor             'slurm'
    clusterOptions       " -n $trimCPU "
    queue                myQueue
    memory               "$trimMemory GB"
    module               trimMod, pearMod, fastxMod
    publishDir           readprepPath, mode: 'copy'

    input:
    set val(id), file(read1), file(read2) from read_files_trimming

    output:
    file '*fasta' into trimmedFastaData
    file "*{trimmed.fastq,SearchReady.fasta}"

    script:
    
    if ( params.singleEnd == '1' ) {
	    """
	    #!/bin/bash

	    fileName=\$(basename "${read1}")
	    baseName="\${fileName%_*_R1.fastq}"

	    truncate -s 0 \${baseName}_SearchReady.fasta      

	    module unload ${nxfJava}
	    module load ${trimJava}
	    
	    java ${params.Java} -jar ${trimJar} SE -threads ${trimThreads} -phred33 \
	    ${read1} \${baseName}.trimmed.fastq ${trimOptions} 2> \${baseName}.trimrun.log 

	    if [ -s \${baseName}.trimmed.fastq ]
	    then
	       fastq_to_fasta -v -i \${baseName}.trimmed.fastq -o \${baseName}_SearchReady.fasta
	    fi

	    module unload ${trimJava}
	    module load ${nxfJava}
	    
	    """	    
    } else { 
    	    """ 
	    #!/bin/bash

	    fileName=\$(basename "${read1}")
	    baseName="\${fileName%_*_R1.fastq}"

	    truncate -s 0 \${baseName}_SearchReady.fasta

	    module unload ${nxfJava}
	    module load ${trimJava}

	    java  ${params.Java} -jar ${trimJar} SE -threads ${trimThreads} -phred33 $read1 \${baseName}.qctrim.R1.fastq ${R1trimOptions} 2> \${baseName}.R1.trimrun.log 
	    java  ${params.Java} -jar ${trimJar} SE -threads ${trimThreads} -phred33 $read2 \${baseName}.qctrim.R2.fastq ${R2trimOptions} 2> \${baseName}.R2.trimrun.log    

	    pear ${stitchOptions} -f \${baseName}.qctrim.R1.fastq -r \${baseName}.qctrim.R2.fastq -o \${baseName}.stitch

	    cat \${baseName}.stitch.assembled.fastq \${baseName}.stitch.unassembled.forward.fastq >> \${baseName}.trimmed.fastq

	    if [ -s \${baseName}.trimmed.fastq ]
	    then
		fastq_to_fasta -v -i \${baseName}.trimmed.fastq -o \${baseName}_SearchReady.fasta
	    fi

	    module unload ${trimJava}
	    module load ${nxfJava}

	    """
    }

}

/*
* Step 2. RunVsearch
* This process runs vsearch between fluidigm reads and a customized database of fluidigm targets.
* The top vsearch hits are selected,then stats are calculated and finally a summary is generated.
* We need to check for empty files after each calculation because nextflow breaks down when there is no data to feed the pipe with
*/

process runVsearch {
	executor             'slurm'
	clusterOptions       " -n $vsearchCPU "
	queue                myQueue
	memory               "$vsearchMemory GB"
	module               vsearchMod
	publishDir           vsearchPath, mode: 'copy'

	input:
	file inputfile from trimmedFastaData


	output:
	file 'out' into summaryChannel    
	file "*"


	"""   
	#!/bin/bash

	fileName=\$(basename "${inputfile}")
	baseName="\${fileName%_SearchReady.fasta}"


	truncate -s 0 \${baseName}_vsearch.tsv
	truncate -s 0 \${baseName}_vsearchSummary.txt
	truncate -s 0 out


	if [ "${params.IsVtarget}"  = "YES" ]
	then
		params.maxhits=\${params.maxhitsV}
	else
		params.maxhits=\${params.maxhitsNOV}	
	fi

	if [ -s "${inputfile}" ]
	then

		echo ${inputfile} file has data. Run vsearch

		vsearch --usearch_global ${inputfile} \
		--db ${params.DBfasta} \
		--notrunclabels    --threads ${vsearchCPU} \
		--id        ${params.percidentity} \
		--maxhits   ${params.maxhits}  \
		--alnout    \${baseName}_vsearch_prefiltered_aln.txt \
		--blast6out \${baseName}_vsearch_prefiltered.tsv


		perl ${params.ScriptFilterVsearch}   ${params.mincov} \
		\${baseName}_vsearch_prefiltered.tsv \
		\${baseName}_vsearch.tsv \
		\${baseName}_vsearch_prefiltered_aln.txt \
		\${baseName}_vsearch.aln.txt

		if [ -s \${baseName}_vsearch.tsv ]
		then
			echo \${baseName}_vsearch.tsv file has data. Run summary
		
			perl ${params.ScriptSummaryVsearch} \
			\${baseName}_vsearch.tsv \
			\${baseName}_vsearchSummary.txt \
			${params.maxhits} ${params.percentile} 
		fi

		if [ -s \${baseName}_vsearchSummary.txt ]
		then
			echo \${baseName}_vsearchSummary.txt file has data. Send to output channel
			
			cat \${baseName}_vsearchSummary.txt > out
			
			header="HEADER_${params.Amplicon}_seqLot${params.SeqLot}_vsearchSummaryReport.txt"
			truncate -s 0  \${header}
			
			echo -e "PARAMETERS:" >> \${header}
			echo -e "Amplicon=${params.Amplicon}  SequencingLot=${params.SeqLot} percIdent=${params.percidentity} maxhits=${params.maxhits} mincov=${params.mincov} percentile=${params.percentile}" >> \${header}
			echo -e "\nCMDS:" >> \${header}
			echo -e "command1 run: vsearch --usearch_global <inputfile> --db <database.fa>  --threads <threads>  --notrunclabels --id ${params.percidentity} --maxhits ${params.maxhitsNOV}  --blast6out <outputfile>" >> \${header}
			echo -e "command2 run: filter_vsearch_output <mincov> <outputfile.tsv> <filtered.outputfile.tsv> <outputfile.aln.txt> <filtered.outputfile.aln.txt>" >> \${header}
			echo -e "\nCOLUMN HEADERS:" >> \${header}
			echo -e "Sample\tAmplicon\tEXTENDED_SEQID\tTOTAL_READS_IN_DEMULTIPLEXED_SAMPLE\tTOTAL_HITS\tMIN_percIdent\tMAX_percIdent\tAVG_percIdent\tMIN_alnLen\tMAX_alnLen\tAVG_alnLen\tMIN_Coverage\tMAX_Coverage\tAVG_Coverage\t${params.percentile}_PERCENTILE" >> \${header}
		fi
	fi
   
    	"""

} 


summaryChannel 
   .collectFile(name: file(params.out))
   .println { "Result saved to file: $it" }
      

workflow.onComplete {
    def subject = "[Task #${params.ticket}] B_Allan FLUIDIGM PIPELINE: Fluidigm read preparation and clustering for amplicon: ${params.Amplicon} batch: ${params.SeqLot}"

    final_log = """
Pipeline execution summary
---------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}
---------------------------
    """
    
    ['mail', '-s', subject, params.email].execute() << "${final_log}\n${run_info}"
}

