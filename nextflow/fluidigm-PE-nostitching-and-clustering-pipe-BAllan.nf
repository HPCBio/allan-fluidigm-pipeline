#!/usr/bin/env nextflow
/*
* USAGE: nextflow run fluidigm-PE-nostitching-and-clustering-pipe-Ballan.nf -qs 6
* Where  "-qs" is similar to "#PBS -t" and will only run a specified # of jobs at a time.
*
* This pipeline performs fluidigm sequence QC, processing and clustering when no stitching is involved
* Fluidigm raw reads should be demultiplexed by amplicon-sample and arranged into separate folders as per fluidigm-sequencing-pipeline-V2
* Only R1 is used for read-processing.
* This script creates hard links to data that exists in nextflow's work directory for everything
*/
/*
 
/*Some hacks*/

import java.text.SimpleDateFormat

// version
version = 0.1

//
// EDIT BLOCK STARTS HERE 
// The user can edit this section for input folder and tool parameters
//

/* Email for sending execution notifications-updates*/
params.email = "grendon@illinois.edu"

/*Path to project folder*/
projectPath = "/home/groups/hpcbio_shared/brianallan_group"

/*Name of inputs*/
params.Amplicon      = "16S8FE3_B_GA1B2"         // The name of the Amplicon as it appears in the sequencing report and in the filenames of the reads
params.SeqLot        = "20170228"                // Default: the date of the sequencing report. This string that indicate the provenance of the Amplicon. 


/*input channels*/
//rawDataPath = "$projectPath/data/rawReads/20170228_rawReads/16S8FE3_B_GA1B2"

rawDataPath = "$projectPath/data/rawReads/fortesting"
readsCh = Channel.fromPath("$rawDataPath/*_R1.fastq")
              .ifEmpty { error "Cannot find any FASTQ file in: ${rawDataPath}" }



/*folder paths for results*/

timestamp = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())
trimPath       = "$projectPath/results/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/trimmomatic"
conversionPath = "$projectPath/results/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/fastaReady"
vsearchPath    = "$projectPath/results/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/vsearchUnfiltered"
summaryPath    = "$projectPath/results/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/vsearchFinalResults"

/*Name of outputs*/
params.out           = "${summaryPath}/${params.Amplicon}_vsearchSummaryReport.txt" // The name of the output report


/*Variables for Trimming with Trimmomatic*/
params.readLen       = '250'                       // read length
params.R1minlen      = '170'                       // minimum sequence length after qctrim
params.R1adaptorLen  = '20'                        // length of adaptor for R1. This value comes from the sequencing report
params.R1maxlen      = '230'                       // maximum sequence length of R1 after qctrim (deafult value= readlen -adaptorlen)   
params.minqual       = '30'                        // minimum base quality. The range goes from 0-40                        






/*Variables for clustering with VSEARCH*/
params.percidentity  = '0.90'                      //percent identity cutoff       
params.maxhitsV      = '10'                        //maximum number of hits to report per query; applicable to all V amplicons              
params.maxhitsNOV    = '10'                        //maximum number of hits to report per query; applicable to all NON-V amplicons 
params.mincov        = '50'                        //minimum coverage cutoff to use when post-processing hits   
params.percentile    = '5'                         //percentile value

/* Path for internal DB, programs and scripts*/
params.ScriptFilterVsearch  = "$projectPath/src/nextflow/filter_vsearch_output.pl" 
params.ScriptSummaryVsearch = "$projectPath/src/nextflow/vsearch_calcSummary.pl" 
params.DBfasta              = "/home/groups/hpcbio_shared/brianallan_group/data/database/20160721_database/AmpliconsDB_with_taxid_cleaned.fa"
params.DBtsv                = "/home/groups/hpcbio_shared/brianallan_group/data/database/20160721_database/AmpliconsDB_with_taxid_cleaned.tsv"
params.adaptors             = "$projectPath/src/nextflow/adaptors.fa"

//
// EDIT BLOCK ENDS HERE 
// The user should no edit the remaining portion of the pipeline
// Proceed at your own risk
//


/*Biocluster options. List memory in gigabytes like in suggestions below*/
myQueue       = 'normal'
trimMemory    = '5'
trimCPU       = '4'
fastqcMemory  = '2'
fastqcCPU     = '1'
fastxMemory   = '2'
fastxCPU      = '1'
vsearchMemory = '2'
vsearchCPU    = '2'


/*Module versions*/
fastxMod      = 'FASTX-Toolkit/0.0.14-IGB-gcc-4.9.4'
vsearchMod    = 'VSEARCH/2.4.3-IGB-gcc-4.9.4'
trimVersion   = '0.36'
trimMod       = "Trimmomatic/0.36-Java-1.8.0_121"
trimHome      = '\$EBROOTTRIMMOMATIC'
trimJar       = "${trimHome}/trimmomatic-${trimVersion}.jar"


trimOptions   = " HEADCROP:${params.R1adaptorLen} ILLUMINACLIP:${params.adaptors}:1:10:7 CROP:${params.R1maxlen}  LEADING:${params.minqual} TRAILING:${params.minqual} SLIDINGWINDOW:4:20"




// Redmine Ticket must be set for project mgmt tracking
params.ticket = 7000
if (params.ticket == 0) exit 1, "Must set Redmine ticket for pipeline summary to be sent"


// Credit to Phil Ewels for this little bit of self-documentation
run_info = """
          
====================================
 B_Allan Fluidigm_NOStitch : Fluidigm qctrim and vsearch clustering when no stitching of reads is involved 
====================================

Pipeline ver : ${version}
Amplicon     : ${params.Amplicon}
SequencingLot: ${params.SeqLot}
Reads        : ${rawDataPath}
Database     : ${params.DBfasta}
VSEARCH      : ${vsearchMod}
Fastx-Toolkit: ${fastxMod}
Trimmomatic  : ${trimMod}
Trim options : ${trimOptions}
Vsearch options: 
percent_identity : ${params.percidentity}
min_coverage     : ${params.mincov}
maxhits_perQuery : ${params.maxhitsNOV} 
Current home : $HOME
Current user : $USER
Current path : $PWD
Script dir   : $baseDir
Working dir  : $workDir
Output dir   : $summaryPath
Project dir  : ${projectPath}
Ticket       : ${params.ticket}


====================================

"""

log.info run_info


/*
* Step 1. Trimmomatic on R1 only
*/

process qctrimOnly {
    executor   'slurm'
    cpus       trimCPU
    queue      myQueue
    memory     "$trimMemory GB"
    module     trimMod
    publishDir trimPath, mode: 'link'
    validExitStatus 0,1

    input:
    file read from readsCh 

    output:
    file '*.trimmed.fastq' into trimmedReads
    file '*.log' into trimmedFastqLogs

    script:
    """
    fileName=\$(basename "${read}")
    baseName="\${fileName%_*_R1.fastq}"
       

    java -jar ${trimJar} SE -threads $trimCPU -phred33 \
    ${read} \${baseName}.trimmed.fastq ${trimOptions} 2> \${baseName}.trimrun.log 
    """
}



/*
* Step 2. FASTQ 2 FASTA conversion on
*/

process runFQ2FAconversion {
    executor   'slurm'
    cpus       fastxCPU
    queue      myQueue
    memory     "$fastxMemory GB"
    module     fastxMod
    publishDir conversionPath, mode: 'link'

    input:
    file trimmedFQ from trimmedReads

    output:
    file '*fasta' into trimmedFastaData

    """
    fileName=\$(basename "${trimmedFQ}")
    baseName="\${fileName%.trimmed.fastq}"
    
    
    fastq_to_fasta -v -i $trimmedFQ -o \${baseName}_SearchReady.fasta
    """

}

/*
* Step 3. Run VSEARCH
*/

process runVsearch {
    executor   'slurm'
    cpus       vsearchCPU
    queue      myQueue
    memory     "$vsearchMemory GB"
    module     vsearchMod
    publishDir vsearchPath, mode: 'link'

    input:
    file inputfile from trimmedFastaData


    output:
    file "*_vsearch.tsv" into vsearchChannel
    file "*" 
    
    
    """   
    fileName=\$(basename "${inputfile}")
    baseName="\${fileName%_SearchReady.fasta}"
  
    vsearch --usearch_global ${inputfile} \
    --db ${params.DBfasta} \
    --notrunclabels    --threads ${vsearchCPU} \
    --id        ${params.percidentity} \
    --maxhits   ${params.maxhitsNOV}  \
    --alnout    \${baseName}_vsearch_prefiltered_aln.txt \
    --blast6out \${baseName}_vsearch_prefiltered.tsv
    
 
    perl ${params.ScriptFilterVsearch}   ${params.mincov} \
    \${baseName}_vsearch_prefiltered.tsv \
    \${baseName}_vsearch.tsv \
    \${baseName}_vsearch_prefiltered_aln.txt \
    \${baseName}_vsearch.aln.txt
    

    """

}

/*
* Step 4. Run VSEARCH-Summary
*/

process runSummary {
    executor   'slurm'
    cpus       1
    queue      myQueue
    memory     "4 GB"
    publishDir summaryPath, mode: 'link'
    
    input:
    file inputfile from vsearchChannel

    output:
    file 'out' into summaryChannel
    file "*"
    
    """
    fileName=\$(basename "${inputfile}")
    baseName="\${fileName%_vsearch.tsv}"
    
    perl ${params.ScriptSummaryVsearch} \
    ${inputfile} \
    \${baseName}_vsearchSummary.txt \
    ${params.maxhitsNOV} ${params.percentile} 

    cat \${baseName}_vsearchSummary.txt > out

    header="HEADER_vsearchSummaryReport_percIdent_${params.percidentity}_maxhits_${params.maxhitsNOV}_mincov_${params.mincov}"
    truncate -s 0 \${header}
    echo -e "command1 run: vsearch --usearch_global <inputfile> --db <database.fa>  --threads <threads>  --notrunclabels --id ${params.percidentity} --maxhits ${params.maxhitsNOV}  --blast6out <outputfile>" >> \${header}
    echo -e "command2 run: filter_vsearch_output <mincov> <outputfile.tsv> <filtered.outputfile.tsv> <outputfile.aln.txt> <filtered.outputfile.aln.txt>" >> \${header}
    echo -e "\n\nSample\tAmplicon\tEXTENDED_SEQID\tTOTAL_READS_IN_DEMULTIPLEXED_SAMPLE\tTOTAL_HITS\tMIN_percIdent\tMAX_percIdent\tAVG_percIdent\tMIN_alnLen\tMAX_alnLen\tAVG_alnLen\tMIN_Coverage\tMAX_Coverage\tAVG_Coverage\t${params.percentile}_PERCENTILE" >> \${header}

    """
} 


summaryChannel 
   .collectFile(name: file(params.out))
   .println { "Result saved to file: $it" }
   

workflow.onComplete {
    def subject = "[Task #${params.ticket}]  B_Allan Fluidigm_NOStitch : Fluidigm qctrim and vsearch clustering when no stitching of reads is involved "
    def recipient = '${params.email}'

    final_log = """
Pipeline execution summary
---------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}

    """
    
    ['mail', '-s', subject, recipient].execute() << "${final_log}\n${run_info}"
}

