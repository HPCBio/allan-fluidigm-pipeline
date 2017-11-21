#!/usr/bin/env nextflow
/*
* USAGE: nextflow run fluidigm-pipe-generic-RrhiOB2F_RrhiOB2R-dataset20170515.nf -qs 6
* Where  "-qs" specifies the maximum # of processes that will run simultaneously.
*
* This pipeline performs fluidigm sequence QC-trim of paired reads, then processing and clustering with vsearch
* Fluidigm raw reads should be demultiplexed by amplicon-sample and arranged into separate folders as per fluidigm-sequencing-pipeline-V2
* This script creates hard links to data that exists in nextflow's work directory for everything
*/
/*
 
/*Some hacks to deal with dates*/

import java.text.SimpleDateFormat

// version
version = 0.1

//
// EDIT BLOCK STARTS HERE 
// The user can edit this section for input folder and tool parameters
//

/* Email to receive execution notifications */
params.email = "grendon@illinois.edu"

/*Path to project folder*/
projectPath = "/home/groups/hpcbio_shared/brianallan_group"

/*Name of inputs*/
params.Amplicon      = "RrhiOB2F_RrhiOB2R"         // The name of the Amplicon as it appears in the sequencing report and in the filenames of the reads
params.SeqLot        = "20170515"                        // This string that indicate the provenance of the Amplicon. Default: the date of the sequencing run


/*input channels */
rawDataPath = "$projectPath/newdata/rawReads/20170515_rawReads/RrhiOB2F_RrhiOB2R"

Channel
    .fromFilePairs("$rawDataPath/*_R{1,2}.fastq", flat: true)
    .ifEmpty {error "Cannot find any reads matching: ${rawDataPath}"}
    .set {readFilesTrim}

readsCh = Channel.fromPath("$rawDataPath/*_R1.fastq")
              .ifEmpty { error "Cannot find any FASTQ file in: ${rawDataPath}" }

/*folder paths for results*/

timestamp      = new SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())
readprepPath   = "$projectPath/newresults/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/readprep"
vsearchPath    = "$projectPath/newresults/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/vsearchUnfiltered"
summaryPath    = "$projectPath/newresults/" + timestamp + "_seqLot${params.SeqLot}_${params.Amplicon}/vsearchFinalResults"

/*Name of outputs*/
params.out           = "${summaryPath}/${params.Amplicon}_vsearchSummaryReport.txt" // The name of the output report


/*Variables for read preparation. parameters for R2 will be ignored if stitching is not performed*/
params.skipStitch    = 'YES'                       // required. possible values YES NO
params.readLen       = '250'                       // read length
params.R1minlen      = '170'                       // minimum sequence length after qctrim of R1
params.R1maxlen      = '230'                       // maximum sequence length of R1 after qctrim (deafult value= readlen -adaptorlen)
params.R1adaptorLen  = '20'                        // length of adaptor for R1. This value comes from the sequencing report
params.R2minlen      = '170'                       // minimum sequence length of R2 after qctrim 
params.R2maxlen      = '230'                       // maximum sequence length of R2 after qctrim (deafult value= readlen -adaptorlen)   
params.R2adaptorLen  = '20'                        // length of adaptor for R2. This value comes from the sequencing report
params.maxSeqLen     = '230'                       // maximum sequence length of STITCHED read
params.minSeqLen     = '170'                       // minimum sequence length of STITCHED read
params.minqual       = '30'                        // minimum base quality. The range goes from 0-40                        
params.minoverlap    = '20'                        // default=20


/*Variables for clustering with VSEARCH*/
params.percidentity  = '0.90'                      // required. percent identity cutoff       
params.maxhitsV      = '50'                        // required. maximum number of hits to report per query; applicable to all V amplicons              
params.maxhitsNOV    = '10'                        // required. maximum number of hits to report per query; applicable to all NON-V amplicons 
params.mincov        = '50'                        // required. minimum coverage cutoff to use when post-processing hits   
params.percentile    = '5'                         // required. percentile value

/* Path for internal DB, programs and scripts*/
params.ScriptFilterVsearch  = "$projectPath/newsrc/nextflow/filter_vsearch_output.pl" 
params.ScriptSummaryVsearch = "$projectPath/newsrc/nextflow/vsearch_calcSummary.pl" 
params.DBfasta              = "/home/groups/hpcbio_shared/brianallan_group/newdata/database/2017-11/DB_version2_combined_and_filtered_2017.fasta"
params.DBtsv                = "/home/groups/hpcbio_shared/brianallan_group/newdata/database/2017-11/DB_version2_combined_and_filtered_2017.tsv"
params.adaptors             = "$projectPath/newsrc/nextflow/adaptors.fa"

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
pearMemory    = '5'
pearCPU       = '4'

/*Module versions*/
pearMod       = 'PEAR/0.9.8-IGB-gcc-4.9.4'
fastxMod      = 'FASTX-Toolkit/0.0.14-IGB-gcc-4.9.4'
vsearchMod    = 'VSEARCH/2.4.3-IGB-gcc-4.9.4'
trimVersion   = '0.36'
trimMod       = "Trimmomatic/0.36-Java-1.8.0_121"
trimHome      = '\$EBROOTTRIMMOMATIC'
trimJar       = "${trimHome}/trimmomatic-${trimVersion}.jar"

/* trim options for stitching case */
R1trimOptions   = " HEADCROP:${params.R1adaptorLen} "
R2trimOptions   = " HEADCROP:${params.R2adaptorLen} "
stitchOptions   = " -v ${params.minoverlap} -t ${params.minSeqLen} -m ${params.maxSeqLen} -q ${params.minqual} "

/* trim options for non-stitching case */
trimOptions     = " HEADCROP:${params.R1adaptorLen} ILLUMINACLIP:${params.adaptors}:1:10:7 LEADING:${params.minqual} TRAILING:${params.minqual} SLIDINGWINDOW:4:20 CROP:${params.R1maxlen} "


// Redmine Ticket must be set for project mgmt tracking
params.ticket = 7000
if (params.ticket == 0) exit 1, "Must set Redmine ticket for pipeline summary to be sent"


// Credit to Phil Ewels for this little bit of self-documentation
run_info = """
          
====================================
 B_Allan FLUIDIGM-RrhiOB2F_RrhiOB2R : Fluidigm qctrim and vsearch clustering for RrhiOB2F_RrhiOB2R dataset-20170515 
====================================

Pipeline ver : ${version}
Amplicon     : ${params.Amplicon}
SequencingLot: ${params.SeqLot}
Reads        : ${rawDataPath}
Database     : ${params.DBfasta}
VSEARCH      : ${vsearchMod}
Fastx-Toolkit: ${fastxMod}
Trimmomatic  : ${trimMod}
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
* Step 1. Trimmomatic, PEAR, concat, convert2fasta using R1 and R2
*/

process qctrimNstitch {
    executor   'slurm'
    cpus       trimCPU
    queue      myQueue
    memory     "$trimMemory GB"
    module     trimMod, pearMod, fastxMod
    publishDir readprepPath, mode: 'link'

    input:
    set val(id), file(read1), file(read2) from readFilesTrim

    output:
    file '*fasta' into trimmedFastaData1
    file "*"

    when:
    params.skipStitch    =~ /NO/    

    script:
    """ 
    #!/bin/bash
    
    fileName=\$(basename "${read1}")
    baseName="\${fileName%_*_R1.fastq}"

    truncate -s 0 \${baseName}_SearchReady.fasta

    
    java -jar ${trimJar} SE -threads $trimCPU -phred33 $read1 \${baseName}.qctrim.R1.fastq ${R1trimOptions} 2> \${baseName}.R1.trimrun.log 
    java -jar ${trimJar} SE -threads $trimCPU -phred33 $read2 \${baseName}.qctrim.R2.fastq ${R2trimOptions} 2> \${baseName}.R2.trimrun.log    

    
    pear ${stitchOptions} -f \${baseName}.qctrim.R1.fastq -r \${baseName}.qctrim.R2.fastq -o \${baseName}.stitch


    cat \${baseName}.stitch.assembled.fastq \${baseName}.stitch.unassembled.forward.fastq >> \${baseName}.concat.fastq

    if [ -s \${baseName}.concat.fastq ]
    then
    	fastq_to_fasta -v -i \${baseName}.concat.fastq -o \${baseName}_SearchReady.fasta
    fi
    
    """
}

/*
* Step 1. Trimmomatic and convert2fasta on R1 only
*/

process qctrimOnly {
    executor   'slurm'
    cpus       trimCPU
    queue      myQueue
    memory     "$trimMemory GB"
    module     trimMod, fastxMod
    publishDir readprepPath, mode: 'link'

    input:
    file read from readsCh 

    output:
    file '*fasta' into trimmedFastaData2
    file "*"

    when:
    params.skipStitch    =~ /YES/    

    script:
    """
    #!/bin/bash
    
    fileName=\$(basename "${read}")
    baseName="\${fileName%_*_R1.fastq}"

    truncate -s 0 \${baseName}_SearchReady.fasta      

    java -jar ${trimJar} SE -threads $trimCPU -phred33 \
    ${read} \${baseName}.trimmed.fastq ${trimOptions} 2> \${baseName}.trimrun.log 

    if [ -s \${baseName}.trimmed.fastq ]
    then
       fastq_to_fasta -v -i \${baseName}.trimmed.fastq -o \${baseName}_SearchReady.fasta
    fi
    
    """
}

/*
* Step 2. Run VSEARCH for stitched reads
*/

process runVsearch {
    executor   'slurm'
    cpus       vsearchCPU
    queue      myQueue
    memory     "$vsearchMemory GB"
    module     vsearchMod
    publishDir vsearchPath, mode: 'link'

    input:
    file inputfile from trimmedFastaData1


    output:
    file "*_vsearch.tsv" into vsearchChannel1
    file "*" 
    
    
    """   
    #!/bin/bash
    
    fileName=\$(basename "${inputfile}")
    baseName="\${fileName%_SearchReady.fasta}"
    
    truncate -s 0 \${baseName}_vsearch.tsv

    if [ -s ${inputfile} ]
    then
       vsearch --usearch_global ${inputfile} \
       --db ${params.DBfasta} \
       --notrunclabels    --threads ${vsearchCPU} \
       --id        ${params.percidentity} \
       --maxhits   ${params.maxhitsV}  \
       --alnout    \${baseName}_vsearch_prefiltered_aln.txt \
       --blast6out \${baseName}_vsearch_prefiltered.tsv
    
 
       perl ${params.ScriptFilterVsearch}   ${params.mincov} \
       \${baseName}_vsearch_prefiltered.tsv \
       \${baseName}_vsearch.tsv \
       \${baseName}_vsearch_prefiltered_aln.txt \
       \${baseName}_vsearch.aln.txt
    
    fi
    
    """

}

/*
* Step 2. Run VSEARCH for non-stitched
*/

process runVsearch {
    executor   'slurm'
    cpus       vsearchCPU
    queue      myQueue
    memory     "$vsearchMemory GB"
    module     vsearchMod
    publishDir vsearchPath, mode: 'link'

    input:
    file inputfile from trimmedFastaData2


    output:
    file "*_vsearch.tsv" into vsearchChannel2
    file "*" 
    
    
    """   
    #!/bin/bash
    
    fileName=\$(basename "${inputfile}")
    baseName="\${fileName%_SearchReady.fasta}"
    
    truncate -s 0 \${baseName}_vsearch.tsv

    if [ -s ${inputfile} ]
    then
       vsearch --usearch_global ${inputfile} \
       --db ${params.DBfasta} \
       --notrunclabels    --threads ${vsearchCPU} \
       --id        ${params.percidentity} \
       --maxhits   ${params.maxhitsV}  \
       --alnout    \${baseName}_vsearch_prefiltered_aln.txt \
       --blast6out \${baseName}_vsearch_prefiltered.tsv
    
 
       perl ${params.ScriptFilterVsearch}   ${params.mincov} \
       \${baseName}_vsearch_prefiltered.tsv \
       \${baseName}_vsearch.tsv \
       \${baseName}_vsearch_prefiltered_aln.txt \
       \${baseName}_vsearch.aln.txt
    
    fi
    
    """

}

/*
* Step 3. Run VSEARCH-Summary for stitched reads
*/

process runSummary {
    executor   'slurm'
    cpus       1
    queue      myQueue
    memory     "4 GB"
    publishDir summaryPath, mode: 'link'
    
    input:
    file inputfile from vsearchChannel1

    output:
    file 'out' into summaryChannel1
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

/*
* Step 3. Run VSEARCH-Summary for non-stitched reads
*/

process runSummary {
    executor   'slurm'
    cpus       1
    queue      myQueue
    memory     "4 GB"
    publishDir summaryPath, mode: 'link'
    
    input:
    file inputfile from vsearchChannel2

    output:
    file 'out' into summaryChannel2
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

summaryChannel1 
   .collectFile(name: file(params.out))
   .println { "Result saved to file: $it" }
   
summaryChannel2 
   .collectFile(name: file(params.out))
   .println { "Result saved to file: $it" }


workflow.onComplete {
    def subject = "[Task #${params.ticket}]  B_Allan FLUIDIGM-RrhiOB2F_RrhiOB2R : Fluidigm qctrim and vsearch clustering for RrhiOB2F_RrhiOB2R dataset-20170515"
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

