#!/usr/bin/env nextflow
/*
 B_ALLAN PREPARE DATABASE
This pipeline runs local blast using a file of seed sequences against NCBI-NT and then it filters the resulting blast report of hits
Finally it produces TWO output files: one in fasta format and one in tabular format. Both output files will be needed by the fluidigm analysis pipeline
*/

// version
version = 0.1


// Credit to Phil Ewels for this segment to generate the help message at the command-line

def helpMessage() {
    log.info"""
    =========================================
     B_ALLAN PREPARE DATABASE v${version} 
    =========================================

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run prepare-database-ballan-v1.nf
    
    To override existing values from the command line, please use these parameters:
    
    --blastdbdir                  Path to directory with blast-formatted NCBI-NT database
    --queryFile                   File of seed sequences in fasta format
    --outPrefix                   Prefix of the output files
    --outdir                      Path of output directory
    
    Other optional parameters:
        -with-dag                     Graphic of the execution report in PDF format. Valid value:a file name
        -with-timeline                Timeline of all processes executed durin the execution. Report in HTML format. Valid value:a file name 
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

params.email         = "grendon@illinois.edu"
params.projectdir    = "/home/a-m/grendon/Fluidigm"
params.srcdir        = "/home/a-m/grendon/Fluidigm/current_src/dbutil"
params.blastdbdir    = "/home/mirror/ncbi/blastdb/20170702/nt"
params.queryFile     = "/home/a-m/grendon/Fluidigm/current_src/dbutil/Seed_Seqs_for_database_version2_2017.fasta"
params.outPrefix     = "DB_version2_2017"

/* Path for folders with temp and final results  */
params.outdir        = "${params.projectdir}/results/makedatabase"
params.tmpdir        = "${params.projectdir}/results/makedatabase_tmp"


/*Biocluster2 options. List memory in gigabytes like in suggestions below*/
myQueue        = 'normal'
blastMemory    = '150'
blastCPU       = '4'           
postMemory     = '2'
postCPU        = '1'

/*Biocluster2 module names and versions*/

blastMod      = 'BLAST+/2.6.0-IGB-gcc-4.9.4'

// variables for files
dbFastaFile      = file(params.blastdbdir)
queryFile        = file(params.queryFile)
seedFile         = file(params.queryFile)

// Sanity check
if (!params.outdir)          exit 1, "Must set --outdir"
if (!params.blastdbdir)      exit 1, "Must set --blastdbdir"
if (!params.outPrefix)       exit 1, "Must set --outPrefix"
if (!params.queryFile)       exit 1, "Must set --queryFile"

// Credit to Phil Ewels for this little bit of self-documentation
run_info = """
          
====================================
B_ALLAN PREPARE DATABASE v${version} 
====================================

Pipeline version            : ${version}                    
Input file w Seed Seq       : ${params.queryFile}
Output prefix               : ${params.outPrefix}
Blast                       : ${blastMod}

Script dir                  : ${params.srcdir}
dir w Blast-formatted DB    : ${params.blastdbdir}
Output dir                  : ${params.outdir} 
Project dir                 : ${params.projectdir}
Current user                : $USER
Current path                : $PWD

====================================

"""

log.info run_info


/*
* Step 1. Blast step
*/

process Blast {
    executor             'slurm'
    queue                myQueue
    memory               "$blastMemory GB"
    cpus                 blastCPU
    module               blastMod
    publishDir           params.tmpdir, mode: 'copy'

    input:
    file DB              from dbFastaFile
    file query           from queryFile

    output:
    file 'OUT.TAB'       into blastChannel
    file "*"

    script:
    
    """
    echo run blastn
  
    blastn -query $query -db $DB -num_threads $blastCPU -max_hsps 1 -perc_identity 90 -max_target_seqs 100 -outfmt 11 > OUT.ASN

    echo produce blast report in tabular format

    blast_formatter -archive OUT.ASN -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qacc sacc qlen slen qcovs staxids sscinames sseq" > OUT.TAB


    """
  
}  

process filterBlast {
    executor             'slurm'
    queue                myQueue
    memory               "$postMemory GB"
    cpus                 postCPU
    publishDir           params.outdir, mode: 'copy'

    input:
    file TAB             from blastChannel
    file seed            from seedFile

    output:
    file "*"

    script:
    
    """

    echo filter and format
    
    perl ${params.srcdir}/filter_blastout.pl -inputfa $seed -blastin $TAB -outfasta ${params.outPrefix}.fasta -outab ${params.outPrefix}.tab


    """
}


workflow.onComplete {
    def subject = "B_ALLAN PREPARE DATABASE v${version}"

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