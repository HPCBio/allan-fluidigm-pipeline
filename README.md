# Description

In a typical metagenomics analysis pipeline with fluidigm data, reads are compared to  sequences in a taxonomic database such as Silva in order to generate a taxonomic abundance table. In our case we need to prepare and use a custom database of tick pathogens in order to produce a presence/absence table with only the species included in that database. 

![Alt text](docs/workflow.jpg?raw=true "Workflow")


## old_scripts
In the `old_scripts` folder...

## nextflow_scripts
In the `nextflow_scripts` folder...

## test_data
In the `test_data` folder...

# Dependencies

This program expects the following tools/languages to be installed as modules and be available in your path:

- Nextflow        tested with version 0.27.3  ( download page https://github.com/nextflow-io/nextflow/releases?after=v0.29.0-RC1 )
- Trimmomatic     tested with version 0.36    ( download page https://github.com/timflutre/trimmomatic)
- PEAR            tested with version   0.9.8 (download page https://pear.php.net/package/PEAR/download)
- FASTX-Toolkit   tested with version  0.0.14 (download page http://hannonlab.cshl.edu/fastx_toolkit/download.html)
- vsearch         tested with version 2.4.3   (download page https://github.com/torognes/vsearch/releases)
- blast           tested with version BLAST+ 2.6.0  (download page ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
- Perl            tested with version 5.24.1  ( download page https://www.perl.org/ )

# Installation instructions

- Install all dependencies first. You may need to have root access to install some of these tools/languages on a cluster.
- Do not forget to launch the 'hello world' nextflow pipeline (as per https://www.nextflow.io/) to make sure it works fine.
- Make a copy of the  `nextflow_scripts` folder 

# Database preparation
- You can use the database we provide in the  `test_data` folder.
- To prepare your own database, put your seed sequences in a fasta file and run this program:

# More preparation steps

- <b>Prepare the configuration file</b>.
The fluidigm pipeline can process different kinds of amplicons such as V1, V3, V4. However, each amplicon has to be run separately because the analysis steps vary slightly. 
Prepare a configuration file for each amplicon. This configuration file is a list of parameter-value pairs that is used to specify read preparation, search and filtering options. For instance, if the reads overlap, then stitching the reads is a step that is performed by this pipeline, otherwise only R1 is used for the analysis. Some examples are provided in `nextflow_scripts/config`.

- <b>Prepare the raw reads</b>.
This pipeline expects paired-ended short reads that have been already demultiplexed.
The pipeline also expects that the filenames of the raw reads have this pattern: <i> amplicon </i> -  <i> sample </i> .R[1,2].fq
Put the demultiplexed raw reads that belong to the same amplicon in the same folder. 
It is important that you <i> do not mix </i> raw reads from different amplicons in the same folder. 
You can use any of the datasets provided in the  `test_data` folder.

# Running the program
To run the fluidigm pipeline type this command: <i> nextflow run -c config fluidigm-template-ballan-v0.3.nf  </i>

# Outputs
Nextflow generates two folders to keep track of execution progress. You can delete them once the execution ends successfully. They are called <i>.nextflow/ </i> and <i>work/ </i>

The actual results are placed in these folders: 
- <b>readprep/</b>  contains the results of QC, filter and trim of the raw reads. If the reads were stitched together, then PEAR results will be placed here too.
- <b>vsearchResults</b> contains the results of searches performed with VSEARCH, all subsequent filtering steps applied to each demultiplexed file as well as the final results of this step which is the file with this pattern 

<i> amplicon </i> * <i>__vsearchSummaryReport.txt</i>

This file is a tabulated table with these columns: Sample	Amplicon	EXTENDED_SEQID	TOTAL_READS_IN_DEMULTIPLEXED_SAMPLE	TOTAL_HITS	MIN_percIdent	MAX_percIdent	AVG_percIdent	MIN_alnLen	MAX_alnLen	AVG_alnLen	MIN_Coverage	MAX_Coverage	AVG_Coverage	5_PERCENTILE

Where:
- EXTENDED_SEQID is the sequence identifier of the hit that concatenates together the accession number and the species name
- TOTAL_READS_IN_DEMULTIPLEXED_SAMPLE is the total read count of the demultiplexed sample
- TOTAL_HITS is the read coverage for this hit
- MIN_percIdent	MAX_percIdent	AVG_percIdent are minimum, maximum and average percent identity of the reads for this hit
- MIN_alnLen	MAX_alnLen	AVG_alnLen are minimum, maximum and average alignment length of the reads for this hit
- MIN_Coverage	MAX_Coverage	AVG_Coverage are minimum, maximum and average coverage of the reads for this hit
- 5_PERCENTILE is the 5th percentile 

# Downstream analysis
This pipeline produces a presence/absence table with only the species present in the custome database. 
It is not an OTU table (table of taxonomic abundace) that could be further analyzed with tools such as QIIME.

# Citation

# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/




