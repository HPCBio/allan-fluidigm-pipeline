# Description

General description goes here...
## old_scripts
In the `old_scripts` folder...

## nextflow_scripts
In the `nextflow_scripts` folder...

## test_data
In the `test_data` folder...

# Dependencies

This program expects the following tools/languages to be installed as modules and be available in your path:

- Nextflow        tested with version 0.27.3 ( download page https://github.com/nextflow-io/nextflow/releases?after=v0.29.0-RC1 )
- Trimmomatic     tested with version 0.36 ( download page https://github.com/timflutre/trimmomatic)
- PEAR            tested with version   0.9.8 (download page https://pear.php.net/package/PEAR/download)
- FASTX-Toolkit   tested with version  0.0.14 (download page http://hannonlab.cshl.edu/fastx_toolkit/download.html)
- vsearch         tested with version 2.4.3 (download page https://github.com/torognes/vsearch/releases)
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

Prepare the configuration file.
The fluidigm pipeline can process different kinds of amplicons such as V1, V3, V4. However, each amplicon has to be run separately because the analysis steps vary slightly. 
Prepare a configuration file for each amplicon. This configuration file is a list of parameter-value pairs that is used to specify read preparation, search and filtering options.  Some examples are provided in `nextflow_scripts/config`.

Prepare the rew reads.
- This pipeline expects paired-ended short reads that have been already demultiplexed.
- Put the raw reads in one folder and give the folder the exact name of the amplicon. 
- You can use any of the datasets provided in the  `test_data` folder.

# Running the program
To run the fluidigm pipeline type this command: <i> nextflow run -c config fluidigm-template-ballan-v0.3.nf  </i>

# Outputs

# Downstream analysis

# Citation

# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/




