#########################################################################
# Author: Gloria Rendon
# Date: July 2016
# script to merge all amplicons into a single database. 
# Warnings: 
# Amplicons may come in several fasta files.
# Definition lines are not necessarily NCBI  or EMBL formatted, i.e. they may not have accession information that can be used 
# There could be duplicate sequences.
# All fasta files to be merged into a file must be in the same folder and have the same file extension 
#########################################################################

#########################################################################
# global variables
#########################################################################

$prefix="";         # the prefix of the fasta files in the folder
$inputfile="";      # the name of each input file with extension $prefix   
$type="";
%FASTA=();          # hash to hold all fasta sequences
$seqorder=0;        # a unique sequence identifier assigned to the sequences
$outfastaDB="";     # string to be used to name output fasta file
$outtableDB="";     # string to be used to name output table file
$count=0;           # counter of sequences read and processed

#########################################################################
# Main program starts here
#########################################################################

if ($#ARGV<2) {
	print "script to read all fasta files found in current folder; then merge into a single database of amplicons.\n\n";
	print "EXAMPLE: perl $0  InputFileExtension OutputDatabase.fa OutputDatabase.tsv \n\n";
	exit 0;
}
$prefix=shift;
$outfastaDB=shift;
$outtableDB=shift;

die("Required parameter: OutputDatabaseName \n") if $outfastaDB eq "";
die("Required parameter: OutputDatabaseTable \n") if $outtableDB eq "";
die("Required parameter: fastaFileExtension \n") if $prefix eq "";

open(OUTTABLEDB,">$outtableDB")  || die("Cannot open $outtableDB \n");
open(OUTFASTADB,">$outfastaDB")  || die("Cannot open $outfastaDB \n");

print "processing fasta files in the folder...\n\n";

$count = ReadSeqs();  #this function reads the files and populates a hash of sequences

$dup=$count - $seqorder;

print "$count sequences PROCESSED\t$seqorder sequences WRITTEN\t$dup duplicate sequences\n\n";
print "generating output files $outfastaDB  $outtableDB\n\n";

PrintDB();  # this function generates the output files


#reset data structures
%FASTA=();
close(OUTTABLEDB);
close(OUTFASTADB);

print "Done processing data. Exiting now\n\n";

exit 0;

#########################################################################
# Functions begin here
#########################################################################


sub PrintDB {
	print "Generating output files $outfastaDB $outtableDB\n\n";
	# outfastaDB is a global variable for output file with renamed sequences
	# outtableDB is a global variable for output file with sequences and the fasta files where they appear
	# FASTA is a global variable for a hash of fasta sequences

	# print header of outtableDB
	print OUTTABLEDB "SEQ_ORDER\tFILENAME(S\tDEFINITION_LINE";

	# print the rest of the file and also the database of sequences
	foreach $seqkey ( sort keys %FASTA ) {
	        $seq=$FASTA{$seqkey}{"SEQUENCE"};
	        $seqnum=$FASTA{$seqkey}{"SEQNUMBER"};
	        $seqfiles=$FASTA{$seqkey}{"SOURCE"};
	        #chomp $seqkey;
	        $seqkey =~ s/\n//g;
	        $seqfiles =~ s/\n//g;
		print OUTFASTADB "$seq";
		print OUTTABLEDB "\nSEQID$seqnum\t$seqfiles\t$seqkey";
	}

	#done. exiting now
	return 0;
}

sub ReadSeqs {
	print "inside ReadSeq function\n\nRead each fasta file in folder and populate FASTA, the hash of fasta sequences\n\n";
	$numfiles=0;
	$path=".";
	opendir( DIR, $path ) || die("Can't open $path: $!");
	print "reading dir $path....\n";

	while ( $inputfile = readdir( DIR ) ) {
		######################################################################
		print "next file in directory: $inputfile\n\nPerforming sanity checks\n\n";
		######################################################################
		
		#print "check if this is a hidden file\n";
		next if ($inputfile =~ /^\.{1,2}$/);  #skip hidden files

		#print "check if the file extension is $prefix\n";		
		next if ($inputfile !~ /$prefix/);    #skip files with extension other than prefix

		######################################################################
		#print "Next loop to read and process $inputfile\n\n";	
		######################################################################

		$numfiles++;
		local $/ = ">";                      # this is the record delimiter
		open(INFA,"<$inputfile") || die("Cannot open $inputfile \n");
		while (<INFA>) {
			$seq_all = $_;
			$seq_all =~ s/>//;           # to correct the parsing result by removing > from the end of the record
			$seq_all = ">".$seq_all;     # to correct the parsing result by prepending > to the record

			if ( $seq_all =~ m/>(.*)\n((.|\n)*)/ ) {
				$defline =$1;
				$seq=$2;
				# let's populate the FASTA hash
				if ( ! defined $FASTA{$defline} ) {
					#print "record does not exists: $defline\n";
					#create the record for this sequence in the FASTA hash
					$seqorder++;
					$newseq=">SEQID$seqorder $defline\n$seq";
					$FASTA{$defline}{"SEQUENCE"}=$newseq;
					$FASTA{$defline}{"SOURCE"}="$inputfile";
					$FASTA{$defline}{"SEQNUMBER"}=$seqorder;
			
				} else {
					#print "record already exists: $defline\n";
					#add this filename to the FASTA hash
					$listFiles=$FASTA{$defline}{"SOURCE"};
					$listFiles=$listFiles." ".$inputfile;
					$FASTA{$defline}{"SOURCE"}=$listFiles;
				}
			} else {
				print "seqnum=$count failed to parse $seq_all in file $inputfile\n";		
			}
			$count++;
			#end processing this sequences, next please
		}  # end while read inputfile
		close(INFA);
		######################################################################
		#print "End loop to read and process $inputfile\n\n";	
		######################################################################		
	} #end while readdir	

	closedir( DIR );
	print "done reading files. $numfiles files processed\n\n";
	#done. exiting now
	return $count;	
}


		
