#!/usr/bin/perl
###############################################
# Author: Gloria Rendon
# Date: Feb 2020
# program config-generator.pl
# Program that reads in a table of fluidigm data, one target per row, 
# and generates a separate nextflow configuration for each target
###############################################

use Getopt::Long;
use warnings;
use Carp;

######################################
# USAGE
######################################

my $USAGE =<<USAGE;

Script Name: config-generator.pl
Purpose: reads in a table of fluidigm data, one target per row, and generates a separate nextflow configuration for each target

USAGE:  perl config-generator.pl[options] 

options:
-h               to print this message
-outprefix       prefix of the config file
-inputfile       input file. It is a TAB delimited table with 27 columns
column1	params.Amplicon
column2	params.SeqLot
column3	params.email
column4	params.projectdir
column5	params.reads
column6	params.srcdir
column7	params.dbdir
column8	params.resultsdir
column9	params.readLen 
column10	params.skipStitch 
column11	params.isVtarget
column12	params.R1minlen 
column13	params.R1maxlen 
column14	params.R1adaptorLen 
column15	params.R2minlen 
column16	params.R2maxlen 
column17	params.R2adaptorLen 
column18	params.maxSeqLen 
column19	params.minSeqLen 
column20	params.minqual 
column21	params.minoverlap
column22	params.percidentity
column23	params.mincov
column24	params.percentile
column25	params.maxhits
column26	databaseFASTA
column27	databaseCSV

Example:

perl config-generator.pl  -outprefix nxf-data_2020_01 -inputfile Configuration_data_2020_01_gr.txt

USAGE

my $commandLine = $0 . " " . (join " ", @ARGV);
my $help;
my $OutputPrefix="";
my $inputfile="";



$result=GetOptions ( 
	"outprefix=s"          => \$OutputPrefix,
	"inputfile=s"          => \$inputfile,
	"h=s"                  =>\$help
 );

if ($commandLine =~ /-h/){
    print "$USAGE\n";
    exit;
}

die("ERROR: parsing failed\nType -h for help message\n\n")                    if ! defined $result;
die("ERROR: prefix must be specified\nType -h for help message\n\n")          if $OutputPrefix eq "";
die("ERROR: inputfile must be specified\nType -h for help message\n\n")       if $inputfile eq "";


###############################################
# read input file and generate one output file per row
###############################################
#
#the table has these columns
#
#col0  = params.Amplicon
#col1  = params.SeqLot
#col2  = params.email
#col3  = params.projectdir
#col4  = params.reads
#col5  = params.srcdir
#col6  = params.dbdir
#col7  = params.resultsdir	
#col8  = params.readLen 
#col9  = params.skipStitch
#col10 = params.isVtarget
#col11 = params.R1minlen 
#col12 = params.R1maxlen 
#col13 = params.R1adaptorLen 
#col14 = params.R2minlen 
#col15 = params.R2maxlen 
#col16 = params.R2adaptorLen 
#col17 = params.maxSeqLen 
#col18 = params.minSeqLen 
#col19 = params.minqual 
#col20 = params.minoverlap
#col21 = params.percidentity
#col22 = params.mincov
#col23 = params.percentile
#col24 = params.maxhits
#col25 = name of database.fasta
#col26 = name of database.csv
###############################################

my @header=();
my $totrec=0;
$date=localtime;
print "\n$date\nstart processing $inputfile\n";

open(IN,"<$inputfile") || die("$inputfile this file does not exist or it cannot be opened!\n");
while (my $line = <IN>) {
	chomp $line;
	print "processing $line";
	
	# processing header line
	if ($totrec < 1) { @header=split(/\t/,$line); $totrec++; next; }
	
	# processing next target	
	my @detline = split(/\t/,$line);
	my $size = @detline;
	next if $size < 26;
	$totrec++;
	
	# generating corresponding config file
	my $filename = $OutputPrefix."-".$detline[0].".conf";
	open(OUT,">$filename")  || die("$filename this file does not exist or it cannot be opened!\n");

	# the first values  must be enclosed in double quotes	
	for ($i=0; $i<=7; $i++) {
		$param = $header[$i];
		$value = "\"".$detline[$i]."\"";
		print OUT "$param = $value\n";
	}
	
	# the other values must be enclosed in single quotes	
	for ($i=8; $i<=23; $i++) {
		$param = $header[$i];
		$value = "\'".$detline[$i]."\'";
		print OUT "$param = $value\n";
	}

	# col24 param-name for maxhits, params.maxhitsV or params.maxhitsNOV, depends on type of target specified in $detline[10] 
	
	if ( $detline[10] eq "YES" || $detline[10] eq "Y") {
		$param = "params.maxhitsV";
		$value = "\'".$detline[24]."\'";
		print OUT "$param = $value\n";
	} else {
		$param = "params.maxhitsNOV";
		$value = "\'".$detline[24]."\'";
		print OUT "$param = $value\n";	
	}

  # col25 and col26 for params.DBfasta and params.DBtsv 
  $param = "params.DBfasta"; $value = "\"".$detline[6]."\/".$detline[25]."\"";
  print OUT "$param = $value\n";
	$param = "params.DBtsv";   $value = "\"".$detline[6]."\/".$detline[26]."\"";
  print OUT "$param = $value\n";

	# done processing this row, Next one please
	close(OUT);
}

close(IN);

$date=localtime;
print "\n$date\ndone processing $inputfile. $totrec records processed.\n\nExiting now\n";
