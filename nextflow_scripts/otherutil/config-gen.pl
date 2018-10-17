#!/usr/bin/perl
###############################################
# Author: Gloria Rendon
# Date: Nov 2017
# program config-gen.pl
# Program that reads in a table of fluidigm data, one target per row, 
# and generates a separate nextflow configuration for each target
###############################################

use Getopt::Long;
use warnings;
use Carp;

my $date=localtime;
print "\n$date\nstarting config-gen.pl\n";

my $OutputPrefix="";
my $inputfile="";


$result=GetOptions ( 
		     "outprefix=s"                             => \$OutputPrefix,
		     "inputfile=s"                             => \$inputfile
 );


die("parsing failed\n")                                        if ! defined $result;
die("prefix must be specified\n")                              if $OutputPrefix eq "";
die("inputfile must be specified\n")                           if $inputfile eq "";


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
	next if $size < 24;
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

	# last param-name for maxhits, params.maxhitsV or params.maxhitsNOV, depends on type of target specified in $detline[10] 
	
	if ( $detline[10] eq "YES" || $detline[10] eq "Y") {
		$param = "params.maxhitsV";
		$value = "\'".$detline[24]."\'";
		print OUT "$param = $value\n";
	} else {
		$param = "params.maxhitsNOV";
		$value = "\'".$detline[24]."\'";
		print OUT "$param = $value\n";	
	}
	
	# done processing this row, Next one please
	close(OUT);
}

close(IN);

$date=localtime;
print "\n$date\ndone processing $inputfile. $totrec records processed.\n\nExiting config-gen.pl.\n";
