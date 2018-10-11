# Author: Gloria Rendon
# Date: April 2011

# Script to clean up a fasta file of DNA or RNA; 
# replacing residues [-O*JRY] for N everywhere in the sequence except in the description line

use warnings;
use strict;
use Carp;


###############################################
# getting the parms from STDIN and validating input
###############################################
my $infile="";
my $outfile="";

if ($#ARGV<1) {
        print "program that replaces - O  J Y R K for N in a fasta file of DNA or RNA sequences\n\n";
	print "EXAMPLE:$0 infile.fasta out.file.fasta\n\n";
	exit 0;
}

$infile = shift;
$outfile = shift;
die("Input file name must be specified\n") if $infile eq "";
die("Output file name must be specified\n") if $outfile eq "";

open(OUT,">$outfile") || die("Cannot open $outfile !\n");
open(FASTA,"<$infile") || die("$infile this file does not exist or it cannot be opened!\n");
$/='>'; # set record delimiter
while (<FASTA>) {
   if (/(.*)\n((.|\n)*)/) 
        {
          my $descline=$1;
          my $seq=$2;
          $seq =~ s/ //g;
          $seq =~ s/[^AGCT]/N/g;
          $seq =~ s/NNN$//;
          print OUT ">$descline\n$seq\n";
	}
}
close(FASTA);
close(OUT);
exit 0;
