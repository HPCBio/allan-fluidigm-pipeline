#!/usr/bin/env perl

#version 0.1

##############################################################################################
# Program: filter_blastout.pl
#
# Description: program that reads in a blast report, filters out redundant or weak hits and generates two files. FASTA and TSV.
#              The output files that this program produces are the input database files for running the fluidigm pipeline
#
# Parameters:  This program needs these parameters
# -inputfa  <name of fasta file with seed sequences>
# -blastin  <name of blast report>
# -outfasta <name of output database file in fasta format>
# -outtab   <name of output database file in TSV format>
##############################################################################################

use Getopt::Long;

%SEQ=();
my $totseq=0;
my $totkept=0;
my $totfiltered=0;
my $seqid=0;

my $path=`pwd`;
my $date=localtime;
print "\n$date\nstarting program filter_blastout.pl from $path\n";

$result=GetOptions ("inputfa=s"    => \$inputfa,
		    "blastin=s"    => \$blastin,
		    "outfasta=s"   => \$outfasta,
		    "outab=s"      => \$outTab);

die("parsing failed\n")        if ! defined $result;

open(INFASTA, "<$inputfa")     or die "cannot read from inputfa $inputfa\n";
close(INFASTA);

open(OUTAB, ">$outTab")        or die "cannot write to outab $outTab  \n";
open(OUTFA, ">$outfasta")      or die "cannot write to outfasta $outfasta\n";
open(BLASTREPORT, "<$blastin") or die "cannot read from blastin $blastin\n";

##############################################################################################
# First part: read report and filter hits using two criteria
# condition1 short alignment: qcovs < 90
# condition2 duplicate sequence: when two or more hits have the same sacc sstart ssend 
#
##############################################################################################
# these are the columns expected in the input file
#col0=qseqid   this is the header of the seed sequence. This field must be formatted as follows: amplicon;species;taxidnnnn
#col1=sseqid
#col2=pident
#col3=length
#col4=mismatch
#col5=gapopen
#col6=qstart
#col7=qend
#col8=sstart
#col9=send
#col10=evalue
#col11=bitscore
#col12=qacc
#col13=sacc
#col14=qlen
#col15=slen
#col16=qcovs
#col17=staxids
#col18=sscinames
#col19=sseq
##############################################################################################


$date=localtime;
print "\n$date\nreading input file...\n";

while (my $line=<BLASTREPORT>){

  	chomp $line;
  	$totseq++;
  	print "processing row\[$totseq\]\n";
	@det= split(/\t/, $line);                                    #parse line
	@qseqid = split(/;/,$det[0]);                                #parse qseqid
	my $hitident = $det[13]."_".$det[8]."_".$det[9];             #make hit_identity by combining 3 fields
	
	# check conditions for filtering this hit
	# cond1: check if alignment is too short

	if ( $det[16] < 90 ) {
		print "row\[$totseq\]; hit with low coverage, removing it\n";
		$totfiltered++;
		next;
	}

	# cond2: check if it is duplicate hit
	if ( defined $SEQ{$hitident} ) {
		print "row\[$totseq\]; duplicate hit, removing it\n";
		$totfiltered++;
		my $val=$SEQ{$hitident}{"amplicon"};
		$val=$val.";".$qseqid[0] if $qseqid[0] !~ $val;      #to remove duplicates in the same aplicons
		$SEQ{$hitident}{"amplicon"}=$val;
	} else {
		# this hit is good, we keep it
		$totkept++;
		$SEQ{$hitident}{"amplicon"}=$qseqid[0];
		$SEQ{$hitident}{"detline"}=$line;	
	}
	
	# we're done with this hit; next one please
}

close(BLASTREPORT);

$date=localtime;
print "\n$date\nfinished reading input file.\nTotal lines processed $totseq\nTotal hits removed $totfiltered\nTotal hits kept $totkept\n\n\n";
print "generating output files...";

##############################################################################################
#
# Second part: generate all output files: FASTA and TAB
#
##############################################################################################


foreach my $key ( sort keys %SEQ ) {
	
	$seqid++;
	my $uniqueid="seqid".$seqid;
	my $amplicons=$SEQ{$key}{"amplicon"};
	
	my $line = $SEQ{$key}{"detline"};
	my @detline= split(/\t/,$line);

	my $sacc   = $detline[13];
	my $sstart = $detline[8];
	my $send   = $detline[9];	
	my $species= $detline[18];
	my $taxids = $detline[17];
	my $seqbody= $detline[19];
	$seqbody =~ s/-//g;                                          #to remove gaps from the sequence
	
	my $seqheader= ">".$uniqueid."|".$species."|taxids ".$taxids."|".$amplicons."|accession ".$sacc."|".$sstart."|".$send;
	my $rowdet   = $uniqueid."\t".$species."\ttaxids ".$taxids."\t".$amplicons."\taccession ".$sacc."\t".$sstart."\t".$send;	

	print OUTFA "$seqheader\n$seqbody\n";
	print OUTAB "$rowdet\n";
}

$date=localtime;
print "\n$date\nfinished writing filtered blast results to database files.\ntotal sequences in database $seqid\n\n";

##############################################################################################
#
# Third part: append seed sequences to both database files
#
##############################################################################################
#
# Expected format of the header line for the seed sequences is three values separated by ; like this:
# >amplicon_name;species_name;taxidnnnnn
#
# example
# >16S8FE3;Ehrlichia_canis;taxid269484

$/='>'; # set record delimiter
open(INFASTA, "<$inputfa")     or die "cannot read from inputfa $inputfa\n";

while (my $sequence = <INFASTA>) {
   	if ($sequence =~ /(.*)\n((.|\n)*)/) {
          	# parse fasta record
		my $descline=$1;
          	my $seqbody=$2;        
          	$seqbody =~ s/>//g;                  # to cleanup the body
          	my @header=split(/;/,$descline);     # to parse the header line

		print "adding $descline\n";
		
		# let us put together the output

		$seqid++;
		my $uniqueid="seqid".$seqid;
		my $amplicons=$header[0];
		my $species=$header[1];
		my $taxids=$header[2];
		$taxids =~ s/taxid//;

		my $seqheader= ">".$uniqueid."|".$species."|taxids ".$taxids."|".$amplicons;
		my $rowdet   = $uniqueid."\t".$species."\ttaxids ".$taxids."\t".$amplicons;	

		print OUTFA "$seqheader\n$seqbody\n";
		print OUTAB "$rowdet\n";
		
	}
}

close(INFASTA);
close(OUTFA);
close(OUTAB);

$date=localtime;
print "\n$date\nfinished adding seed sequences to database files.\ntotal sequences in database $seqid\n\nDone xiting now";
exit 0;
