# Author: Gloria Rendon
# Date: March 2016
#
# Script to remove lines w low_coverage alignments from  the vsearch reports -- written in blast_output_tabular format
# and create a new file with  filtered results
#
# to run this program type: filter_vsearch_output.pl <min_coverage> <vsearch_report.infile.tsv> <vsearch_report.infile.aln.txt> <vsearch_report.outputfile>
#
#######################
# infile and outfile are vsearch reports; i.e. files in tab-delimited format with columns in this order:
#
# col0=queryId
# col1=subjectId
# col2=percIdentity
# col3=alnLength
# col4=mismatchCount
# col5=gapOpenCount
# col6=queryStart
# col7=queryEnd
# col8=subjectStart
# col9=subjectEnd
# col10=eVal
# col11=bitScore
#######################
# outfile contains all the columns of infile PLUS a column for coverage:
#
# col0=queryId
# col1=subjectId
# col2=percIdentity
# col3=alnLength
# col4=mismatchCount
# col5=gapOpenCount
# col6=queryStart
# col7=queryEnd
# col8=subjectStart
# col9=subjectEnd
# col10=eVal
# col11=bitScore
# col12=coverage
#######################
# the formula for applying the filtering is as follows:
#
# qlen= queryEnd - queryStart
# slen= subjectEnd - subjectStart
# coverage= ( alnLength * 2) / ( qlen + slen )
# if coverage < mincov then the line is filtered out and not written to the output file
#######################

########################
# global variables
########################

$invsearchfile="";
$inalnfile="";
$outvsearchfile="";
$outalnfile="";
$mincov="";
$tcount=0;
$delcount=0;
%REMOVED=();

########################
# reading input from stdin and sanity check
########################

if ($#ARGV<4) {
	print "script to filter vsearch reports using min_coverage as the cutoff\n\n";
	print "EXAMPLE: $0 min_coverage input.vsearch.tsv output.vsearch.tsv input.vsearch.aln.txt output.vsearch.aln.txt\n\n";
	exit 0;
}

$mincov = shift;
$invsearchfile= shift;
$outvsearchfile= shift;
$inalnfile= shift;
$outalnfile= shift;

print "values read\nmin_coverage=$mincov\n";
print "input.vsearch.tsv=$invsearchfile\noutput.vsearch.tsv=$outvsearchfile\n";
print "input.vsearch.aln.txt=$inalnfile\noutput.vsearch.aln.txt=$outalnfile\n\n";

die("No output vsearch filename was specified\n") if $outvsearchfile eq "";
die("No input vsearch filename was specified\n") if $invsearchfile eq "";
die("No input aln filename was specified\n") if $inalnfile eq "";
die("No output aln filename was specified\n") if $outalnfile eq "";
die("No value for min_coverage was specified\n") if $mincov eq "";



#######################
# reading infile, one line at a time, and applying the filtering procedure
#######################
print "processing $invsearchfile and generating $outvsearchfile\n\n";

open(OUT,">$outvsearchfile") || die("Cannot open $outvsearchfile !\n");
open(IN,"<$invsearchfile") || die("$invsearchfile this file does not exist or it cannot be opened!\n");

while (my $line = <IN>) {
    chomp $line;
    @vals = split(/\t/,$line);
    $alnLength=$vals[3];
    $queryEnd=$vals[7];
    $queryStart=$vals[6];    
    $subjectEnd=$vals[9];
    $subjectStart=$vals[8];   
    $qlen= $queryEnd - $queryStart;
    $slen= $subjectEnd - $subjectStart;
    $coverage = int( ( $alnLength * 2 * 100 ) / ( $qlen + $slen ) );
    if ( $coverage < $mincov ) {
       #print "\[$tcount\]: alnlen=$alnLength qlen=$qlen slen=$slen coverage=$coverage mincov=$mincov\tthis record will be filtered out of output report\n";
       #add this record to hash of REMOVED alignments
       $id=$vals[0]."\t".$vals[1];
       $REMOVED{$id}=1;
       $delcount++;       
    } else {
       #print "\[$tcount\]: alnlen=$alnLength qlen=$qlen slen=$slen coverage=$coverage mincov=$mincov\tthis record will be included in the output report\n";
       print OUT "$line\t$coverage\n";
    }
    $tcount++;
}

print "done processing $invsearchfile\ttotal records processed:$tcount\ttotal records filtered out:$delcount\n\n";
close(OUT);
close(IN);

#######################
# reading alnfile, several lines at a time for each alignment result
#######################
print "processing $inalnfile and generating $outalnfile\n\n";
$numaln=0;
$numdel=0;
open(OUTALN,">$outalnfile") || die("Cannot open $outalnfile !\n");
open(ALN,"<$inalnfile") || die("$inalnfile this file does not exist or it cannot be opened!\n");
$/ = "\nQuery >";  #this is the record separator
while ($record = <ALN>) {
   #print "\nparse rec=$record";
   if ( $record =~ m/((.|\n)*)\n Query(.*)nt >(.*)\nTarget(.*)nt >(.*)\n/ ) {
       $query=$4; $target=$6;
       #print "\nquery=$query\ttarget=$target";
       $this_id=$query."\t".$target;
       if ( defined $REMOVED{$this_id} ) {
          print "this rec\[$this_id\] was removed from the tsv file and is marked for removal from the output aln file too\n";
          $numdel++;
       } else {
          print "this rec\[$this_id\] was NOT removed from the tsv file and will be included in the output aln file\n";       
         $record =~ s/\n\nQuery >/\n/;
         $record = "\nQuery >".$record;
         print OUTALN "$record";
       }   
   }
   $numaln++;    
}


print "done processing $infile\ttotal records processed:$numaln\ttotal records filtered out:$numdel\n\n";
close(OUTALN);
close(ALN);
%REMOVED=();
exit 0;