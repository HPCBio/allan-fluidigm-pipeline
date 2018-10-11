# Author: Gloria Rendon
# Date: August 2016

# Script to parse the vsearch reports written in blast_output_tabular format
# an create a tab delimited file with summary results.


print "\n\n##############################################################################";
print "\n############   STARTS VSEARCH_CALCSUMMARY.SH           #########";
print "\n##############################################################################\n\n"; 



########################
# global variables
########################


$inputfile="";
$outputfile="";
$percentile="";
$maxhits="";
$suffix="vsearch.tsv";

my %PERCARRAY=();         # hash for percentId values for a hit
my %COVARRAY=();          # hash for coverage scores for a hit
my %ALNARRAY=();          # hash for alignment values for a hit
my %SEENPERC=();          # percentId values for the entire report
my %WRITTEN=();
my $numreads=0;           # number of reads for the entire report 


my $samplename="";
my $ampliconname="";

print "\n##############################################################################";
print "\n############    reading input from stdin and sanity check\n";
print "\n##############################################################################\n";

#if ($#ARGV<3) {
#	print "script that parses all vsearch reports for amplicons, then generates a summary of results.\n\n";
#	print "Each amplicon has its own folder of results\n";
#	print "EXAMPLE: $0 <vsearchIn> <vsearchOut> <maxhits> <percentile>\n\n";
#	exit 0;
#}


$inputfile=shift;          #filename of vsearch input file
$outputfile = shift;       #filename of the vsearch summary report
$maxhits = shift;          #threshold for maximum number of lines to report per sample-amplicon pair
$percentile = shift;       #only one percentile value will be caculated and reported


die("No inputfile was specified\n")            if $inputfile  eq "";
die("No outputfile  was specified\n")          if $outputfile eq "";
die("No value for maxhits was specified\n")    if $maxhits    eq "";
die("No value for percentile was specified\n") if $percentile eq "";

print "\n##############################################################################";
print "\n############   input arguments have been parsed successfully\n";
print "\n##############################################################################";
print "\n############   Sanity check: return if inputfile is empty\n";
print "\n##############################################################################";



open (OUTSUM,">>$outputfile") || die("cannot open outputfile $outputfile"); #we need to open it to append because it already has a header

print "\n############   prepare outputfile= $outputfile\n";

		

########################
# Processing inputfile
########################


print "\n############   processing $inputfile\n";

if ( $inputfile =~ m/(.*)-(.*)_vsearch.tsv/ ) { $samplename=$1; $ampliconname=$2; }

print "\n############   inputfile=$inputfile sample=$samplename amplicon=$ampliconname\n";
	
open(VSEARCH,"<$inputfile") || die("cannot read file $inputfile\n");

	
#######################
#
# vsearch report is a blast-like output text file in tab-delimited format with columns in this order:
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
# col12=coverage  THIS FIELD IS NOT PART OF THE ORIGINAL VSEARCH OUTPUT. IT IS A CALCULATED FIELD. 
#######################
	
	while ( my $vline = <VSEARCH> ) {
		#print "processing $vline";
		chomp $vline;
		@vcols=split(/\t/,$vline); #parsing this line, see description of columns above

		# here we simply populate the hashes with stats, but we do not write summaries yet
                if ( defined $PERCARRAY{$vcols[1]} ) {
			#this hit already exists, update values of hashes w stats 
  			
  			push @{$PERCARRAY{$vcols[1]}}, $vcols[2];  # grow the the array of percent identities
  			push @{$COVARRAY{$vcols[1]}},  $vcols[12]; # grow the the array of coverage scores
  			push @{$ALNARRAY{$vcols[1]}},  $vcols[3];  # grow the the array of alignment lengths
  			
                } else {
			#this hit is new, initialize hashes
			
  			$PERCARRAY{$vcols[1]}=[$vcols[2]];         # initialize the array for percent identities
  			$COVARRAY{$vcols[1]}=[$vcols[12]];         # initialize the array for coverage scores
  			$ALNARRAY{$vcols[1]}=[$vcols[3]];          # initialize the array for alignment lengths
  			
                }
		$numreads++;                                       # count reads
		$SEENPERC{$vcols[2]}{$vcols[1]}=1;                 # keep record of percIdent for this report
	} #end while read
	
	close(VSEARCH);
	# done reading all results, lets generate ouptput lines
	# remember, we will only generate up to $maxhits per sample-amplicon pair

	
	# generate the lines for non-empty vsearch reports in this loop
	
	# outer loop by perIdent, to generate only the top n hits where n=maxhits
	# inner loop by hits stored in MAXPERC hash
	
        foreach my $perc ( sort keys %SEENPERC ) {
        	foreach my $hit (  sort  keys %{ $SEENPERC{$perc} } ) {
        		if ( ! defined $WRITTEN{$hit} ) {
        			# we mark hit off in the hash and then use it to calculate which ones are missing
        			$WRITTEN{$hit}=1;
        		
				# we generate a new line in the summary if we have not reached the limit

				last if $numhits > $maxhits;

				# calculate the stats 
				
				my @percIdentity=();      # reset array of percent identity values				
				my @align=();             # reset array of alignment lengths
				my @coverage=();          # reset array of coverage scores

				
				@percIdentity = sort {$a <=> $b} @{$PERCARRAY{$hit}};
				@align        = sort {$a <=> $b} @{$ALNARRAY{$hit}};
				@coverage     = sort {$a <=> $b} @{$COVARRAY{$hit}};

				my $size = @percIdentity;
				my $last=$size-1;
				my $numhits=$size;

				my $max_percIdent = $percIdentity[$last];
				my $max_aln       = $align[$last];
				my $max_cov       = $coverage[$last];
				
				my $min_percIdent = $percIdentity[0];
				my $min_aln       = $align[0];
				my $min_cov       = $coverage[0];
				
				my $sum_percIdent = 0;  map { $sum_percIdent += $_ } @percIdentity;
				my $sum_aln = 0;        map { $sum_aln += $_ }       @align;
				my $sum_cov = 0;        map { $sum_cov += $_ }       @coverage;				
				
				my $avg_percIdent= ( $sum_percIdent / $numhits );
				my $avg_aln      = ( $sum_aln / $numhits );
				my $avg_cov      = ( $sum_cov / $numhits );
				
				# calculate the perentile				

				$percentileIndex=int( ($percentile/100) * $size );
				
				#print "\npercIdent elms= @percIdentity\n";
				#print "\ncoverage elms= @coverage\n";
				#print "\nalignLength elm= @align\n";
				
				#print "\nsum_percIdent=$sum_percIdent\nsum_cov=$sum_cov\nsum_aln=$sum_aln";
				#print "percentile specified: $percentile\tnumhits: $numhits\tindex of percentile: $percentileIndex\tvalue of percentile: $percIdentity[$percentileIndex]\n";
				
				# put together the line and write it to output file

				$detline="$ampliconname\t$samplename\t$hit\t$numreads\t$numhits\t$min_percIdent\t$max_percIdent\t$avg_percIdent\t$min_aln\t$max_aln\t$avg_aln\t$min_cov\t$max_cov\t$avg_cov\t$percIdentity[$percentileIndex]";

				$detline =~ s/\r\n\t/\t/g;   # to delete some pesky DOS characters

				print OUTSUM "$detline\n";

			} # end if witten
		} # end foreach hit SEENPERC

		last if $numhits > $maxhits;        	
        } #end foreach percSeen


########################
# END Processing inputfile
########################


print "\n############  end processing $inputfile\n";
