# Author: Gloria Rendon
# Date: August 2016

# Script to parse the vsearch reports written in blast_output_tabular format
# an create a tab delimited file with summary results.


print "\n\n##############################################################################";
print "\n############   STARTS VSEARCH_CALCSUMMARY_BY_SAMPLE_FORV2.SH           #########";
print "\n##############################################################################\n\n"; 



########################
# global variables
########################

%PRESENT=();
%XREF=();
%AMPLICON=();
%SAMPLE=();
%OTHER=();

$inputfileNames="";
$outputfilePrefix="";
$outputfile="";
$percentile="";
$maxhits="";
$count=0;
$rootdir="";
$sampleDir="";
$entry="";
$type="";
$sampleid="";
$tableDB="database.tsv";
$suffix="vsearch.tsv";
$ampliconfile="BARCODE_names.txt";


print "\n##############################################################################";
print "\n############    reading input from stdin and sanity check\n";
print "\n##############################################################################\n";

#if ($#ARGV<3) {
#	print "script that parses all vsearch reports for amplicons, then generates a summary of results.\n\n";
#	print "Each amplicon has its own folder of results\n";
#	print "EXAMPLE: $0 <outputdir> <sampleid> <outputfilePrefix> <maxhits>\n\n";
#	exit 0;
#}

$sampleDir=shift;          #inputdir where the raw vsearch results are
$sampleid=shift;           #sample name
$outputfile = shift;       #filename of the summary report
$maxhits = shift;          #threshold for maximum number of lines to report per sample-amplicon pair
$percentile = shift;       #only one percentile value will be caculated and reported

die("No outputdir was specified\n")            if $sampleDir  eq "";
die("No sampleid was specified\n")             if $sampleid  eq "";
die("No $outputfile  was specified\n")         if $outputfile eq "";
die("No value for maxhits was specified\n")    if $maxhits eq "";
die("No value for percentile was specified\n") if $percentile eq "";

print "\n##############################################################################";
print "\n############   input arguments have been parsed successfully\n";
print "\n##############################################################################";
print "\nparam1 sampledir=$sampleDir";
print "\nparam2 sampleid=$sampleid";
print "\nparam3 outputfile=$outputfile";
print "\nparam4 maxhits=$maxhits";
print "\nparam5 percentile=$percentile";

open (OUTSUM,">>$outputfile") || die("cannot open outputfile $outputfile"); #we need to open it to append because it already has a header

print "\n############   prepare outputfile= $outputfile\n";


#set +x;
#echo -e "\n\n##############################################################################"      
#echo -e "\n\n############   Main Program starts here                      #################\n\n"
#echo -e "\n\n##############################################################################">&2; set -x;  



########################
# hashing master table of amplicon ids
########################

print "\n############   hashing $ampliconfile\n";


open(AMPLICONFILE,"<$ampliconfile") || die("Cannot open file $ampliconfile \n");
while ( my $line = <AMPLICONFILE> ) {
	$line =~ s/\R//;
	chomp $line;
	$len= length $line;
	if ( $len > 0 ) {
	       $AMPLICON{$line}=1;
	       print "hashing amplicon=$line\n";
	       $count++;
	}
}
close(AMPLICONFILE);

print "\n############   done reading list of amplicons in $ampliconfile\t$count records read\n";

sleep(5);
$count=0;
		

########################
# hashing master table with meta level info for the sequences, it is a tab delimited file with EXACTLY four columns in this order
#
# col0=SEQIDnnnn  where nnnn start at 1, this string is the unique identofier that we use internally for this sequence	
# col1=string with extended sequence identifier arranged like this: taxid_xxxxxx_gi_yyyyyyy_accessionInfo_SpeciesScientificName	if there was a match to a taxid
# col1=with the string manualAnnot	if there was NOT a match to a taxid, we will use the original defline for tax identification purposes    
# col2=list of file names containing this sequence	
# col3=original definition line 
########################

print "\n############   hashing $tableDB\n";

open(TABLEXREF,"<$tableDB") || die("Cannot open file $tableDB \n");
while ( my $line = <TABLEXREF> ) {
	chomp $line;
	if ( $count == 0 ) {
	       # skipping header line
	       $count++; 
	} else { 
	       @det=split(/\t/,$line);
	       $seqidExpanded=$det[0]."_".$det[1]             if $det[1] =~ /taxid/;
	       $seqidExpanded=$det[0]."_manualAnnot_".$det[3] if $det[1] !~ /taxid/;               
	       $XREF{$seqidExpanded}{"DEFLINE"}=$det[3];
	       $XREF{$seqidExpanded}{"FILENAMES"}=$det[2];
	       #print "reading $line\nhashing $seqidExpanded\n";               
	       $count++; 
	}
}
close(TABLEXREF);

print "\n############   done reading $tableDB table of metalevel info about fasta sequences used to perform the vsearches\t$count records read\n";

sleep(10);
$count=0;

########################
# reading directory to process the vsearch reports
########################
#$path="$sampleDir";   ##paths are a pain in perl; we're using a hack inside vsearch_by_sample_forV2.sh -- the script that invokes this one
$path=".";
opendir( DIR, $path ) || die("Can't open $path: $!");

print "\n############   reading sample directory $path\n";

sleep(5);

while ( $entry = readdir( DIR ) ) {
	print "next file in directory: $entry\n";
	next if ($entry =~ /^\.{1,2}$/);

	if ( $entry =~ /_vsearch.tsv/ ) {
		print "\n############   processing vsearch file $entry\n";
		processVsearch($entry,$sampleid);
		print "\n############   done processing vsearch file $entry\n";	
	}

}

closedir( DIR );

print "\n############   DONE reading directory $path. vsearch summary report generated for $sample and written to $outputfile\n";
sleep(5);

########################
# DONE reading directory and processing vsearch reports for this sample
########################

print "\n############   Now we need to add lines with zero hits for the amplicons for which there were no results\n";

processMissing($sampleid);


print "\n\n##############################################################################";
print "\n############   ENDS VSEARCH_CALCSUMMARY_BY_SAMPLE_FORV2.SH           #########";
print "\n##############################################################################\n\n" ;

close(OUTSUM);

exit 0;


########################
# additional functions begin here
########################
sub processVsearch {
	my ($filename,$samplename)=@_;  # make a local copy of the arguments passed to this call of the function
	
	print "\n############   inside processVsearch with these arguments: filename=$filename sample=$samplename\n";
	
	open(VSEARCH,"<$filename") || die("cannot read file $filename\n");
        $ampliconname="";
	
        # parse filename to get amplicon and fill hash PRESENT
       
	foreach my $ampl ( sort keys %AMPLICON ) {
		my $apliconplus="-".$ampl."_vsearch.tsv";
		if  ( $filename =~ /$apliconplus/ ) {
			$ampliconname=$ampl;
			print "amplicon=$ampl\n";
			$PRESENT{$ampl}=1;
			last;
		}
	}


	#skip this process if file is empty
	print OUTSUM "$ampliconname\t$samplename\t0\n" if -z $filename;
	return if -z $filename;

	# file is non-empty, continue with process
	
	######################
	# declare local variables
	######################
	
	my %PERCARRAY=();         # hash for percentId values for a hit
	my %COVARRAY=();          # hash for coverage scores for a hit
	my %ALNARRAY=();          # hash for alignment values for a hit
	my %SEENPERC=();          # percentId values for the entire report
	my %WRITTEN=();
	my $numreads=0;           # number of reads for the entire report 

 	
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
        
        return;
} # end sub


sub processMissing {

	my ($samplename)=@_;
	print "\n############   inside processMissing $samplename ";
	
	# next look to check which amplicons were processed already
	
	foreach my $ampl ( sort keys %AMPLICON ) {
		print "checking if we have results for $ampl\n"; 
		if ( ! defined $PRESENT{$ampl} ) {
			print "no hits for sample=$samplename amplicon=$ampl\n";
			print OUTSUM "$ampl\t$samplename\t0\n";
		}
	
	}
	
	return;
}
