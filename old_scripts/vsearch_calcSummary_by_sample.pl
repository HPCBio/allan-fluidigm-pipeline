# Author: Gloria Rendon
# Date: March 2016

# Script to parse the vsearch reports written in blast_output_tabular format
# an create a tab delimited file with summary results.

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
$maxhits="";
$count=0;
$rootdir="";
$sampleDir="";
$entry="";
$type="";
$sampleid="";
$tableDB="database.tsv";
$suffix="vsearch.tsv";
$ampliconfile="AMPLICON_names.txt";


########################
# reading input from stdin and sanity check
########################

if ($#ARGV<4) {
	print "script that parses all vsearch reports for amplicons, then generates a summary of results.\n\n";
	print "Each amplicon has its own folder of results\n";
	print "EXAMPLE: $0 <outputdir> <sampleid> <inputfileNames> <outputfilePrefix> <maxhits>\n\n";
	exit 0;
}

$sampleDir=shift;
$sampleid=shift;
$inputfileNames = shift;   #file with list of vsearch reports for a sample
$outputfile = shift;       #filename of the summary report
$maxhits = shift;          #threshold for maximum number of lines to report per sample-amplicon pair

die("No outputdir was specified\n")            if $sampleDir  eq "";
die("No sampleid was specified\n")             if $sampleid  eq "";
die("No inputfileNames was specified\n")       if $inputfileNames  eq "";
die("No $outputfile  was specified\n")         if $outputfile eq "";
die("No value for maxhits was specified\n")    if $maxhits eq "";

print "input arguments have been parsed successfully\n";

# files to be processed

$fileNames=`cat $inputfileNames`;
$fileNames=~ s/:/\n/g;
print "files to process: $fileNames\n\n";

print "open outputfile= $outputfile\n";
open (OUTSUM,">>$outputfile") || die("cannot open outputfile $outputfile");

#set +x;
#echo -e "\n\n##############################################################################"      
#echo -e "\n\n############   DO NOT EDIT THE SCRIPT BEYOND THIS POINT      #################\n\n"
#echo -e "\n\n##############################################################################">&2; set -x; 

#set +x;
#echo -e "\n\n##############################################################################"      
#echo -e "\n\n############   Main Program starts here                      #################\n\n"
#echo -e "\n\n##############################################################################">&2; set -x;  



########################
# hashing master table of amplicon ids
########################

print "hashing $ampliconfile\n";
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
print "done reading list of amplicons in $ampliconfile\t$count records read\n";
sleep(5);
$count=0;
		

########################
# hashing master table with meta level info for the sequences, it is a tab delimited file with EAXCTLY four columns in this order
#
# col0=SEQIDnnnn  where nnnn start at 1, this string is the unique identofier that we use internally for this sequence	
# col1=string with extended sequence identifier arranged like this: taxid_xxxxxx_gi_yyyyyyy_accessionInfo_SpeciesScientificName	if there was a match to a taxid
# col1=with the string manualAnnot	if there was NOT a match to a taxid, we will use the original defline for tax identification purposes    
# col2=list of file names containing this sequence	
# col3=original definition line 
########################

print "hashing $tableDB\n";
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
print "done reading $tableDB table of metalevel info about fasta sequences used to perform the vsearches\t$count records read\n";
sleep(10);
$count=0;

########################
# reading directory to process the vsearch reports
########################
#$path="$sampleDir";
$path=".";
opendir( DIR, $path ) || die("Can't open $path: $!");
print "reading sample directory $path\n";
sleep(5);

while ( $entry = readdir( DIR ) ) {
	print "next file in directory: $entry\n";
	next if ($entry =~ /^\.{1,2}$/);

	if ( $entry =~ /_vsearch.tsv/ ) {
		print "processing vsearch file $entry";
		processVsearch($entry,$sampleid);
		print "done processing vsearch file $entry";	
	}

}

closedir( DIR );
print "DONE reading directory $path. vsearch summary report generated for $sample and written to $outputfile\n";
sleep(5);

########################
# DONE reading directory and processing vsearch reports for this sample
########################

print "\n\nNow we need to add lines with zero hits for the amplicons for which there were no results\n";

processMissing($sampleid);

print "\n\nDone. Exiting now\n\n\n";

exit 0;


########################
# additional functions begin here
########################
sub processVsearch {
	my ($filename,$samplename)=@_;  # make a local copy of the arguments passed to this call of the function
	print "\n\ninside processVsearch with these arguments: filename=$filename sample=$samplename\n";
	open(VSEARCH,"<$filename") || die("cannot read file $filename\n");
        $ampliconname="";
	
        # parse filename to get amplicon and fill hash PRESENT
        
        if  ( $filename =~ /(.*)_$samplename/ ) {  
		$ampliconname=$1;
		print "amplicon=$ampliconname\n";
		$PRESENT{$ampliconname}=1;
	}

	#skip this process if file is empty
	print OUTSUM "$samplename\t$ampliconname\t0\n" if -z $filename;
	return if -z $filename;

	# file is non-empty, continue with process
	
	######################
	# declare local variables
	######################
	
	my %NUMHITS=();           # numhits for the same subject id
	my %PERCSUM=();           # percentId values for the same subject id
	my %PERCMAX=();
	my %PERCMIN=();
	my %COVMAX=();
	my %COVMIN=();
	my %ALNMAX=();
	my %ALNMIN=();
	my %COVSUM=();            # coverage values for the same subject id
	my %ALNSUM=();            # alnLength values for the same subject id
	my %SEENPERC=();          # percentId values for the entire report
	my %WRITTEN=();
	my $numreads=0;           # number of reads for the entire report 
	my $numhits=0;            # keep track of hits reported
	
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

		# here we simply populate the hashes, but we do not write summaries yet
                if ( defined $PERCSUM{$vcols[1]} ) {
			#this subject already exists, update values of hashes w stats
			
			$PERCMAX{$vcols[1]} = $vcols[2] if $vcols[2] > $PERCMAX{$vcols[1]};
			$PERCMIN{$vcols[1]} = $vcols[2] if $vcols[2] < $PERCMIN{$vcols[1]};
                        $PERCSUM{$vcols[1]} += $vcols[2];

                        $ALNMAX{$vcols[1]}  = $vcols[3] if $vcols[3] > $ALNMAX{$vcols[1]};                        
                        $ALNMIN{$vcols[1]}  = $vcols[3] if $vcols[3] < $ALNMIN{$vcols[1]};
                        $ALNSUM{$vcols[1]}  += $vcols[3];

			$COVMAX{$vcols[1]}  = $vcols[12] if $vcols[12] > $COVMAX{$vcols[1]};                        
			$COVMIN{$vcols[1]}  = $vcols[12] if $vcols[12] < $COVMIN{$vcols[1]}; 
                        $COVSUM{$vcols[1]}  += $vcols[12]; 
                                                
  			$NUMHITS{$vcols[1]}++;                      
                } else {
			#this subject is new, initialize hashes
			
                        $PERCSUM{$vcols[1]}=$vcols[2];         
                        $PERCMIN{$vcols[1]}=$vcols[2];          
                        $PERCMAX{$vcols[1]}=$vcols[2];          

                        $COVSUM{$vcols[1]} =$vcols[12]; 
                        $COVMIN{$vcols[1]} =$vcols[12]; 
                        $COVMAX{$vcols[1]} =$vcols[12]; 

                        $ALNSUM{$vcols[1]} =$vcols[3];
                        $ALNMIN{$vcols[1]} =$vcols[3];                        
                        $ALNMAX{$vcols[1]} =$vcols[3]; 
                        
  			$NUMHITS{$vcols[1]} =1;                                   
                }
		$numreads++;                                     #count reads
		$SEENPERC{$vcols[2]}{$vcols[1]}=1;               #keep record of percIdent for this report
	} #end while read
	
	close(VSEARCH);
	# done reading all results, lets generate ouptut lines
	# remember, we will only generate up to $maxhits per sample-amplicon pair

	
	# generate the lines for non-empty vsearch reports in this loop
	
	# outer loop by perIdent, to generate only the top n hits where n=maxhits
	# inner loop by hits stored in MAXPERC hash
	
        foreach my $perc ( sort keys %SEENPERC ) {
        	foreach my $hit (  sort  keys %{ $SEENPERC{$perc} } ) {
        		if ( ! defined $WRITTEN{$hit} ) {
        			# we mark it off
        			$WRITTEN{$hit}=1;
        		
				# we generate a new line in the summary if we have not reached the limit
				$numhits++;
				last if $numhits > $maxhits;

				# calculate the stats 
				
				my $min_percIdent = $PERCMIN{$hit};
				my $min_cov       = $COVMIN{$hit};
				my $min_aln       = $ALNMIN{$hit};

				my $max_percIdent = $PERCMAX{$hit};
				my $max_cov       = $COVMAX{$hit};
				my $max_aln       = $ALNMAX{$hit};
				
				my $mean_percIdent= ( $PERCSUM{$hit} / $NUMHITS{$hit} );
				my $mean_cov      = ( $COVSUM{$hit} / $NUMHITS{$hit} );
				my $mean_aln      = ( $ALNSUM{$hit} / $NUMHITS{$hit} );

				# put together the line and write it to output file

				$detline="$samplename\t$ampliconname\t$hit\t$numreads\t$NUMHITS{$hit}\t$min_percIdent\t$max_percIdent\t$mean_percIdent\t$min_aln\t$max_aln\t$mean_aln\t$min_cov\t$max_cov\t$mean_cov";
				$detline =~ s/\n\t/\t/g;
				print OUTSUM "$detline\n";
			} # end if witten
		} # end foreach hit SEENPEC

		last if $numhits > $maxhits;        	
        } #end foreach percSeen
        
        return;
} # end sub


sub processMissing {

	my ($samplename)=@_;
	print "inside processMissing $samplename ";
	
	# next look to check which amplicons were processed already
	
	foreach my $ampl ( sort keys %AMPLICON ) {
		print "checking if we have results for $ampl\n"; 
		if ( ! defined $PRESENT{$ampl} ) {
			print "no hits for sample=$samplename amplicon=$ampl\n";
			print OUTSUM "$samplename\t$ampl\t0\n";
		}
	
	}
	
	return;
}