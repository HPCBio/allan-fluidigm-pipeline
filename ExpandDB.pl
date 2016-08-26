#########################################################################
# Author: Gloria Rendon
# Date: July 2016
# script to expand database of fasta sequences with these headers: SEQID_TAXID_ACCESSION_SCCINAME DEFLINE 
#########################################################################

#########################################################################
# global variables
#########################################################################
%FASTA=();
%XTABLE=();
$infastaDB="";
$intableDB="";
$outfastaDB="";
$outtableDB="";
$inBlast="";
$countSeqs=0;
$countRows=0;
$countHits=0;
$prevseq="";
$maxhits=5;
$counthits=0;
@array=();

        
#########################################################################
# Main program starts here
#########################################################################

if ($#ARGV<5) {
	print "script to expand database of fasta sequences with headers like this: SEQID_TAXID_ACCESSION DEFLINE\n\n";
	print "EXAMPLE: perl $0  blast.report indatabase.fa indatabase.tsv  outdatabase.fa outdatabase.tsv redeoDB\n\n";
	exit 0;
}

$inBlast=shift;
$infastaDB=shift;
$intableDB=shift;
$outfastaDB=shift;
$outtableDB=shift;
$redoDB=shift;  # file with sequences for which taxid assignment could not be made with blast results

die("missing parameter blast.report ")    if $inBlast eq "";
die("missing parameter indatabase.fa ")   if $infastaDB eq "";
die("missing parameter outdatabase.fa ")  if $outfastaDB eq "";
die("missing parameter indatabase.tsv ")  if $intableDB eq "";
die("missing parameter outdatabase.tsv ") if $outtableDB eq "";


open(OUTTABLEDB,">$outtableDB")  || die("Cannot open $outtableDB \n");
open(OUTFASTADB,">$outfastaDB")  || die("Cannot open $outfastaDB \n");
open(REDODB,">$redoDB")  || die("Cannot open $redoDB \n");

print "Inside ExpandDN.pl with $inBlast $infastaDB $intableDB $outfastaDB $outtableDB  $redoDB...\n\n";
$countSeqs=ReadInFASTA($infastaDB);

print "done processing $infastaDB\t$countSeqs records processed\n";
$countRows=ReadinTABLE($intableDB);

print "done processing $intableDB\t$countRows records processed\n";
$countHits=ReadInBLAST($inBlast);

print "done processing $inBlast\t$countHits records processed\n";

#reset data structures
%FASTA=();
close(OUTTABLEDB);
close(OUTFASTADB);
close(REDODB);

print "Done processing data. Exiting now\n\n";

exit 0;

#########################################################################
# Functions begin here
#########################################################################

sub ReadInBLAST {
	my ($inBlast)=@_;
	print "inside ReadInBLAST with file=$inBlast\n\n";
	my $redocount=0;
	my $count=0;
	my $newrecs=0;
        my $skipped=0;
        my $hitline="";
        my $seqid="";
        my $rule=0;
	######################################################################
	# print "Start loop to read and process $inBlast\n\n";	
	# table with 17 columns. NO HEADER LINE
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
	# col12=qlen
	# col13=slen
	# col14=sacc
	# col15=staxids
	# col16=sscinames
	######################################################################
	
        open(BLAST,"<$inBlast")  || die("Cannot open $inBlast \n");
        while ( $hitline= <BLAST>) {
		# new blast hit. let's parse it and then apply classification rules        
        	$hitline=~ s/\n//;
        	@detline=split(/\t/,$hitline);
        	$seqid=$detline[0];
        	$seqid =~ s/ //;
        	
        	if ( ! defined $FASTA{$seqid} ) {
			print "no record found for $FASTA{$seqid}\n";
        		$skipped++;

        	} elsif ( $seqid eq $prevseq ) {
        		# same sequence, add the row to the array
        		$counthits++;
        		$array[$counthits]=$hitline;        		
        	} else {
        		# new sequence. two cases:  first batch  or new batch of blast hits
        		if ( $count == 0 ) {
        			# first record, no need to wrap up previous seqid. add the row to the array
				$array[$counthits]=$hitline;        		
        		} else {
        			# NOT the first record, wrap up previous seqid 
        			# first: let's get the new taxid info
        			($rule,$xtraID)=selectTAXID(\@array,$prevseq);
        			$xtraID =~ s/\r\n//g;
        			
        			# second: let's get the other pieces too
        			$prevDefline=$FASTA{$prevseq}{"DEFLINE"};
        			$prevDefline=~ s/\r\n//g;
        			$prevFiles=$FASTA{$prevseq}{"FILES"};
        			$prevSeqBody=$FASTA{$prevseq}{"SEQUENCE"};
        			
        			# third: let's put all the pieces together
				$newHeader=$prevseq."_taxid_".$xtraID if $rule != 5;
				$newHeader=$prevseq."_manualAnnot_".$xtraID       if $rule == 5;
				
				$newtableDet=$prevseq."\ttaxid_".$xtraID."\t".$prevFiles if $rule != 5;
				$newtableDet=$prevseq."\tmanualAnnot_".$xtraID."\t".$prevFiles  if $rule == 5;
								
				#write updated record to table				
				print OUTTABLEDB "$newtableDet\n";
				
				#write updated record to fasta database
				print OUTFASTADB ">$newHeader\n$prevSeqBody" if $prevSeqBody ne "";			
				
				#write redo file with two columns sequence-accession so that it can be done manually 
				
				$accession="";

				if ( $rule == 5 ) {
					# we are left with sequences that did not get taxid assignment
					# lets parse the defline to grab the accession info only
					if    ( $prevDefline =~/\|gb\|(.*)\|/ )  { $accession=$1; }
					elsif ( $prevDefline =~/\|ref\|(.*)\|/ ) { $accession=$1; }
					elsif ( $prevDefline =~/\|dbj\|(.*)\|/ ) { $accession=$1; }				
					elsif ( $prevDefline =~/\|emb\|(.*)\|/ ) { $accession=$1; }
					else  { 
						# nothing has worked. let's try parse the defline and compare each chunk
						@parts=split(/\|/,$prevDefline);
						foreach $chunk (@parts) {
							print "matching to chunk=$chunk\n"; 
							if ( $chunk =~ /[A-Z0-9]+\.[0-9]{1}/ ) { print "matched to $chunk\n"; $accession=$chunk; last; }
							elsif ( $chunk =~ /(NC|NR)_[0-9]+\.[0-9]{1}/ ) { print "matched to $chunk\n"; $accession=$chunk; last; }
							else { print "no match. continue\n"; }
						}
						if ( $accession eq "" ) { 
							# still no luck. we give up parsing this line. you're on your own
							$accession=$prevDefline; 
						}
					}
					print REDODB "$prevseq\t$accession\n";	
					$redocount++;
				}
				
				# reset variables for the next batch of blast results
				$counthits=0;
				@array=();
				$xtraID="";
				$rule=0;
				$newHeader="";
				$newtableDet="";
				$prevSeqBody="";
				$newrecs++;
				
				# add the row to the array
				$array[$counthits]=$hitline;

			}
        	}
		$count++;
        	$prevseq=$seqid;        
        }
        close(BLAST);
        print "$skipped records skipped\n$redocount redo records\n$newrecs total calls to selectTAXID funtion\n";
	return $count;
}

sub selectTAXID {
	my ( $array_ref, $seqid ) = @_;
	my @hits = @{ $array_ref };
	my $row="";
	my $answer="";
	my $order=0;
	my $bestline="";
	my $bestsciname="";
	
	#print "\n\ninside selectTAXID.\n\n";
		
	my $defline=$FASTA{$seqid}{"DEFLINE"};
	$defline =~ s/_/ /g;
	
	#print "\n\nold defline=$seqid  $defline\n";

	######################################################################
        # go over elements of @hits to select the best candidate
	# array has up to five rows. Each row has  17 columns
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
	# col12=qlen
	# col13=slen
	# col14=sacc
	# col15=staxids
	# col16=sscinames
	######################################################################


	foreach $row (@hits) {
		#print "hits\[$order\]=$row\n"; 
		
		#apply classification heuristics 
		@blastout=split(/\t/,$row);
		$species=$blastout[16];
		$species=~ s/ /_/g;
				
		if ( ($blastout[2] == 100 ) && ( $defline  =~ m/$blastout[14]/ ) && ( $defline =~ m/$blastout[16]/ ) ) {
			# rule1: percident=100 && definition line matches accession and species name
			$answer=$blastout[15]."_".$blastout[1]."_".$species; #TAXID_SUBJECTID_SSCINAME
			#print "rule1=$answer\t$defline\n";
			return (1,$answer);
		} elsif ( ( $defline  =~ m/$blastout[14]/ ) && ( $defline =~ m/$blastout[16]/ )  ) {
			# rule2: definition line matches accession and species name
			$answer=$blastout[15]."_".$blastout[1]."_".$species; #TAXID_SUBJECTID_SSCINAME
			#print "rule2=$answer\t$defline\n";
			return (2,$answer);			
		} elsif ( $defline =~ m/$blastout[14]/ ) {
			# rule3: accession matches definition line	
			$answer=$blastout[15]."_".$blastout[1]."_".$species; #TAXID_SUBJECTID_SSCINAME
			#print "rule3=$answer\t$defline\n";			
			return (3,$answer);	
		#} elsif ( $defline =~ m/$blastout[16]/ ) {
		#	# rule4: it matches species name
		#	$bestline=$defline;
		#	$bestsciname=$blastout[15]."_".$blastout[1]."_".$species; #TAXID_SUBJECTID_SSCINAME
		} elsif ( $order == 0 ) {
			# rule5: anything else, return the manual annotation that appears in the defline
			$bestline=$defline;
			$bestline=~ s/ /_/g;
			$bestline=~ s/\r\n//g;
			$besthit=$blastout[15]."_".$blastout[1]."_".$species; #TAXID_SUBJECTID_SSCINAME
		}
		$order++;
	} # end foreach

        #if ( $bestsciname ne "" ) {
	#	#print "rule4=$besthit\t$bestline\n";	
	#	return (4,$bestsciname); #rule 4 it matches species name only
	#}
	
	#if we reached this point is because there was no match; we need to trust the original annotation
	#print "rule5=$besthit\t$bestline\n";
	return (5,$bestline);                           #rule 5 nothing else matched
	
}


sub ReadinTABLE {
	my ($intableDB)=@_;
	print "inside ReadinTABLE with file=$intableDB\n\n";
	my $count=0;

	######################################################################
	# print "Start loop to read and process $intableDB\n\n";	
	# table header line
	# SEQ_ORDER	FILENAME(S)	DEFINITION_LINE
	######################################################################

	open(INTABLEDB,"<$intableDB")  || die("Cannot open $intableDB \n");
	while (<INTABLEDB>) {
		if ( $count == 0 ) {
		    #skip header row
		    $count++
		} else {
			$row=$_;
			chomp $row;
			( $seqid, $files, $defline )=split(/\t/,$row);
			if ( length($files) > 0 ) { 
				#print "row \[$count\]=$row";
				#print "parsed row \[$count\]:$seqid\t$files\$defline";

				# sanity check. db.fa and db.tsv must match
				$defline=~ s/\n//g;

				die("mismatch in seqids")   if  ! defined $FASTA{$seqid};
				die("missing defline")      if  ! defined $FASTA{$seqid}{"DEFLINE"};		
				#die("mismatch in deflines") if  $FASTA{$seqid}{"DEFLINE"}  ne $defline;
				$FASTA{$seqid}{"FILES"}=$files;
			}
			$count++;
		}
	
	
	} #end while read

	######################################################################
	#print "End loop to read and process $infastaDB\n\n";	
	######################################################################	
	
        close(INTABLEDB);
	return $count;
} #end sub

sub ReadInFASTA {
	my ($infastaDB)=@_;
	print "inside ReadInFASTA with file=$infastaDB\n\n";
	my $count=0;
	local $/ = ">";                      # this is the record delimiter
	open(INFA,"<$infastaDB")  || die("Cannot open $infastaDB \n");

	######################################################################
	# print "Start loop to read and process $infastaDB\n\n";	
	######################################################################
	while (<INFA>) {
		my $seq_all = $_;
		$seq_all =~ s/>//;           # to correct the parsing result by removing > from the end of the record
		$seq_all = ">".$seq_all;     # to correct the parsing result by prepending > to the record
		my $header="";
		my $seq="";
		if ( $seq_all =~ m/>(.*)\n((.|\n)*)/ ) {
			$header=$1;
			$seq=$2;
			chop $header;
			if ( length($seq) > 0 ) {
				($seqid,$defline)=split(/ /,$header);
				$defline=~ s/\r\n//g;
				# let's populate the FASTA hash
				if ( ! defined $FASTA{$seqid} ) {
					#print "record does not exists: $seqid\n";
					$FASTA{$seqid}{"SEQUENCE"}=$seq;
					$FASTA{$seqid}{"DEFLINE"}=$defline;				
				} else {
					print "duplicate sequence. Skipping it\n";
				}
			} else {
				print "empty sequence";
			}
			
		} else {
			print "failed to parse $seq_all in file $inputfile\n";		
		}
		$count++;
		#end processing this sequences, next please
	}  # end while read inputfile

	######################################################################
	#print "End loop to read and process $infastaDB\n\n";	
	######################################################################	
	
	close(INFA);
	return $count;	
}


		
