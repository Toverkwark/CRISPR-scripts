use Getopt::Long;
use Term::ANSIColor;
use strict;
require "Damerau.pl";
sub MatchBarcode($@);
sub ScoreTwoStrings($$);

print "Usage:perl $0 -input -output -report\n-input\tName of input file\n-output\tName of output file. Default is inputfile.stripped\n-report\tName of report file. Default is inputfile.report\n";

my $StartTime=time;

#Define screen specific settings
my $BarcodeOffset = 0; #Position of start of barcode
my $BarcodeLength = 6; #Number of nucleotides that the barcode is long
my $ExpectedInsertLength = 20; #Number of nucleotides of the insert between leading and trailing sequence
my $ExpectedLeadingSequence = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
my $ExpectedTrailingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC"; #Sequence that is expected to come after the gRNA/shRNA sequence
my $ErrorThresholdLeading = 6; #This number of mutations or indels can be present in the leading  sequences
my $ErrorThresholdTrailing = 6; #This number of mutations or indels can be present in the trailing sequences
my @Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
my %Results  = ();
my %LeadingErrors;
my %TrailingErrors;
my %InsertLengths;
my $InputFile;
my $OutputFile;
my $ReportFile;
my $RecordsAnalyzed;

GetOptions(
	"input=s"  => \$InputFile,
	"output=s" => \$OutputFile,
	"report=s" => \$ReportFile,
);

if ( !$OutputFile ) {
	$OutputFile = $InputFile . ".stripped";
}

if ( !$ReportFile ) {
	$ReportFile = $InputFile . ".report";
}

open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
open( OUTPUT, ">", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";
open( REPORT, ">", $ReportFile ) or die "ERROR in $0:Report file $ReportFile is not accessible.\n";
open( NOTANALYZED, ">", $InputFile . ".notanalyzed") or die "ERROR in $0:File " . $InputFile . ".analyzed is not accessible.\n";

while ( defined( my $line = <INPUT> ) ) {
	my $BarcodeMatched;
	my $PromoterFound;
	my $PromoterExactlyMatched;
	my $TracrFound;
	my $TracrFoundAtCorrectPosition;
	
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 10000 ) ) {
		print "Analyzing record $RecordsAnalyzed\n";
	}
	my $line2    = <INPUT>;
	my $line3    = <INPUT>;
	my $line4    = <INPUT>;
	chomp($line2);
	my $Sequence = $line2;
	#Get the barcode. See if it exists. If not, try to map it with maximally 1 nucleotide replacement and only 1 match existing.
	#if that can be found, make sure to also write the output file sequences with the mapped barcode
	my $Barcode = substr( $Sequence, $BarcodeOffset, $BarcodeLength );
	if ( grep( /$Barcode/, @Barcodes ) ) {
		$BarcodeMatched=1;
		$Results{$Barcode}->[0]++; #Barcodes found
		$Results{$Barcode}->[1]++; #Exact barcodes found
	}
	else {
		my $MatchedBarcode = MatchBarcode( $Barcode, @Barcodes );
		if ($MatchedBarcode) {
			$Barcode = $MatchedBarcode;
			$BarcodeMatched=1;
			$Results{$Barcode}->[0]++;
			#Change the barcode in the outputfile to the one it was matched to
			$line2 = substr($Sequence,0,$BarcodeOffset) . $Barcode . substr($Sequence,$BarcodeLength+$BarcodeOffset,length($Sequence-$BarcodeOffset)-$BarcodeLength);
		}
	}

	#Try to find the leading sequence and record any mistakes there. Determine an offset in case it is found
	my $LeadingOffset = 0;
	my $LeadingSequenceFound=0;
	my %DamerauResults;
	my $NotConvergedOffset=1;
	if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset), length($ExpectedLeadingSequence)) eq $ExpectedLeadingSequence ) {
		$Results{$Barcode}->[2]++; #Leading sequences found
		$Results{$Barcode}->[3]++; #Exact leading sequences found
		#print "Exact leading sequence found\n";
		$LeadingSequenceFound=1;
	}
	else {
		$NotConvergedOffset=1;
		while ($NotConvergedOffset) {
			undef %DamerauResults;
			DetermineDamerauLevenshteinDistance($ExpectedLeadingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset)),\%DamerauResults);
			if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdLeading) {
				$NotConvergedOffset=0;
				#Record the errors found in the Leading sequence
				foreach my $Change (keys $DamerauResults{'Changes'}) {
					if($Change==(length($ExpectedLeadingSequence)+$LeadingOffset)) {
						if($DamerauResults{'Changes'}->{$Change} eq 'Insertion') {
							$LeadingOffset++;
							$NotConvergedOffset=1;
						}
						if($DamerauResults{'Changes'}->{$Change} eq 'Deletion') {
							$LeadingOffset--;
							$NotConvergedOffset=1;
						}
					}
				}
			}
			else {
				$NotConvergedOffset=0;		
			}
		}		
	}
	
	if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdLeading) {
		$LeadingSequenceFound=1;
		$Results{$Barcode}->[2]++; #Leading sequences found
		#print "Leading sequence found with changes:";
		foreach my $Change (keys $DamerauResults{'Changes'}) {
			$LeadingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
			#print $DamerauResults{'Changes'}->{$Change} . " @" . $Change . " - ";
		}
		#print "\n";
	} 
	#print "No leading sequence found\n" if !$LeadingSequenceFound;
	
	#Try to find the trailing sequence and record any mistakes there. Determine an offset in case it is found
	my $InsertLength=$ExpectedInsertLength;
	my $TrailingSequenceFound=0;
	$NotConvergedOffset=1;
	while ($NotConvergedOffset) {
		undef %DamerauResults;
		#print "Comparing:\n1:$ExpectedTrailingSequence\n2:" . substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)) . "\n";
		if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)) eq $ExpectedTrailingSequence ) {
			$Results{$Barcode}->[4]++; #Trailing sequences found
			$Results{$Barcode}->[5]++; #Exact trailing sequences found
			#print "Exact trailing sequence found\n";
			$TrailingSequenceFound=1;
			$NotConvergedOffset=0;	
		}
		else {
			DetermineDamerauLevenshteinDistance($ExpectedTrailingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)),\%DamerauResults);
			if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdTrailing) {
				$NotConvergedOffset=0;
				#Record the errors found in the Trailing sequence
				foreach my $Change (keys $DamerauResults{'Changes'}) {
					if($Change==0 && $DamerauResults{'Changes'}->{$Change} eq 'Insertion') {
						$InsertLength++;
						$NotConvergedOffset=1;
					}
					if($Change==1 && $DamerauResults{'Changes'}->{$Change} eq 'Deletion') {
						$InsertLength--;
						$NotConvergedOffset=1;
					}
				}
			}
			else {
				$NotConvergedOffset=0;		
			}
		}		
	}
	
	if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdTrailing) {
		$TrailingSequenceFound=1;
		$Results{$Barcode}->[4]++; #Trailing sequences found
		#print "Trailing sequence found with changes:";
		foreach my $Change (keys $DamerauResults{'Changes'}) {
			$TrailingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
			#print $DamerauResults{'Changes'}->{$Change} . " @" . $Change . " - ";
		}
		#print "\n";
	} 
	#print "No trailing sequence found. LeadingOffset:$LeadingOffset ExpectedInsertLength:$ExpectedInsertLength \n" if !$TrailingSequenceFound;
	#print "Insert length:$InsertLength\n";
	
	#Only continue if in any way a leading and trailing sequence are found
	if($LeadingSequenceFound && $TrailingSequenceFound) {
		$InsertLengths{$InsertLength}++;
		if($InsertLength==$ExpectedInsertLength) {
			my $LeadingSequence=substr( $Sequence, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset));
			my $InsertSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)),$InsertLength);
			my $TrailingSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)+$InsertLength));
			print color("red") . $LeadingSequence;
			print color("green") . $InsertSequence;
			print color("blue") . $TrailingSequence;
			print color("reset") . "\n";
			print OUT "$Barcode\t$InsertSequence";
		}
	}
	else {
		print NOTANALYZED "$Sequence\n";
	}
}


#Output the tracrRNA positional histogram to the report file
print REPORT "Insert sizes:\n";
foreach my $InsertLength (sort {$a <=> $b} keys %InsertLengths) {
	print REPORT "Size\t" . ($InsertLength) . "\t" . $InsertLengths{$InsertLength} . "\n";
}
print REPORT "\n\n";


#print REPORT "Barcode\tReads\tOf which by match\tPromoter seeds found\tCompletely Correct Promoters\tTracr found\tTracr found+Correct Promoter\tTracr found 20nt after promoter\tCompletely correct\n";
#
#
#foreach $Barcode ( @Barcodes ) {
#	print REPORT "$Barcode\t";
#	for ( $counter = 0 ; $counter <= 7 ; $counter++ ) {
#		print REPORT sprintf( "%i", $Results{$Barcode}->[$counter] ) . "\t";
#		$Totals[$counter] += $Results{$Barcode}->[$counter];
#	}
#	if ( $Results{$Barcode}->[0] > 0 ) {
#		print REPORT sprintf( "%4.2f%%", 100 * ( $Results{$Barcode}->[7] / $Results{$Barcode}->[0] ) ) . "\n";
#	}
#	else {
#		print REPORT sprintf( "%4.2f%%", 0 ) . "\n";
#	}
#}
#print REPORT sprintf( "TOTAL\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", $Totals[0], $Totals[1], $Totals[2], $Totals[3], $Totals[4], $Totals[5], $Totals[6], $Totals[7]);
#
#print REPORT sprintf(
#	"\n\t\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%",
#	100 * $Totals[1] / $Totals[0],
#	100 * $Totals[2] / $Totals[0],
#	100 * $Totals[3] / $Totals[0],
#	100 * $Totals[4] / $Totals[0],
#	100 * $Totals[5] / $Totals[0],
#	100 * $Totals[6] / $Totals[0],
#	100 * $Totals[7] / $Totals[0],
#);
print REPORT "\n\nAnalyzed reads\t" . sprintf( "%i", $RecordsAnalyzed );
close(INPUT)  or die "Could not close input file $InputFile.\n";
close(OUTPUT) or die "Could not close output file $OutputFile.\n";
close(REPORT) or die "Could not close report file $ReportFile.\n";
close(NOTANALYZED) or die "Could not close file $InputFile.notanalyzed\n";

my $EndTime=time;
print "Time spent analyzing:" . ($EndTime-$StartTime) . " seconds, or " . ($RecordsAnalyzed/($EndTime-$StartTime)) . " sequences per second\n";

sub MatchBarcode($@) {
	my ( $Barcode, @BarcodeList ) = @_;
	my $MatchesFound = 0;
	my $MatchedBarcode;
	foreach my $BarcodeFromList (@BarcodeList) {
		if ( ScoreTwoStrings( $Barcode, $BarcodeFromList ) <= 1 ) {
			$MatchedBarcode = $BarcodeFromList;
			$MatchesFound++;
		}
	}
	if ( $MatchesFound == 1 ) {
		return $MatchedBarcode;
	}
	else {
		return 0;
	}
}

sub ScoreTwoStrings($$) {
	my ( $Barcode, $BarcodeFromList ) = @_;
	my $Deviations = 0;
	for ( my $i = 0 ; $i < ( length $Barcode ) ; $i++ ) {
		if ( ( substr $Barcode, $i, 1 ) ne ( substr $BarcodeFromList, $i, 1 ) ) {
			$Deviations++;
		}
	}
	return $Deviations;
}
