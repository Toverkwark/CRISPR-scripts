use Getopt::Long;
require "Damerau.pl";
sub MatchBarcode($@);
sub ScoreTwoStrings($$);

print "Usage:perl $0 -input -output -report\n-input\tName of input file\n-output\tName of output file. Default is inputfile.stripped\n-report\tName of report file. Default is inputfile.report\n";

#Define screen specific settings
my $BarcodeOffset = 0; #Position of start of barcode
my $BarcodeLength = 6; #Number of nucleotides that the barcode is long
my $ExpectedLeadingSequence = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
my $ExpectedTrailingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC"; #Sequence that is expected to come after the gRNA/shRNA sequence
my $ErrorThresholdForAnalysis = 3; #This number of mutations of indels can be present in both leading and trialing sequences and it the sequence will still be taken along in the analysis
my $ErrorThresholdForFiltering = 1;#This number of mutations of indels can be present in both leading and trialing sequences and it the sequence will still be accepted in the output filtering

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

@Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
%Results  = ();
open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
open( OUTPUT, ">", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";
open( REPORT, ">", $ReportFile ) or die "ERROR in $0:Report file $ReportFile is not accessible.\n";

while ( defined( my $line = <INPUT> ) ) {
	my $BarcodeMatched;
	my $PromoterFound;
	my $PromoterExactlyMatched;
	my $TracrFound;
	my $TracrFoundAtCorrectPosition;
	
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 100000 ) ) {
		print "Analyzing record $RecordsAnalyzed\n";
	}
	$line2    = <INPUT>;
	$line3    = <INPUT>;
	$line4    = <INPUT>;
	chomp($line2);
	$Sequence = $line2;

	#Get the barcode. See if it exists. If not, try to map it with maximally 1 nucleotide replacement and only 1 match existing.
	#if that can be found, make sure to also write the output file sequences with the mapped barcode
	$Barcode = substr( $Sequence, $BarcodeOffset, $BarcodeLength );
	if ( grep( /$Barcode/, @Barcodes ) ) {
		$BarcodeMatched=1;
		$Results{$Barcode}->[0]++;
	}
	else {
		my $MatchedBarcode = MatchBarcode( $Barcode, @Barcodes );
		if ($MatchedBarcode) {
			$Barcode = $MatchedBarcode;
			$BarcodeMatched=1;
			$Results{$Barcode}->{'BarcodesFound'}++;
			$Results{$Barcode}->{'BarcodesFoundByMatching'}++;
			#Change the barcode in the outputfile to the one it was matched to
			$line2 = substr($Sequence,0,$BarcodeOffset) . $Barcode . substr($Sequence,$BarcodeLength+$BarcodeOffset,length($Sequence-$BarcodeOffset)-$BarcodeLength);
		}
	}

	#Try to find the leading sequence and record any mistakes there. Determine an offset in case it is found
	my $Offset = 0;
	if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset), length($ExpectedLeadingSequence)) eq $ExpectedLeadingSequence ) {
		$Results{$Barcode}->{'LeadingSequencesFound'}++;
		$Results{$Barcode}->{'ExactLeadingSequencesFound'}++;
	}
	else {
		my %DamerauResults;
		my $NotConvergedOffset=1;
		while ($NotConvergedOffset) {
			DetermineDamerauLevenshteinDistance($ExpectedLeadingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset), length($ExpectedLeadingSequence),\%DamerauResults);
			if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdForAnalysis) {
				#Record the errors found in the Leading sequence
				#Determine the offset for the inserts
				#Output the results only if threshold is matched
				#Determine if offset is converged
			}
		}		
	}
	
	#Extract all construct features assuming all offsets are correct
	$promoter = substr( $Sequence, 6 + $Offset,  24);
	$gRNA  = substr( $Sequence, 30 + $Offset, 20 );
	$tracr = substr( $Sequence, 50 + $Offset);

	#Check if the promoter is intact
	if($promoter eq 'CCCTATCAGTGATAGAGACTCGAG') {
		$Results{$Barcode}->[3]++;
		$PromoterExactlyMatched=1;
	}

	#Make histogram of where tracrs are found, which will indicate amount of 20nt inserts
	#Only proceed with those reads where the tracr is in the correct position
	for(my $i=1;$i<=length($Sequence);$i++) {
		if(substr($Sequence,$i,8) eq 'GTTTTAGA') {
			$TracrFound=1;
			$Results{$Barcode}->[4]++;
			$Results{$Barcode}->[5]++ if ($PromoterExactlyMatched);
			$Histo{$i-50-$Offset}++;
			if($i==50+$Offset) {
				$TracrFoundAtCorrectPosition=1;
				$Results{$Barcode}->[6]++;
			}
		}	
	}

	#If all criteria are met, count this read as a correct read
	if ( $BarcodeMatched && $PromoterFound && $PromoterExactlyMatched && $TracrFound && $TracrFoundAtCorrectPosition) {
		$Results{$Barcode}->[7]++;
		print OUTPUT $Barcode . $promoter . $gRNA . $tracr . "\n";
	}
}


#Output the tracrRNA positional histogram to the report file
print REPORT "tracrRNAs found relative to expected:\n";
foreach $HistoDistance (sort {$a <=> $b} keys %Histo) {
	print REPORT "Distance\t" . ($HistoDistance) . "\t" . $Histo{$HistoDistance} . "\n";
}
print REPORT "\n\n";


print REPORT "Barcode\tReads\tOf which by match\tPromoter seeds found\tCompletely Correct Promoters\tTracr found\tTracr found+Correct Promoter\tTracr found 20nt after promoter\tCompletely correct\n";


foreach $Barcode ( @Barcodes ) {
	print REPORT "$Barcode\t";
	for ( $counter = 0 ; $counter <= 7 ; $counter++ ) {
		print REPORT sprintf( "%i", $Results{$Barcode}->[$counter] ) . "\t";
		$Totals[$counter] += $Results{$Barcode}->[$counter];
	}
	if ( $Results{$Barcode}->[0] > 0 ) {
		print REPORT sprintf( "%4.2f%%", 100 * ( $Results{$Barcode}->[7] / $Results{$Barcode}->[0] ) ) . "\n";
	}
	else {
		print REPORT sprintf( "%4.2f%%", 0 ) . "\n";
	}
}
print REPORT sprintf( "TOTAL\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", $Totals[0], $Totals[1], $Totals[2], $Totals[3], $Totals[4], $Totals[5], $Totals[6], $Totals[7]);

print REPORT sprintf(
	"\n\t\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%",
	100 * $Totals[1] / $Totals[0],
	100 * $Totals[2] / $Totals[0],
	100 * $Totals[3] / $Totals[0],
	100 * $Totals[4] / $Totals[0],
	100 * $Totals[5] / $Totals[0],
	100 * $Totals[6] / $Totals[0],
	100 * $Totals[7] / $Totals[0],
);
print REPORT "\n\nAnalyzed reads\t" . sprintf( "%i", $RecordsAnalyzed );
close(INPUT)  or die "Could not close input file $InputFile.\n";
close(OUTPUT) or die "Could not close output file $OutputFile.\n";
close(REPORT) or die "Could not close report file $ReportFile.\n";

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
