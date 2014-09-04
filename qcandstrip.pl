use Getopt::Long;
use Term::ANSIColor;
use strict;
require "Damerau.pl";
sub MatchBarcode($@);
sub ScoreTwoStrings($$);

print "Usage:perl $0 -input -output -report -library\n-input\tName of input file\n-output\tName of output file. Default is inputfile.stripped\n-report\tName of report file. Default is inputfile.report\n-library\tName of library file to which inserts are mapped\n";

my $StartTime=time;

#Define screen specific settings
my $BarcodeOffset = 0; #Position of start of barcode
my $BarcodeLength = 6; #Number of nucleotides that the barcode is long
my $ExpectedInsertLength = 20; #Number of nucleotides of the insert between leading and trailing sequence
my $ExpectedLeadingSequence = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
my $ExpectedTrailingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC"; #Sequence that is expected to come after the gRNA/shRNA sequence
my $ErrorThresholdLeading = 20; #This number of mutations or indels can be present in the leading  sequences
my $ErrorThresholdTrailing = 20; #This number of mutations or indels can be present in the trailing sequences
my @Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
my %Results  = ();
my %LeadingErrors;
my %TrailingErrors;
my %InsertLengths;
my %QualitiesByBarcode;
my %InsertCounts;
my $InputFile;
my $OutputFile;
my $ReportFile;
my $LibraryFile;
my $RecordsAnalyzed;

GetOptions(
	"input=s"  => \$InputFile,
	"output=s" => \$OutputFile,
	"report=s" => \$ReportFile,
	"library=s" => \$LibraryFile
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
open( LIBRARY, $LibraryFile ) or die "ERROR in $0:Library file $LibraryFile is not accessible.\n";

#Start by reading in the library file
#The format of this file should be [ID] tab [SEQUENCE]
print "Reading library file\n";
my %Library;
my $InsertsFound;
while ( defined( my $Line = <LIBRARY> ) ) {
	$InsertsFound++;
	chomp($Line);
	my @values = split( /\t/, $Line );
	$Library{$values[1]} = $values[0];
}
close(LIBRARY) or die "Could not close file $LibraryFile\n";
print "$InsertsFound inserts found in library file $LibraryFile\n";

while ( defined( my $line = <INPUT> ) ) {
	my $BarcodeFound;
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 100 ) ) {
		print "Analyzing record $RecordsAnalyzed\n";
	}
	my $line2    = <INPUT>;
	my $line3    = <INPUT>;
	my $line4    = <INPUT>;
	chomp($line2);
	my $Sequence = $line2;
	my $Qualities = $line4;
	#Get the barcode. See if it exists. If not, try to map it with maximally 1 nucleotide replacement and only 1 match existing.
	my $Barcode = substr( $Sequence, $BarcodeOffset, $BarcodeLength );
	if ( grep( /$Barcode/, @Barcodes ) ) {
		$BarcodeFound=1;
		$Results{$Barcode}->[0]++; #Barcodes found
		$Results{$Barcode}->[1]++; #Exact barcodes found
	}
	else {
		my $MatchedBarcode = MatchBarcode( $Barcode, @Barcodes );
		if ($MatchedBarcode) {
			$Barcode = $MatchedBarcode;
			$BarcodeFound=1;
			$Results{$Barcode}->[0]++;
		}
	}

	#Only continue to analyze this read if a valid barcode has been found
	if($BarcodeFound) {
		
		#Store the quality information of this read in a per barcode manner
		my $CharNumber=0;
		foreach my $Char (split //, $Qualities) {
 			$QualitiesByBarcode{$Barcode}->{$CharNumber} += ord($Char) unless ord($Char)==10;
 			$CharNumber++;
		}
		
		#Try to find the leading sequence and record any mistakes there. Determine an offset in case it is found
		my $LeadingOffset = 0;
		my $LeadingSequenceFound=0;
		my %DamerauResults;
		my $NotConvergedOffset=1;
		if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset), length($ExpectedLeadingSequence)) eq $ExpectedLeadingSequence ) {
			$Results{$Barcode}->[2]++; #Leading sequences found
			$Results{$Barcode}->[3]++; #Exact leading sequences found
			$LeadingSequenceFound=1;
		}
		else {
			$NotConvergedOffset=1;
			while ($NotConvergedOffset) {
				undef %DamerauResults;
				#DetermineDamerauLevenshteinDistance($ExpectedLeadingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset)),\%DamerauResults);
				if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdLeading) {
					$NotConvergedOffset=0;
					#Test if there is an insertion or deletion at the end of the leader sequence, which would indicate we need to analyze again +1 or -1, respectively
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
			#Store errors in leader sequence
			foreach my $Change (keys $DamerauResults{'Changes'}) {
				$LeadingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
			}
		} 
		
		#Try to find the trailing sequence and record any mistakes there. Determine an offset in case it is found
		my $InsertLength=$ExpectedInsertLength;
		my $TrailingSequenceFound=0;
		$NotConvergedOffset=1;
		while ($NotConvergedOffset) {
			undef %DamerauResults;
			if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)) eq $ExpectedTrailingSequence ) {
				$Results{$Barcode}->[4]++; #Trailing sequences found
				$Results{$Barcode}->[5]++; #Exact trailing sequences found
				#print "Exact trailing sequence found\n";
				$TrailingSequenceFound=1;
				$NotConvergedOffset=0;	
			}
			else {
				#DetermineDamerauLevenshteinDistance($ExpectedTrailingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)),\%DamerauResults);
				if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdTrailing) {
					$NotConvergedOffset=0;
					#Test if there is an insertion or deletion at the end of the leader sequence, which would indicate we need to analyze again +1 or -1, respectively
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
			#Store errors in trailing sequence
			foreach my $Change (keys $DamerauResults{'Changes'}) {
				$TrailingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
			}
		} 
		
		#Only continue if in any way a leading and trailing sequence are found
		if($LeadingSequenceFound && $TrailingSequenceFound) {
			$InsertLengths{$InsertLength}++;
			$Results{$Barcode}->[6]++;
			if($InsertLength==$ExpectedInsertLength) {
				$Results{$Barcode}->[7]++;			
#				my $LeadingSequence=substr( $Sequence, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset));
#				my $LeadingQualities=substr( $Qualities, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset));
				my $InsertSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)),$InsertLength);
#				my $InsertQualities=substr($Qualities,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)),$InsertLength);
#				my $TrailingSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)+$InsertLength));
#				my $TrailingQualities=substr($Qualities,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)+$InsertLength));
#				print color("red") . $LeadingSequence;
#				print color("green") . $InsertSequence;
#				print color("blue") . $TrailingSequence;
#				print color("reset") . "\n";
				if($Library{$InsertSequence}) {
					$Results{$Barcode}->[8]++;
					$InsertCounts{$InsertSequence}->{$Barcode}++;
				}
			}
		}
		else {
			print NOTANALYZED "$Sequence\n";
		}
	}		
}

#Output the insert size distribution
print "Writing insert size distribution to report\n";
print REPORT "Insert sizes\n";
foreach my $InsertLength (sort {$a <=> $b} keys %InsertLengths) {
	print REPORT "Size\t" . ($InsertLength) . "\t" . $InsertLengths{$InsertLength} . "\n";
}
print REPORT "\n";

#Output the leading error information
print "Writing leading error information to report\n";
print REPORT "Leading errors\nPosition\tInsertions\tDeletions\tMutations\n";
foreach my $Position (sort {$a <=> $b} keys %LeadingErrors) {
	my $TotalErrors=0;
	print REPORT "$Position\t";
	my @ErrorTypes = ('Insertion', 'Deletion', 'Mutation');
	foreach my $ErrorType (@ErrorTypes) {
		$TotalErrors=$TotalErrors+$LeadingErrors{$Position}->{$ErrorType};
		print REPORT $LeadingErrors{$Position}->{$ErrorType} . "\t";
	}
	print REPORT $TotalErrors . "\n";
}
print REPORT "\n";

#Output the trailing error information
print "Writing trailing error information to report\n";
print REPORT "Trailing errors\nPosition\tInsertions\tDeletions\tMutations\n";
foreach my $Position (sort {$a <=> $b} keys %TrailingErrors) {
	my $TotalErrors=0;
	print REPORT "$Position\t";
	my @ErrorTypes = ('Insertion', 'Deletion', 'Mutation');
	foreach my $ErrorType (@ErrorTypes) {
		$TotalErrors=$TotalErrors+$TrailingErrors{$Position}->{$ErrorType};
		print REPORT $TrailingErrors{$Position}->{$ErrorType} . "\t";
	}
	print REPORT $TotalErrors . "\n";
}
print REPORT "\n";

#Print report of read filtering per barcode
print "Writing filtering steps per barcode\n";
my @Totals;
print REPORT "Barcode\tBarcodes Found\tPerfect Barcodes Found\tLeader found\tPerfect leaders found\tTrailers found\tPerfect trailers found\tLeader and trailes found\tCorrect insert length\tInserts in library\n";
foreach my $Barcode ( @Barcodes ) {
	print REPORT "$Barcode\t";
	for ( my $counter = 0 ; $counter <= 8 ; $counter++ ) {
		print REPORT sprintf( "%i", $Results{$Barcode}->[$counter] ) . "\t";
		$Totals[$counter] += $Results{$Barcode}->[$counter];
	}
	if ( $Results{$Barcode}->[0] > 0 ) {
		print REPORT sprintf( "%4.2f%%", 100 * ( $Results{$Barcode}->[8] / $Results{$Barcode}->[0] ) ) . "\n";
	}
	else {
		print REPORT sprintf( "%4.2f%%", 0 ) . "\n";
	}
}
print REPORT sprintf( "TOTAL\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", $Totals[0], $Totals[1], $Totals[2], $Totals[3], $Totals[4], $Totals[5], $Totals[6], $Totals[7], $Totals[8]);

print REPORT sprintf(
	"\n\t\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%",
	100 * $Totals[1] / $Totals[0],
	100 * $Totals[2] / $Totals[0],
	100 * $Totals[3] / $Totals[0],
	100 * $Totals[4] / $Totals[0],
	100 * $Totals[5] / $Totals[0],
	100 * $Totals[6] / $Totals[0],
	100 * $Totals[7] / $Totals[0],
	100 * $Totals[8] / $Totals[0],
);
print REPORT "\n\nAnalyzed reads\t" . sprintf( "%i", $RecordsAnalyzed ) . "\n\n";

#Print quality information
print "Writing quality information to report\n";
print REPORT "Quality Information\n";
foreach my $Barcode ( @Barcodes ) {
	print REPORT "$Barcode\t";
	foreach my $CharNumber (sort {$a<=>$b} keys %{$QualitiesByBarcode{$Barcode}}) {
		print REPORT sprintf("%4.2f", $QualitiesByBarcode{$Barcode}->{$CharNumber} / $Results{$Barcode}->[0]) . "\t";
	}
	print REPORT "\n";
}

#Output the read counts of library hits
print "Writing filtered library insert counts\n";
print OUTPUT "LibraryID\tSequence";
foreach my $Barcode (@Barcodes) {
	print OUTPUT "\t$Barcode";
}
print OUTPUT "\n";
foreach my $InsertSequence (sort {$Library{$a} cmp $Library{$b}} keys %Library) {
	print OUTPUT $Library{$InsertSequence} . "\t" . $InsertSequence;
	foreach my $Barcode (@Barcodes) {
		if($InsertCounts{$InsertSequence}->{$Barcode}) {
			print OUTPUT "\t" . $InsertCounts{$InsertSequence}->{$Barcode};	
		}
		else {
			print OUTPUT "\t0";
		}
	}
	print OUTPUT "\n";
}

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
