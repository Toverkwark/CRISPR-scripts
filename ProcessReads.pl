sub ProcessReads($$$$$$$$$) {
	my ($InputFile,$BarcodeLength,$BarcodeOffset,$ExpectedInsertLength,$ExpectedLeadingSequence,$ExpectedTrailingSequence,$ErrorThresholdLeading,$ErrorThresholdTrailing,$Library) = @_;
	my %LeadingErrors;
	my %TrailingErrors;
	my %InsertLengths;
	my %QualitiesByBarcode;
	my $BarcodeFound=0;
	my $BarcodeFoundExact=0;
	my $LeadingSequenceFound=0;
	my $LeadingSequenceFoundExact=0;
	my $TrailingSequenceFound=0;
	my $TrailingSequenceFoundExact=0;
	my $LeadingAndTrailingFound=0;
	my $InsertCorrectLength=0;
	my $InsertInLibrary=0;
	my $InsertSequence;
	my @Barcodes;
	@Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
	

	open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
	
	while ( defined( my $line = <INPUT> ) ) {
		$RecordsAnalyzed++;
		if ( !( $RecordsAnalyzed % 100 ) ) {
			print "Analyzing record $RecordsAnalyzed of inputfile $InputFile\n";
		}
		my $Sequence = <INPUT>;
		my $line3 = <INPUT>;
		my $Qualities = <INPUT>;
		chomp($Sequence);
		chomp($Qualitiges);
		#Get the barcode. See if it exists. If not, try to map it with maximally 1 nucleotide replacement and only 1 match existing.
		my $Barcode = substr( $Sequence, $BarcodeOffset, $BarcodeLength );
		if ( grep( /$Barcode/, @Barcodes ) ) {
			$BarcodeFound=1;
			$BarcodeFoundExact=1;
		}	
		else {
			my $MatchedBarcode = MatchBarcode( $Barcode, @Barcodes );
			if ($MatchedBarcode) {
				$Barcode = $MatchedBarcode;
				$BarcodeFound=1;
			}
		}

		#Only continue to analyze this read if a valid barcode has been found
		if($BarcodeFound) {	
			#Try to find the leading sequence and record any mistakes there. Determine an offset in case it is found		
			my $LeadingOffset = 0;
			my %DamerauResults;
			my $NotConvergedOffset=1;
			if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset), length($ExpectedLeadingSequence)) eq $ExpectedLeadingSequence ) {
				$LeadingSequenceFound=1;
				$LeadingSequenceFoundExact=1;
			}
			else {
				$NotConvergedOffset=1;
				while ($NotConvergedOffset) {
					undef %DamerauResults;
					DetermineDamerauLevenshteinDistance($ExpectedLeadingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset), (length($ExpectedLeadingSequence)+$LeadingOffset)),\%DamerauResults);
					if($DamerauResults{'AccuratelyDetermined'} && $DamerauResults{'Distance'}<=$ErrorThresholdLeading) {
						$NotConvergedOffset=0;
						#Test if there is an insertion or deletion at the end of the leader sequence, which would indicate we need to analyze again +1 or -1, respectively
						foreach my $Change (keys $DamerauResults{'Changes'}) {
							if($Change==length($ExpectedLeadingSequence)) {
								if($DamerauResults{'Changes'}->{$Change} eq 'Insertion') {
									$LeadingOffset--;
									$NotConvergedOffset=1;
								}
								if($DamerauResults{'Changes'}->{$Change} eq 'Deletion') {
									$LeadingOffset++;
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
				#Store errors in leader sequence
				foreach my $Change (keys $DamerauResults{'Changes'}) {
					lock(%LeadingErrors);
					$LeadingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
				}
			} 
			
			#Try to find the trailing sequence and record any mistakes there. Determine an offset in case it is found
			my $InsertLength=$ExpectedInsertLength;
			$NotConvergedOffset=1;
			my $Trimming=0;
			while ($NotConvergedOffset) {
				undef %DamerauResults;
				if ( substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), length($ExpectedTrailingSequence)) eq $ExpectedTrailingSequence ) {
					$TrailingSequenceFound=1;
					$TrailingSequenceFoundExact=1;
					$NotConvergedOffset=0;	
				}
				else {
					DetermineDamerauLevenshteinDistance($ExpectedTrailingSequence,substr( $Sequence, ($BarcodeLength+$BarcodeOffset+$LeadingOffset+$InsertLength+length($ExpectedLeadingSequence)), (length($ExpectedTrailingSequence)+$Trimming)),\%DamerauResults);
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
							if($Change==length($ExpectedTrailingSequence) && $DamerauResults{'Changes'}->{$Change} eq 'Insertion') {
								$Trimming--;
								$NotConvergedOffset=1;
							}
							if($Change==length($ExpectedTrailingSequence) && $DamerauResults{'Changes'}->{$Change} eq 'Deletion') {
								$Trimming++;
								$NotConvergedOffset=1;
							}
							#Break off searching for trailer if needing to search past the sequence length
							if($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)+$InsertLength+$Trimming+length($ExpectedTrailingSequence)>150) {
								$NotConvergedOffset=0;
								$DamerauResults{'AccuratelyDetermined'}=0;
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
				#Store errors in trailing sequence
				foreach my $Change (keys $DamerauResults{'Changes'}) {
					lock(%TrailingErrors);
					$TrailingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
				}
			} 
			$Results{$Barcode}->[0]++;
			$Results{$Barcode}->[1]++ if ($BarcodeFoundExact);
			
			#Store the quality information of this read in a per barcode manner
			my $CharNumber=0;
			foreach my $Char (split //, $Qualities) {
 				lock(%QualitiesByBarcode);
 				$QualitiesByBarcode{$Barcode}->{$CharNumber} += ord($Char) unless ord($Char)==10;
 				$CharNumber++;
 			}
 			if ($LeadingSequenceFound) {
 				$Results{$Barcode}->[2]++;
 				$Results{$Barcode}->[3]++ if ($LeadingSequenceFoundExact);	
 			} 
 			if ($TrailingSequenceFound) {
 				$Results{$Barcode}->[4]++;
 				$Results{$Barcode}->[5]++ if ($TrailingSequenceFoundExact);
 			}
 			
 			if($LeadingSequenceFound && $TrailingSequenceFound) {
				lock(%InsertLengths);
				$InsertLengths{$InsertLength}++;
				$Results{$Barcode}->[6]++;
				if($InsertLength==$ExpectedInsertLength) {
					$Results{$Barcode}->[7]++;			
					$InsertSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)),$InsertLength);
					lock(%$Library);
					if($$Library{$InsertSequence}) {
						$Results{$Barcode}->[8]++;
						$InsertCounts{$InsertSequence}->{$Barcode}++;
					}	
				}
			}
			else {
				$NotAnalyzed=$NotAnalyzed . "$Sequence\n";
			}
		}			
	}
	close(INPUT)  or die "Could not close input file $InputFile.\n";
	return();
}

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

1;