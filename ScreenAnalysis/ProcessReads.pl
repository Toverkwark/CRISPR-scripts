sub ProcessReads($$$$$$$$$) {
	my ($InputFile,$BarcodeLength,$BarcodeOffset,$ExpectedInsertLength,$ExpectedLeadingSequence,$ExpectedTrailingSequence,$ErrorThresholdLeading,$ErrorThresholdTrailing,$Library,@Barcodes) = @_;
	my %LeadingErrors;
	my %TrailingErrors;
	my %InsertLengths;
	my %InsertCounts;
	my %PerfectInsertCounts;
	my %QualitiesByBarcode;
	my %Results;
	my $RecordsAnalyzed=0;
	my $NotAnalyzed;
		
	open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
	open(MAPPED, ">", ($InputFile . ".mapped")) or die ("Could not open outputfile " . ($InputFile . ".mapped") . "\n");
	open(PERFECTMAPPED, ">", ($InputFile . ".perfect.mapped")) or die ("Could not open outputfile " . ($InputFile . ".perfect.mapped") . "\n");

	while ( defined( my $line = <INPUT> ) ) {
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
		$RecordsAnalyzed++;
		if ( !( $RecordsAnalyzed % 100 ) ) {
			print "Analyzing record $RecordsAnalyzed of inputfile $InputFile\n";
		}
		my $Sequence = <INPUT>;
		my $line3 = <INPUT>;
		my $Qualities = <INPUT>;
		chomp($Sequence);
		chomp($Qualities);
		
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
				$InsertLengths{$InsertLength}++;
				$Results{$Barcode}->[6]++;
				$Results{$Barcode}->[9]++ if ($LeadingSequenceFoundExact && $TrailingSequenceFoundExact);
				#Turn off insert length assertion for pGIN libraries with variable insert sizes
				#if($InsertLength==$ExpectedInsertLength) {
					$Results{$Barcode}->[7]++;
					$Results{$Barcode}->[10]++ if ($LeadingSequenceFoundExact && $TrailingSequenceFoundExact);			
					$InsertSequence=substr($Sequence,($BarcodeLength+$BarcodeOffset+$LeadingOffset+length($ExpectedLeadingSequence)),$InsertLength);
					if($$Library{$InsertSequence}) {
						$Results{$Barcode}->[8]++;
						$InsertCounts{$InsertSequence}->{$Barcode}++;
						print MAPPED $$Library{$InsertSequence} . "\t" . $Barcode . "\n";
						if ($LeadingSequenceFoundExact && $TrailingSequenceFoundExact) {
							$Results{$Barcode}->[11]++;
							$PerfectInsertCounts{$InsertSequence}->{$Barcode}++;
							print PERFECTMAPPED $$Library{$InsertSequence} . "\t" . $Barcode ."\n"; 
						} else {
							print PERFECTMAPPED "ERROR:Leading and/or trailing sequence not perfect\n";
						}
					}
					else {
						print MAPPED "ERROR:Insert not mapped to library\n";
						print PERFECTMAPPED "ERROR:Insert not mapped to library\n";
					}
				#}
				#else {
				#	print MAPPED "ERROR:Insert length incorrect\n";
				#	print PERFECTMAPPED "ERROR:Insert length incorrect\n";
				#}
			}
			else {
				$NotAnalyzed=$NotAnalyzed . "$Sequence\n";
				print MAPPED "ERROR:No leading or trailing sequence found\n";
				print PERFECTMAPPED "ERROR:No leading or trailing sequence found\n";
			}
		}
		else {
			print MAPPED "ERROR:No barcode found\n";
			print PERFECTMAPPED "ERROR:No barcode found\n";
		}			
	}
	close(INPUT)  or die "Could not close input file $InputFile.\n";
	
	#Write tmp output to output file
	open(OUTPUT, ">", ($InputFile . ".tmp")) or die ("Could not open outputfile " . ($InputFile . ".tmp") . "\n");
	#Write Results
	print OUTPUT "***RESULTS***\n";
	foreach my $Barcode (keys %Results) {
		foreach my $Entry (keys $Results{$Barcode}) {
			print OUTPUT "$Barcode\t$Entry\t" . $Results{$Barcode}->[$Entry]  . "\n";
		}
	}
	
	print OUTPUT "***LEADING ERRORS***\n"; 
	foreach my $Position (keys %LeadingErrors) {
		foreach my $Change (keys $LeadingErrors{$Position}) {
			print OUTPUT "$Position\t$Change\t" . $LeadingErrors{$Position}->{$Change} . "\n";
		}
	}
	
	print OUTPUT "***TRAILING ERRORS***\n"; 
	foreach my $Position (keys %TrailingErrors) {
		foreach my $Change (keys $TrailingErrors{$Position}) {
			print OUTPUT "$Position\t$Change\t" . $TrailingErrors{$Position}->{$Change} . "\n";
		}
	}
	
	print OUTPUT "***INSERT LENGTHS***\n";
	foreach my $Length(keys %InsertLengths) {
		print OUTPUT "$Length\t" . $InsertLengths{$Length} . "\n";
	}
	
	print OUTPUT "***INSERT COUNTS***\n";
	foreach my $Sequence (keys %InsertCounts) {
		foreach my $Barcode (keys $InsertCounts{$Sequence}) {
			print OUTPUT "$Sequence\t$Barcode\t" . $InsertCounts{$Sequence}->{$Barcode} . "\n";
		}
	}
	
	print OUTPUT "***PERFECT INSERT COUNTS***\n";
	foreach my $Sequence (keys %PerfectInsertCounts) {
		foreach my $Barcode (keys $PerfectInsertCounts{$Sequence}) {
			print OUTPUT "$Sequence\t$Barcode\t" . $PerfectInsertCounts{$Sequence}->{$Barcode} . "\n";
		}
	}
	
	print OUTPUT "***QUALITIES BY BARCODE***\n";
	foreach my $Barcode (keys %QualitiesByBarcode) {
		foreach my $CharNum (keys $QualitiesByBarcode{$Barcode}) {
			print OUTPUT "$Barcode\t$CharNum\t" . $QualitiesByBarcode{$Barcode}->{$CharNum} . "\n";
		}
	}
	
	print OUTPUT "***RECORDS ANALYZED***\n";
	print OUTPUT $RecordsAnalyzed . "\n";
	print OUTPUT "***NOT ANALYZED***\n";
	print OUTPUT $NotAnalyzed;
	close(OUTPUT) or die "Could not close outputfile " . ($InputFile . ".tmp") . "\n";
	close(MAPPED) or die "Could not close outputfile " . ($InputFile . ".mapped") . "\n";
	close(PERFECTMAPPED) or die "Could not close outputfile " . ($InputFile . ".perfect.mapped") . "\n";
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
