#!/usr/bin/perl
#----------------------------------------------------------------------
#Structure of ReadResults:
#Initially filled in ParseInputFile:
#{$Barcode}->{TargetSequence}->[0]=Number of reads found with that particular TargetSequence in that Barcode
#{$Barcode}->{TargetSequence}->[1]=Number of Ns found in sequences with that particular TargetSequence in that Barcode

#Filled in FindExactMappingsAndCreateDistanceInputFile
#{$Barcode}->{TargetSequence}->[2]="mapped" "unmapped" depending on whether an exact match of the TargetSequence exists in the library

#Filled in ProcessDistances
#{$Barcode}->{TargetSequence}->[3]=The distance to the library gRNA sequence in {$Barcode}->{TargetSequence}->[4]
#{$Barcode}->{TargetSequence}->[4]=The library gRNA sequence this particular TargetSequence in this barcode is mapped to
#!!!ATTENTION!!! Since distance2Bins is run for every bin vs. library gRNA combination, it is in theory possible that these
#values are set multiple times for the same bin, in which case [4] and [5] will be arrays of distances/gRNAs

#In ProduceResultsData, The following is changed
#{$Barcode}->{TargetSequence}->[2] will be set to "ambiguous" when the sequence has more than one gRNA it most closely maps to.
#{$Barcode}->{TargetSequence}->[2] will be set to "distance" when the sequence has exactly one gRNA it most closely maps to.
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#Structure of ResultsData
#{$Barcode}->{$LibrarySequence}->[X] Number of reads found for that library sequence with distance X.
#----------------------------------------------------------------------

use Getopt::Long;
use File::Basename;
use warnings;
use strict;

sub ParseInputFile($$$$$$);
sub ReadgRNALibFile($);
sub FindExactMappingsAndCreateDistanceInputFile($$$);
sub CalculateDistances($$$$);
sub ProcessDistances($$);
sub ProduceResultsData($$);
sub PrintResults($$$$);
sub PrintAmbiguousData($$);
sub PrintNucleotideNDistribution($$$);

#Set parameters
my $InputFile;
my $gRNALibFile;
my $OutputDirectory	= "./";
my $gRNAOffset		= 30;	# Offset of gRNA cassette
my $gRNALength		= 18;	# Length of gRNA target sequences
my $BarcodeOffset	= 0;	# Offset of subexperiment barcode
my $BarcodeLength	= 6;	# Length of subexperiment barcode
my $DistanceThreshold	= 2;	# Threshold for minimum distance (Hamming or Damerau-Levenshtein)
my $DistanceCalculator	= "/media/Data/LIB007/barcodedistanceshort";    # Location of the file containing C distance calculations (Hamming and Damerau-Levenshtein)

#Initialize
system("clear");
print "Script usage:gRNAMapper.pl -InputFile -gRNALibFile [-OutputDirectory -DistanceThreshold]\n";
print "Started analysis at \t" . localtime() . "\n";

GetOptions(
	"InputFile=s"       => \$InputFile,
	"gRNALibFile=s"    => \$gRNALibFile,
	"OutputDirectory=s" => \$OutputDirectory,
	"DistanceThreshold=i" => \$DistanceThreshold
);

die "gRNAMapper.pl -InputFile -gRNALibFile [-OutputDirectory -DistanceThreshold]" if !$InputFile;

if ( !( $OutputDirectory =~ /\/$/ ) ) {
	$OutputDirectory .= "/";
}
open( DIR, $OutputDirectory ) or die "Output directory is not accessible\n";
close(DIR);

my $OutputPrefix = $OutputDirectory . basename $InputFile;
print $OutputPrefix . "\n";

#Parse the input file to fill the data structures
my ( $ReadResults, $NucleotideNDistribution ) = ParseInputFile( $InputFile, $gRNAOffset, $gRNALength, $BarcodeOffset, $BarcodeLength, $gRNALibFile );

#Read in the gRNA Library
my $gRNALibrary = ReadgRNALibFile($gRNALibFile);

#Then, find those library gRNAs that map exactly to the read results and set their mapping state to "mapped" and
#create a distance input file for those reads that could not be exactly mapped to an gRNA from the library. Set those initially to "unmapped"
FindExactMappingsAndCreateDistanceInputFile( $OutputPrefix, $gRNALibrary, $ReadResults );

#Then use that file to create a file with reads that are less far than a threshold distance from any gRNA from the library (those could be multiple!)
#Distance is currently calculated using Hamming distance, which just assumes substitutions. To also assume insertions/deletions, recompile the barcodedistance.c file
#after changing the distance calculation function into the Damerau-Levenshtein one
CalculateDistances( $DistanceCalculator, $gRNALibFile, $DistanceThreshold, $OutputPrefix );

#The the created distance file is used to further fill the ReadResults structure with Library gRNA sequence and mapping distance for those reads that mapped
#with a lower than threhold distance. NOTE!!! Since one read could in principle map to multiple library gRNA sequences, [4] and [5] entries of ReadResults
#may contain multiple values that are not in any particular order!!!
ProcessDistances( $ReadResults, $OutputPrefix );

#Next, Result structures are produced. Mapping to library gRNAs is consolidated. If there is one library gRNA mapping with shortest distance, mapping
#state will be set to "distance", if there's multiple with shortest distance, it will be set to "ambiguous" and it will remain "unmapped" otherwise.
#The structure of ResultsData is ResultsData{$Barcode}->{$Sequence}->[X]->[N] with X being number of reads for that sequence found with distance X
#N=0 contains gRNA hits, N=1 contains all the UMIs
my ( $ResultsData, $AmbiguousData ) = ProduceResultsData( $ReadResults, $gRNALibrary );

#Write outputfiles
PrintResults( $ResultsData, $gRNALibrary, $OutputPrefix, $DistanceThreshold );
PrintAmbiguousData( $AmbiguousData, $OutputPrefix );
PrintNucleotideNDistribution( $NucleotideNDistribution, $OutputPrefix, $gRNALength );

#THE END
print "Finished analysis at \t" . localtime() . "\n";

sub ParseInputFile($$$$$$$) {
	print "\n\n***ParseInputFile***\n";
	my ( $InputFile, $gRNAOffset, $gRNALength, $BarcodeOffset, $BarcodeLength, $gRNALibFile ) = @_;
	print STDERR "Parsing input file $InputFile\n";
	my %Results      = ();
	my %MismatchBias = ();
	my $NumberOfRecordsProcessed = 0;
	open( IN, $InputFile ) or die "Could not access file $InputFile\n";
	while ( defined( my $Line = <IN> ) ) {
		$NumberOfRecordsProcessed++;
		if ( !( $NumberOfRecordsProcessed % 100000 ) ) {
			print "Parsing record $NumberOfRecordsProcessed\n";
		}
		chomp($Line);
		my $Sequence = $Line;
		# Extract barcode
		my $Barcode = substr( $Sequence, $BarcodeOffset, $BarcodeLength );
		#Obtain the gRNA target sequence
		my $TargetSequence = substr( $Sequence, $gRNAOffset, $gRNALength );
		# increment hit counter for that target sequence
		$Results{$Barcode}->{$TargetSequence}->[0]++;
		# store position of Ns for each barcode
		# and number of Ns for each target sequence
		my $NCounter = 0;
		while ( $TargetSequence =~ /N/g ) {
			$NCounter++;
			$Results{$Barcode}->{$TargetSequence}->[1]++;
			$MismatchBias{$Barcode}->{ pos($TargetSequence) }++;
		}
	}
	close(IN) or die "Could not close file $InputFile\n";
	print "Processed $NumberOfRecordsProcessed records for file $InputFile\n";
	return ( \%Results, \%MismatchBias );
}

sub ReadgRNALibFile($) {
	print "\n\n***ReadgRNALibFile***\n";

	#Read in an gRNA Library file
	#The format of this file should be [ID] tab [SEQUENCE]
	my ($gRNALibFile) = @_;
	print STDERR "Reading library file $gRNALibFile\n";
	my %gRNALibrary = ();
	my $gRNAsFound;
	open( gRNA, $gRNALibFile ) or die "Could not read $gRNALibFile.\n";
	while ( defined( my $Line = <gRNA> ) ) {
		$gRNAsFound++;
		chomp($Line);
		my @values = split( /\t/, $Line );
		$gRNALibrary{ (substr($values[1],0,length($values[1])-2)) } = $values[0];
		print "" . (substr($values[1],0,length($values[1])-2)) . "\n\n";
	}
	close(gRNA) or die "Could not close file $gRNALibFile\n";
	print "$gRNAsFound gRNAs found in library file $gRNALibFile\n";
	return \%gRNALibrary;
}

sub FindExactMappingsAndCreateDistanceInputFile($$$) {
	print "\n\n***FindExactMappingsAndCreateDistanceInputFile***\n";
	my ( $OutputPrefix, $gRNALibrary, $ReadResults ) = @_;
	my $DistanceInputFile = $OutputPrefix . ".distance.in.csv";
	print STDERR "Writing file for distance compute $DistanceInputFile\n";
	my $HitCounter    = 0;
	my $RecordCounter = 0;
	open( DISTANCEINPUT, ">", $DistanceInputFile )
	  or die "Could not open $DistanceInputFile for writing.\n";
	foreach my $Barcode ( keys %$ReadResults ) {

		foreach my $TargetSequence ( keys %{ $ReadResults->{$Barcode} } ) {
			$RecordCounter++;
			if ( !exists $gRNALibrary->{$TargetSequence} ) {
				$ReadResults->{$Barcode}->{$TargetSequence}->[2] = "unmapped";
				print DISTANCEINPUT "$Barcode\t$TargetSequence\t" . $ReadResults->{$Barcode}->{$TargetSequence}->[0] . "\n";
			}
			else {
				$ReadResults->{$Barcode}->{$TargetSequence}->[2] = "mapped";
				$HitCounter++;
			}
		}
	}
	close(DISTANCEINPUT) or die "Could not close $DistanceInputFile\n";
	my $LibraryHitRatio = "";
	if($RecordCounter > 0) {
		$LibraryHitRatio = sprintf( "%3.2f", 100 * ( $HitCounter / $RecordCounter ) );	
	}
	print "Found $HitCounter exact library matches out of $RecordCounter target sequences extracted (" . $LibraryHitRatio . "%).\n";
}

sub CalculateDistances($$$$) {
	print "\n\n***CalculateDistances***\n";
	my ( $DistanceCalculator, $gRNALibFile, $DistanceThreshold, $OutputPrefix ) = @_;
	my $DistanceInputFile  = $OutputPrefix . ".distance.in.csv";
	my $DistanceOutputFile = $OutputPrefix . ".distance.out.csv";
	my @args               = ( $DistanceCalculator, $DistanceInputFile, $gRNALibFile, $DistanceOutputFile, $DistanceThreshold );
	print "Running distance calculation command: $DistanceCalculator $DistanceInputFile $gRNALibFile $DistanceOutputFile $DistanceThreshold\n";
	my $DistanceCalculationResult = system(@args);
	print "Finished $DistanceInputFile with exit code\t" . $DistanceCalculationResult . "\n";
}

sub ProcessDistances($$) {
	print "\n\n***ProcessDistances***\n";
	my ( $ReadResults, $OutputPrefix ) = @_;
	my $DistanceOutputFile = $OutputPrefix . ".distance.out.csv";
	open( DIST, $DistanceOutputFile )
	  or die "Could not read distance file $DistanceOutputFile\n";
	while ( defined( my $Line = <DIST> ) ) {
		chomp($Line);

		# values[0] -> Barcode
		# values[1] -> TargetSequence
		# values[2] -> Number of hits
		# values[3] -> Library sequence
		# values[4] -> Library sequence identifier
		# values[5] -> Distance
		my @values = split( /\t/, $Line );
		push( @{ $ReadResults->{ $values[0] }->{ $values[1] }->[3] }, $values[5] );    # Distance
		push( @{ $ReadResults->{ $values[0] }->{ $values[1] }->[4] }, $values[3] );    # gRNA
	}
	close(DIST) or die "Could not close file $DistanceOutputFile\n";
}

sub ProduceResultsData($$) {
	print "\n\n***ProduceResultsData***\n";
	my ( $ReadResults, $gRNALibrary ) = @_;
	my $ResultsData;
	my %AmbiguousData = ();
	foreach my $Barcode ( keys %$ReadResults ) {
		foreach my $TargetSequence ( keys %{ $ReadResults->{$Barcode} } ) {

			# if a read is a direct library hit, add to results
			if ( $gRNALibrary->{$TargetSequence} ) {
				$ResultsData->{$TargetSequence}->{$Barcode}->[0] = $ReadResults->{$Barcode}->{$TargetSequence}->[0];
			}

			# otherwise, investigate indirect hits:
			else {

				#If there is distance information available for this read, it means it has mapped below threshold to an gRNA of the library. So let's determine which gRNA it should be mapped to
				if ( $ReadResults->{$Barcode}->{$TargetSequence}->[3] ) {

					#In principle, each read could have mapped to more than one library gRNA, in which case there are more than one distances in $ReadResults->{$Barcode}->{$TargetSequence}->[3]
					#If of these multiple distances, there is one shortest, then we'll assign that library gRNA to our read, otherwise we'll make it ambiguous
					my @Distances = @{ $ReadResults->{$Barcode}->{$TargetSequence}->[3] };
					my %Hash      = ();
					foreach my $Distance (@Distances) {
						$Hash{$Distance}++;
					}

					#The hash hash contains an array of how many times (value) a certain distance (key) was present
					#Sorting this hash according to keys would result in a @sorted vector with distances. $sorted[0] would thus be the shortest occurring distance
					my @Sorted = ( sort keys %Hash );
					for ( my $i = 0 ; $i < @Distances ; $i++ ) {
						if ( $Distances[$i] == $Sorted[0] ) {

							#If there is more than one shortest distance a read is mapped to, make it ambiguous
							#By doing this in a loop over all distances, there will be one entry in the ambiguousData structure for every gRNA/read combination
							if ( $Hash{ $Sorted[0] } > 1 ) {
								push( @{ $AmbiguousData{$Barcode}->{$TargetSequence} }, ( $Distances[$i] . "\t" . $ReadResults->{$Barcode}->{$TargetSequence}->[4]->[$i] ) );
								$ReadResults->{$Barcode}->{$TargetSequence}->[2] = "ambiguous";
							}

							#If there is only one gRNA that maps closest to the read, then assign the read to that gRNA
							else {
								my $MappedDistance  = $ReadResults->{$Barcode}->{$TargetSequence}->[3]->[$i];
								my $LibrarySequence = $ReadResults->{$Barcode}->{$TargetSequence}->[4]->[$i];
								$ResultsData->{$LibrarySequence}->{$Barcode}->[$MappedDistance] += $ReadResults->{$Barcode}->{$TargetSequence}->[0];
								$ReadResults->{$Barcode}->{$TargetSequence}->[2] = "distance";
							}
						}
					}
				}

				#If there is no distance information available for this read, it will remain unmapped
				else {
					$ResultsData->{$TargetSequence}->{$Barcode}->[0] = $ReadResults->{$Barcode}->{$TargetSequence}->[0];
				}
			}
		}
	}
	return ( $ResultsData, \%AmbiguousData );
}

sub PrintResults($$$$) {
	print "\n\n***PrintResults***\n";
	my ( $ResultsData, $gRNALibrary, $OutputPrefix, $DistanceThreshold) = @_;
	my $ResultsFile  = $OutputPrefix . ".results.csv";
	my %SummaryStats = ();
	open( RES, ">", $ResultsFile )
	  or die "Could not create $ResultsFile results file\n";

	#Discover all barcodes
	my %Barcodes = ();
	foreach my $TargetSequence ( keys %$ResultsData ) {
		foreach my $Barcode ( keys %{ $ResultsData->{$TargetSequence} } ) {
			$Barcodes{$Barcode}++;
		}
	}

	#Write header
	print RES "Library_ID\tSequence";
	foreach my $Barcode ( sort keys %Barcodes ) {
		print RES "\t$Barcode #Exact Reads\t$Barcode #Total Reads\t$Barcode fraction perfect reads";
	}
	print RES "\n";

	#Write data
	foreach my $TargetSequence ( keys %$ResultsData ) {
		my $String = "";
		if ( $gRNALibrary->{$TargetSequence} ) {
			$String .= $gRNALibrary->{$TargetSequence} . "\t";
		}
		else {
			$String .= "**Unmapped**\t";
		}
		$String .= "$TargetSequence";

		foreach my $Barcode ( sort keys %Barcodes ) {
			my $TotalHits = 0;
			my $ExactHits = 0;
			if ( $ResultsData->{$TargetSequence}->{$Barcode}->[0] ) {
				$String .= "\t" . $ResultsData->{$TargetSequence}->{$Barcode}->[0];
				$ExactHits = $ResultsData->{$TargetSequence}->{$Barcode}->[0];
			}
			else {
				$String .= "\t0";
			}

			for ( my $Distance = 0 ; $Distance <= $DistanceThreshold ; $Distance++ ) {
				if ( $ResultsData->{$TargetSequence}->{$Barcode}->[$Distance] ) {
					$TotalHits += $ResultsData->{$TargetSequence}->{$Barcode}->[$Distance];
				}
			}
			if ( $TotalHits > 0 ) {
				$String .= "\t$TotalHits\t " . $ExactHits / $TotalHits;
			}
			else {
				$String .= "\t0\t0";
			}
		}

		# print total hits
		$String .= "\n";
		#Skip unmapped reads
		if ($gRNALibrary->{$TargetSequence}) {
			print RES $String;
		}
	}
	close(RES) or die "Could not close file $ResultsFile\n";
}

sub PrintAmbiguousData($$) {
	print "\n\n***PrintAmbiguousData***\n";
	my ( $AmbiguousHits, $OutputPrefix ) = @_;
	my $OutputFile = $OutputPrefix . ".ambiguous.csv";
	open( AMB, ">", $OutputFile ) or die "Could not create $OutputFile ambiguous hits file\n";
	print AMB "Barcode\tSequence\tDistance\tgRNA\n";
	foreach my $Barcode ( sort keys %$AmbiguousHits ) {
		foreach my $TargetSequence ( keys %{ $AmbiguousHits->{$Barcode} } ) {
			foreach my $AmbiguousHit ( @{ $AmbiguousHits->{$Barcode}->{$TargetSequence} } ) {
				print AMB "$Barcode\t$TargetSequence\t$AmbiguousHit\n";
			}
		}
	}
	close(AMB) or die "Could not close file $OutputFile\n";
}

sub PrintNucleotideNDistribution($$$) {
	print "\n\n***PrintNucleotideNDistribution***\n";
	my ( $NucleotideNDistribution, $OutputPrefix, $gRNALength ) = @_;
	my $NucleotideNDistributionFile = $OutputPrefix . ".nucleotideNdistribution.csv";
	open( my $OutputFileHandle, ">", $NucleotideNDistributionFile )
	  or die "$0: Could not open $NucleotideNDistributionFile $!";
	foreach my $Barcode ( sort keys %$NucleotideNDistribution ) {
		for ( my $Position = 1 ; $Position <= $gRNALength ; $Position++ ) {
			$NucleotideNDistribution->{$Barcode}->{$Position}
			  ? ( print $OutputFileHandle "$Barcode\t$Position\t" . $NucleotideNDistribution->{$Barcode}->{$Position} . "\n" )
			  : ( print $OutputFileHandle "$Barcode\t$Position\t0\n" );
		}
	}
	close($OutputFileHandle)
	  or die "could not close $NucleotideNDistributionFile.\n";
}
