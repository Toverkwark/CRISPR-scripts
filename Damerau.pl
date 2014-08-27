sub DetermineDamerauLevenshteinDistance($$$);
use strict;

my $First =  "abcccdefghijkl";
my $Second = "abcdefghijk";
my %ReturnHash;
DetermineDamerauLevenshteinDistance( $First, $Second, \%ReturnHash );
if ($ReturnHash{'AccuratelyDetermined'}) {
	print "Total Distance:" . $ReturnHash{'Distance'} . "\n";
	foreach my $Change (keys $ReturnHash{'Changes'}) {
		print "$Change\t" . $ReturnHash{'Changes'}->{$Change} . "\n";
	}
}
else {
	print "Could not accurately determine match\n";	
} 

sub DetermineDamerauLevenshteinDistance ($$$) {
	my ( $e, $f, $ReturnHash ) = @_;
	#print "$e\n$f\n";
	my @First = split( "", $e );
	my @Second = split( "", $f );

	my $Rows = length($e);
	my $Columns = length($f);
	my @d;
	my $Row;
	my $Column;
	
	for ( my $k = 0 ; $k <= $Rows ; $k++ ) {
		$d[$k][0] = $k;
	}
		for ( my $k = 0 ; $k <= $Columns ; $k++ ) {
		$d[0][$k] = $k;
	}

	for ( $Row = 1 ; $Row <= $Rows ; $Row++ ) {
		for ( $Column = 1 ; $Column <= $Columns ; $Column++ ) {
			my $cost = 1;
			$cost = 0 if ( $First[$Row-1] eq $Second[$Column-1]);
			$d[$Row][$Column] = $d[$Row][$Column-1]+1;
			$d[$Row][$Column] = ($d[$Row-1][$Column]+1) if (($d[$Row-1][$Column]+1) < $d[$Row][$Column]);
			$d[$Row][$Column] = ($d[$Row-1][$Column-1] + $cost) if (($d[$Row-1][$Column-1] + $cost) < $d[$Row][$Column]);

		}
		
	}

 	#Read back mutations
 	my $Distance = $d[$Rows][$Columns];
 	my $CurrentDistance=$Distance;
 	my $Up;
 	my $Left;
 	my $Diagonal;
 	$Row=$Rows;
 	$Column=$Columns;
 	$$ReturnHash{'AccuratelyDetermined'}=1;
 	do {
 		$Up=$d[$Row-1][$Column];
 		$Left=$d[$Row][$Column-1];
 		$Diagonal=$d[$Row-1][$Column-1];
 		if($Diagonal==$CurrentDistance && $Left>=$Diagonal && $Up>=$Diagonal) {
 			$Row--;
 			$Column--;
 		}
 		else {
 			if($Diagonal<$CurrentDistance && $Diagonal<$Up && $Diagonal<$Left) {
 				$$ReturnHash{'Changes'}->{$Row}="Mutation";
 				#print "mutation at nucleotide " . $Row . "\n";	
 				$Row--;
 				$Column--;
 				$CurrentDistance--;
 			}
 			else {
 				if($Left<$CurrentDistance && $Left<$Diagonal && $Left<$Up) {
 					$$ReturnHash{'Changes'}->{$Row}="Insertion";
 					#print "insertion at nucleotide " . $Row . "\n";
 					$Column--;
 					$CurrentDistance--;
 				}
 				else {
 					if($Up<$CurrentDistance && $Up<$Diagonal && $Up <$Left) {
 						$$ReturnHash{'Changes'}->{$Row}="Deletion";
 						#print "deletion at nucleotide " . $Row . "\n";
 						$Row--;
 						$CurrentDistance--;
 					}
 					else {
 						$$ReturnHash{'AccuratelyDetermined'}=0;
 						#print "Could not accurately determine mutation\n";
 						$Row--;
 						$Column--;
 					}
 				}	
 			}
 		}
 	} while ($Row>0 && $Column>0);
 		
 	while ($Row>0) {
 		$$ReturnHash{'Changes'}->{$Row}="Deletion";
 		#print "deletion at nucleotide " . $Row . "\n";
 		$Row--
 	}
 	while ($Column>0) {
 		$$ReturnHash{'Changes'}->{$Row}="Insertion";				
 		print "insertion at nucleotide " . ($Row+1) . "\n";
 		$Column--
 	} 
 	
 	$$ReturnHash{'Distance'}=$Distance;
	#print "Total distance is $Distance\n";
}
