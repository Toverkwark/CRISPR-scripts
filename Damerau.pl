sub DetermineDamerauLevenshteinDistance($$);
use strict;

my $First =  "abcdefghijkl";
my $Second = "abcdefghixjkl";
DetermineDamerauLevenshteinDistance( $First, $Second );

sub DetermineDamerauLevenshteinDistance ($$) {
	my ( $e, $f ) = @_;
	print "$e\n$f\n";
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
	#Print the table and fill info for finding back mutations
	my @IdenticalMinimalValues;
	my @Indices;
	my @Distances;
	my %MinimalDistanceTable;
	for ($Row=0;$Row<=$Rows;$Row++) {
		my %ColumnValues;
		for($Column=0;$Column<=$Columns;$Column++) {
			$ColumnValues{$Column}=$d[$Row][$Column];
			my $result = sprintf( "%02d", $d[$Row][$Column] );
			print " " . $result;
		}
		my @SortedColumnValues=sort {$a<=>$b} @{$d[$Row]};
		my $NumberOfIdenticalMinimalValues=0;
		foreach (@SortedColumnValues) {
			if ($_ == $SortedColumnValues[0]) {
				$NumberOfIdenticalMinimalValues++;
			}
			else {
				last;
			}
		}
		print "\n";
		my @SortedIndices=sort { @{$d[$Row]}[$a] <=> @{$d[$Row]}[$b]} 0..$#{$d[$Row]};
		$MinimalDistanceTable{$Row}->{'NumberOfIdenticalMinimalValues'}=$NumberOfIdenticalMinimalValues;
		$MinimalDistanceTable{$Row}->{'Index'}=$SortedIndices[0];
		$MinimalDistanceTable{$Row}->{'Distance'}=$SortedColumnValues[0];
	}

 	#Read back mutations
 	my $Distance = $d[$Rows][$Columns];
 	my $CurrentDistance=$Distance;
 	my $Up;
 	my $Left;
 	my $Diagonal;
 	$Row=$Rows;
 	$Column=$Columns;
 	print "\n\n";
 	do {
 		$Up=$d[$Row-1][$Column];
 		$Left=$d[$Row][$Column-1];
 		$Diagonal=$d[$Row-1][$Column-1];
 		print "Row:$Row\tColumn:$Column\tUp:$Up\tLeft:$Left\tDiagonal:$Diagonal\t$CurrentDistance\n";
 		
 		if($Diagonal<$Left && $Diagonal<$Up) {
 			if($Diagonal<$CurrentDistance) {
 				print "mutation at nucleotide " . $Row . "\n";
 				$CurrentDistance--;
 			}
 			$Row--;
 			$Column--;
 		}
 		else {
 			if($Left<$Up && $Diagonal!=$Left) {
 				print "insertion at nucleotide " . $Row . "\n";
 				$Column--;
 				$CurrentDistance--;
 			}
 			else {
 				if($Diagonal!=$Up && $Diagonal!=$Left && $Up!=$Left) {
 					print "deletion at nucleotide " . $Row . "\n";
 					$Row--;
 					$CurrentDistance--;
 				}
 				else {
 					print "complex mutation detected\n";
 					$Row=1;
 					$Column=1
 				}
 			}
		}	 		
 	} while ($Row>1 && $Column>1);
 	print "Row:$Row\tColumn:$Column\tUp:$Up\tLeft:$Left\tDiagonal:$Diagonal\t$CurrentDistance\n";
 		
 	while ($Row>1) {
 		print "deletion at nucleotide " . ($Row-1) . "\n";
 		$Row--
 	}
 	while ($Column>1) {
 		print "insertion at nucleotide " . $Row . "\n";
 		$Column--
 	} 
 	
	print "Total distance is $Distance\n";
}
