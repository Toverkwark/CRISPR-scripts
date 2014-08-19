sub DetermineDamerauLevenshteinDistance($$);
use strict;

my $a = "abcdefghijkl";
my $b = "abxxdefxghij";
DetermineDamerauLevenshteinDistance( $a, $b );

sub DetermineDamerauLevenshteinDistance ($$) {
	my ( $e, $f ) = @_;
	print "$e\n$f\n";
	my @a = split( "", $e );
	my @b = split( "", $f );

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
			$cost = 0 if ( $a[$Row-1] eq $b[$Column-1]);
			my $Change_a = $d[$Row][$Column-1]+1;
			my $Change_b  = $d[$Row-1][$Column]+1;
			my $Change_c  = $d[$Row-1][$Column-1] + $cost;
			if ( ($Change_c < $Change_a) && ($Change_c < $Change_b) ) {
				$d[$Row][$Column] = $Change_c;
			}
			else {
				if ( ($Change_b < $Change_c) ) {
					$d[$Row][$Column] = $Change_b;
				}
				else {
					$d[$Row][$Column] = $Change_a;
				}
			}
		}
	}
	#Print de tabel
	for ($Row=0;$Row<=$Rows;$Row++) {
		for($Column=0;$Column<=$Columns;$Column++) {
			my $result = sprintf( "%02d", $d[$Row][$Column] );
			print " " . $result;
		}
		print "\n";
	}

	#Lees terug wat de veranderingen zijn
	my $MinimalScore  = $d[$Rows][$Columns];
	my $CurrentColumn = $Columns-1;
	for ( $Row = ($Rows-1); $Row >= 0 ; $Row-- ) {
		my @EventsDetectedOnRow;
		for ( $Column = 1 ; $Column <= $Columns ; $Column++ ) {
			if ( $d[$Row][$Column] < $MinimalScore ) {
				push( @EventsDetectedOnRow, $Column);
			}
		}
		if ( scalar(@EventsDetectedOnRow) == 1 ) {
			if ( $EventsDetectedOnRow[0] < $CurrentColumn ) {
				print "Insertion found at position " . ( $Row + 1 ) . "\n";
				$CurrentColumn--;
			}
			else {
				if ( $EventsDetectedOnRow[0] > $CurrentColumn ) {
					print "Deletion found at position " . ( $Row + 1 ) . "\n";
					$CurrentColumn++;
				}
				else {
					print "Mutation found at position " . ( $Row + 1 ) . "\n";
					$CurrentColumn--;
				}
			}
		}
		if(scalar (@EventsDetectedOnRow)>0) {
			$MinimalScore=$MinimalScore-1;	
		}
		$CurrentColumn--;
	}
	my $distance = $d[$Rows][$Columns];
	print "Total distance is $distance\n";
}