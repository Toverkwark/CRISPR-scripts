sub DetermineDamerauLevenshteinDistance($$);
#use strict;

my $First =  "abcdefghijkl";
my $Second = "abcdefghijl";
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
 	
 	print "\n\n";
 	foreach $Row (sort {$b <=> $a} keys %MinimalDistanceTable) {
 		print $MinimalDistanceTable{$Row}->{'NumberOfIdenticalMinimalValues'} . "\t";
 		print $MinimalDistanceTable{$Row}->{'Index'} . "\t";
 		print $MinimalDistanceTable{$Row}->{'Distance'} . "\n";
 	}
 	print "\n\n";
 		
	my $distance = $d[$Rows][$Columns];
	print "Total distance is $distance\n";
}
