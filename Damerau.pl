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
 				$Row--;
 				$Column--;
 				$CurrentDistance--;
 			}
 			else {
 				if($Left<$CurrentDistance && $Left<$Diagonal && $Left<$Up) {
 					$$ReturnHash{'Changes'}->{$Row}="Insertion";
 					$Column--;
 					$CurrentDistance--;
 				}
 				else {
 					if($Up<$CurrentDistance && $Up<$Diagonal && $Up <$Left) {
 						$$ReturnHash{'Changes'}->{$Row}="Deletion";
 						$Row--;
 						$CurrentDistance--;
 					}
 					else {
 						$$ReturnHash{'AccuratelyDetermined'}=0;
 						$Row--;
 						$Column--;
 					}
 				}	
 			}
 		}
 	} while ($Row>0 && $Column>0);
 		
 	while ($Row>0) {
 		$$ReturnHash{'Changes'}->{$Row}="Deletion";
 		$Row--
 	}
 	while ($Column>0) {
 		$$ReturnHash{'Changes'}->{$Row}="Insertion";				
 		$Column--
 	} 
 	
 	$$ReturnHash{'Distance'}=$Distance;
}
1;