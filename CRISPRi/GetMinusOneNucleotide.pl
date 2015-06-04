require 'FetchGenomicSequence.pl';
sub FetchGenomicSequence($$$);

open(IN,'library.csv');
<IN>;
while(defined(my $Line=<IN>)) {
	chomp($Line);
	@values=split(/\t/,$Line);
	$Sense=$values[5];
	$Chr=$values[7];
	$Length=$values[6];
	$Sequence=$values[4];
	$Orientation=$values[8];
	$TSS=$values[10];
	$Position=$values[11];
	if($Chr ne "#N/A") {
		if(($Orientation eq "-" && $Sense eq "Sense") || ($Orientation eq "+" && $Sense eq "Antisense")) {
			if($Orientation eq "-") {
				$PrintText=FetchGenomicSequence($Chr,$TSS-$Position+1,$TSS-$Position+4+length($Sequence));
			} else {
				$PrintText=FetchGenomicSequence($Chr,$TSS+$Position,$TSS+$Position+3+length($Sequence));
			}
			$PrintText=~tr/ACTGactg/TGACTGAC/;
			$PrintText=reverse($PrintText);
			print $PrintText . "\t$Sequence\t$MinusOne\n";
		}
		else {	
			if($Orientation eq "-") {
				$PrintText=FetchGenomicSequence($Chr,$TSS-$Position-length($Sequence)-2,$TSS-$Position+1);
			} else {
				$PrintText=FetchGenomicSequence($Chr,$TSS+$Position-3-length($Sequence),$TSS+$Position);
			}
			print $PrintText . "\t$Sequence\n";
		}
	}
}
close(IN);
