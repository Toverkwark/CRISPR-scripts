#!/bin/perl
open(IN,"../../combined_reads");
open(OUT,'>',"weeded");
while(defined($Line=<IN>)) {
	chomp($Line);
	my @Values=split(/\t/,$Line);	
	if($Values[0] ne 'ERBB2-2' && $Values[0] ne 'FGR-1' && $Values[0] ne 'HCK-3' && $Values[0] ne 'MAP3K19-2' && $Values[0] ne 'PAK7-2' && $Values[0] ne 'PRKD2-2' && $Values[0] ne 'PRKD3-2' && $Values[0] ne 'STK35-1' && $Values[0] ne 'TSSK4-2' && $Values[0] ne 'TSSK4-4' && $Values[0] ne 'ULK2-3' && $Values[1] ne 'Vector' && $Values[3] ne 'Vector') {	
		if(substr($Values[1],0,6) eq 'HGLibA') {
			$Values[1]='NonTargetingControl'
		}
		if(substr($Values[3],0,6) eq 'HGLibA') {
			$Values[3]='NonTargetingControl'
		}
		foreach $value (@Values) {
			print OUT $value . "\t";
		}
		print OUT "\n";
	}
}
close(IN);
close(OUT);
