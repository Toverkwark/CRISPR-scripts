open (IN,'/home/NKI/b.evers/RefSeq/refGene.txt') or die "Could not open refGene file\n";
while (defined ($Line=<IN>)) {
	chomp($Line);
	@RefSeqValues = split (/\t/,$Line);
	$ID=$RefSeqValues[1];
	$RefSeqs{$ID}->[0]++;
	$RefSeqs{$ID}->[1] = $RefSeqs{$ID}->[1] . $Line . "\n";
}

open (OUT,,">",'/home/NKI/b.evers/RefSeq/refGene.unique.txt') or die "Could not open outputfile\n";
foreach $ID (keys %RefSeqs) {
	if($RefSeqs{$ID}->[0]>1) {
		$RefSeqsWithMoreThanOneHit++;
		print $ID . "\t has " . $RefSeqs{$ID}->[0] . " occurences in the refGene.txt list:\n";
#		print $RefSeqs{$ID}->[1];
	}
	else {
		print OUT $RefSeqs{$ID}->[1];
	}
}
print "In total, " . $RefSeqsWithMoreThanOneHit . " refseq IDs had multiple lines in the refGene file\n";

close (IN);
close (OUT);
