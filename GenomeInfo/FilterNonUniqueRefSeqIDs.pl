open (IN,'hg19.old') or die "Could not open refGene file\n";
while (defined ($Line=<IN>)) {
	chomp($Line);
	@RefSeqValues = split (/\t/,$Line);
	#Filter out entries with aberrant chromosome entries	
	if ($RefSeqValues[2] =~ m/chr([XY]|\d{1,2})\Z/) {
		$ID=$RefSeqValues[1];
		$RefSeqs{$ID}->[0]++;
		$RefSeqs{$ID}->[1] = $RefSeqs{$ID}->[1] . $Line . "\n";
	}
}

open (OUT,,">",'hg19.new') or die "Could not open outputfile\n";
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
