$NumberOfCores=32;
my $RefSeqFile = '/home/NKI/b.evers/bastiaan/gRNAs/CRISPR/kinome/genelist.out';
open (IN, $RefSeqFile) or die;
$CurrentCore=0;
for (my $i=0;$i<$NumberOfCores;$i++) {
	$OutputFileHandle = 'OUT.' . $i;
	$OutputFile = 'kinome/KinomeCRISPRSearchChunk.' . $i;
	#print "Opening $OutputFileHandle with file $OutputFile\n";
	open ($OutputFileHandle, ">", $OutputFile) or die "ERROR: Cannot open outputfile $OutputFile\n";
}

while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @RefSeqValues = split( /\t/, $Line );
	$RefSeqID = $RefSeqValues[1];
#	unlink "/home/NKI/b.evers/bastiaan/gRNAs/kinome/qualityfiles/$RefSeqID.qualities.4";
	unless (-e '/home/NKI/b.evers/bastiaan/gRNAs/kinome/qualityfiles/$RefSeqID.qualities.4') {
		$CurrentCore = 0 if ($CurrentCore==$NumberOfCores);
		$OutputFileHandle = 'OUT.' . $CurrentCore;
		select $OutputFileHandle;
		print "./RunRefSeqKinome.sh $RefSeqID\n";
		$CurrentCore++;
	}
}
close (IN) or die "ERROR: Cannot close inputfile";
for (my $i=0;$i<$NumberOfCores;$i++) {
	$OutputFileHandle = 'OUT.' . $i;
	$OutputFile = 'epigenetics/EpigeneticsCRISPRSearchChunk.' . $CurrentCore;
	close ($OutputFileHandle) or die "ERROR: Cannot close outputfile $OutputFile\n";
}
