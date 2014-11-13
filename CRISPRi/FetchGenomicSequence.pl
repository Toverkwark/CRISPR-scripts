sub FetchGenomicSequence($$$) {
	my($Chromosome, $StartSite, $EndSite) = @_;
	my $GenomicSequence;
	open (my $InputHandle, "/home/bastiaan/data/genomestuff/hg19/chr$Chromosome.stripped.fa") or die "Could not open genomic sequence file hg19/chr$Chromosome.stripped.fa\n";
	if($EndSite>=$StartSite)
	{
		seek($InputHandle,$StartSite-1,1);
		read $InputHandle,$GenomicSequence,$EndSite-$StartSite+1;
	}
	close ($InputHandle);
	$GenomicSequence =~ tr/a-z/A-Z/;
	return $GenomicSequence;
}
1;
