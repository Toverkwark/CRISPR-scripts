sub FetchGenomicSequence($$$) {
	my($Chromosome, $StartSite, $EndSite) = @_;
	my $GenomicSequence;
	open (IN, "/media/Data/iKRUNC/hg19/chr$Chromosome.stripped.fa") or die "Could not open genomic sequence file hg19/chr$Chromosome.stripped.fa\n";
	if($EndSite>=$StartSite)
	{
		seek(IN,$StartSite-1,1);
		read IN,$GenomicSequence,$EndSite-$StartSite+1;
	}
	close (IN);
	$GenomicSequence =~ tr/a-z/A-Z/;
	return $GenomicSequence;
}
1;
