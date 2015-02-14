use lib '..';
use LocalSettings;
my %LocalSettings=getconfig();
my $HumanGenome=$LocalSettings{'HumanGenome'};
sub FetchGenomicSequence($$$) {
	my($Chromosome, $StartSite, $EndSite) = @_;
	my $GenomicSequence;
	open (my $InputHandle, "$HumanGenome/chr$Chromosome.stripped.fa") or die "Could not open genomic sequence file $HumanGenome/chr$Chromosome.stripped.fa\n";
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
