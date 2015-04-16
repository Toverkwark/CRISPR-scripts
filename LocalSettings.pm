package LocalSettings;
require Exporter;
use strict;
our @ISA=qw(Exporter);
our @EXPORT=qw(getconfig);
my %confighash;
my $Location='Home';
if($Location eq 'Home') {
	%confighash=(
		"NumberOfCoresToUse" => "7",
		"Bowtie" => "/home/bastiaan/data/genomestuff/bowtie2-2.1.0/bowtie2",
		"HumanGenome" => "/home/bastiaan/data/genomestuff/hg19",
		"IndexedHumanGenome" => "/home/bastiaan/data/genomestuff/hg19-index",
	);	
}
else {
	if($Location eq 'Work') {
		use lib "/home/NKI/b.evers/perl5/lib/perl5";
		%confighash=(
			"NumberOfCoresToUse" => "3",	
			"Bowtie" => "/media/Data/iKRUNC/bowtie2-2.1.0/bowtie2",
			"HumanGenome" => "/media/Data/iKRUNC/hg19",
			"IndexedHumanGenome" => "/media/Data/iKRUNC/hg19-index",
		);	
	}
}

sub getconfig {
	return %confighash;
}