package LocalSettings;
require Exporter;
use strict;
our @ISA=qw(Exporter);
our @EXPORT=qw(getconfig);
my %confighash=(
	"Bowtie" => "/home/bastiaan/data/genomestuff/bowtie2-2.1.0/bowtie2",
	"HumanGenome" => "/home/bastiaan/data/genomestuff/hg19",
	"IndexedHumanGenome" => "/home/bastiaan/data/genomestuff/hg19-index",
	);
	
sub getconfig {
	return %confighash;
}
