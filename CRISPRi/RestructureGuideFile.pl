se Getopt::Std;
use warnings;
use strict;

my %opts;
my %TargetSites;

getopt( 'i', \%opts );
open(IN, my $InputFile=$opts{'i'}) or die "ERROR in $0:Could not open inputfile\n";

#Loop over all guides in inputfile
while (defined (my $Line=<IN>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
}

close(IN) or die "Could not close inputfile\n";
	
#Preparing output structure
my %OutputStructure;
foreach my $Gene (keys %TargetSites) {
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		my @TargetSitesOfRefSeq=@{$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}};
		foreach my $TargetSiteOfRefSeq (@TargetSitesOfRefSeq) {
			if ($TargetSiteOfRefSeq->{'Selected'}) {
				my $GuideSequence=$TargetSiteOfRefSeq->{'GuideSequence'};
				push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'RefSeqs'}}, $RefSeq);
				push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'TSSs'}}, $TargetSites{$RefSeq}->{'TSS'});
				push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'Orientation'}}, $TargetSites{$RefSeq}->{'Orientation'});
				$OutputStructure{$Gene}->{$GuideSequence}->{'Position'}=$TargetSiteOfRefSeq->{'Position'};
				$OutputStructure{$Gene}->{$GuideSequence}->{'Sense'}=$TargetSiteOfRefSeq->{'Sense'};
				$OutputStructure{$Gene}->{$GuideSequence}->{'GuideLength'}=$TargetSiteOfRefSeq->{'GuideLength'};
				$OutputStructure{$Gene}->{$GuideSequence}->{'Extra'}=$TargetSiteOfRefSeq->{'Extra'};
			}
		}
	}
}