use Getopt::Std;
use warnings;
use strict;
require 'FetchGenomicSequence.pl';

sub FetchGenomicSequence($$$);

#Cas9 cuts here:
#XXXXXXXXXXXXXXXXX|XXXNGG

my $TSSUpstream=50; #Use this number of nt upstream of any TSS for protospacer searching
my $TSSDownstream=300; #Use this number of nt downstream of any TSS for protospacer searching
my $MinimalGuideLength=19;
my $MaximalGuideLength=22;

my %opts;
my %TargetSites;

getopt( 'i', \%opts );
open(IN, my $InputFile=$opts{'i'}) or die "ERROR in $0:Could not open inputfile\n";

#Loop over all refseq IDs in inputfile and fetch search sequences
while (defined (my $Line=<IN>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
	my $Gene=$LineValues[0];
	my $RefSeq=$LineValues[1];
	if(substr($RefSeq,0,2) eq 'NM') {
		my $RefSeqInfo = `grep -P "$RefSeq\t" ../refseq/hg19.txt`;
		die "ERROR in $0: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
		my @RefSeqValues = split( /\t/, $RefSeqInfo );
		my $Chromosome = substr( $RefSeqValues[2], 3 );
		#The refGene.txt file uses zero based chromosomal locations for the start sites
		my $GeneStart;
		my $GeneEnd;
		my $Sequence;
		my $Orientation = 0;
		if($RefSeqValues[3] eq '+') {
			$Orientation = 1;
			$GeneStart    = $RefSeqValues[4]+1;
			$GeneEnd      = $RefSeqValues[5];
			$Sequence=FetchGenomicSequence($Chromosome,$GeneStart-17-$TSSUpstream,$GeneStart+17+$TSSDownstream) . "\n";
		}
		else {
			$GeneStart    = $RefSeqValues[5];
			$GeneEnd      = $RefSeqValues[4]+1;
			$Sequence=FetchGenomicSequence($Chromosome,$GeneStart-17-$TSSDownstream,$GeneStart+17+$TSSUpstream) . "\n";
		}
		$TargetSites{$Gene}->{$RefSeq}->{'Sequence'}=$Sequence;
		$TargetSites{$Gene}->{$RefSeq}->{'TSS'}=$GeneStart;
		$TargetSites{$Gene}->{$RefSeq}->{'Orientation'}=$Orientation;
		$TargetSites{$Gene}->{$RefSeq}->{'Chromosome'}=$Chromosome;
	}
}
close(IN) or die "Could not close inputfile\n";
	
#Perform identification of valid protospacers	
foreach my $Gene (keys %TargetSites) {
	print "Working on gene $Gene\n";
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		my @ValidTargetSites;
		print "\tRefSeq $RefSeq\n";
		my $SearchSequence = $TargetSites{$Gene}->{$RefSeq}->{'Sequence'};
		print $SearchSequence . "\n";
		#Loop over guide length in search of a Weissman-like guide
		for(my $GuideLength=$MinimalGuideLength;$GuideLength<=$MaximalGuideLength;$GuideLength++) {
			#First search in the sense direction, run until the 12th nt from the end
			my $FWSearchSequence=substr($SearchSequence,0,length($SearchSequence)-11);
			while ( $FWSearchSequence =~ /(?<=(.{($GuideLength+1)})(GG))/g ) {
					my $TargetSite=substr($1,0,$GuideLength);
					push(@ValidTargetSites,substr($1,0,$GuideLength)) if(substr($TargetSite,0,1) eq 'G' || $GuideLength==20);
			}
			
			#Then search in the anstisense direction. Start at the 12th nt and run until the end
			my $RVSearchSequence=substr($SearchSequence,11);	
			while ( $RVSearchSequence =~ /(?=((CC).{$GuideLength+1}))/g ) {
				my $TargetSite=substr($1,1,$GuideLength);
				$TargetSite =~ tr/ACTG/TGAC/;
				$TargetSite = reverse($TargetSite);
				push(@ValidTargetSites,$TargetSite) if(substr($TargetSite,0,1) eq 'G' || $GuideLength==20);
			}
			
			$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}=\@ValidTargetSites;	
		}
	}
}

#Output to file
my $OutputFile=$InputFile . '.protospacers';
open (OUT, ">", $OutputFile) or die "Cannot open outpufile\n";
			
foreach my $Gene (keys %TargetSites) {
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		foreach my $TargetSite (@{$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}}) {
			unless ($TargetSite =~ /TTTTT/) {
				print OUT $Gene . "\t" . $RefSeq . "\t";
				print OUT $TargetSites{$Gene}->{$RefSeq}->{'Chromosome'} . "\t";
				print OUT $TargetSites{$Gene}->{$RefSeq}->{'TSS'} . "\t";
				print OUT $TargetSites{$Gene}->{$RefSeq}->{'Orientation'} . "\t";
				print OUT $TargetSite . "\n";
			}
		}
	}
}

close (OUT) or die "ERROR in $0: Cannot close outputfile\n";