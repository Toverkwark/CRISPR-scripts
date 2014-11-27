use Getopt::Std;
use warnings;
use strict;
require 'FetchGenomicSequence.pl';
sub FetchGenomicSequence($$$);

#Cas9 cuts here:
#XXXXXXXXXXXXXXXXX|XXXNGG
#In this Script, the last G of the NGG PAM sequence is considered the actual target in terms of locations in the genome

my $TSSUpstream=50; #Use this number of nt upstream of any TSS for protospacer searching
my $TSSDownstream=300; #Use this number of nt downstream of any TSS for protospacer searching
my $MinimalGuideLength=19;
my $MaximalGuideLength=22;
my $PAMLength=3;

my %opts;
my %TargetSites;

getopt( 'i', \%opts );
open(IN, my $InputFile=$opts{'i'}) or die "ERROR in $0:Could not open inputfile\n";

#Loop over all refseq IDs in inputfile and fetch search sequences
print "Parsing input file, looking up TSS surroundings\n";
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
			$Sequence=FetchGenomicSequence($Chromosome,$GeneStart-$MaximalGuideLength-$PAMLength-$TSSUpstream+1,$GeneStart+$MaximalGuideLength+$PAMLength+$TSSDownstream-1) . "\n";
		}
		else {
			$GeneStart    = $RefSeqValues[5];
			$GeneEnd      = $RefSeqValues[4]+1;
			$Sequence=FetchGenomicSequence($Chromosome,$GeneStart-$MaximalGuideLength-$PAMLength-$TSSDownstream+1,$GeneStart+$MaximalGuideLength+$PAMLength+$TSSUpstream-1) . "\n";
			$Sequence =~ tr/ACTG/TGAC/;
			$Sequence = reverse($Sequence);
		}
		$TargetSites{$Gene}->{$RefSeq}->{'Sequence'}=$Sequence;
		$TargetSites{$Gene}->{$RefSeq}->{'TSS'}=$GeneStart;
		$TargetSites{$Gene}->{$RefSeq}->{'Orientation'}=$Orientation;
		$TargetSites{$Gene}->{$RefSeq}->{'Chromosome'}=$Chromosome;
	}
}
close(IN) or die "Could not close inputfile\n";
	
#Perform identification of valid protospacers	
print "Identifying valid protospacers\n";
foreach my $Gene (keys %TargetSites) {
#	print "Working on gene $Gene\n";
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		my @ValidGuides;
#		print "\tRefSeq $RefSeq\n";
		my $SearchSequence = $TargetSites{$Gene}->{$RefSeq}->{'Sequence'};
		#Loop over guide length in search of guides
		for(my $GuideLength=$MinimalGuideLength;$GuideLength<=$MaximalGuideLength;$GuideLength++) {
			#First search in the sense direction
			for(my $Position=$MaximalGuideLength-$GuideLength;$Position<length($SearchSequence)-($GuideLength+$PAMLength-1)-($MaximalGuideLength+$PAMLength);$Position++) {
				my $Guide=substr($SearchSequence,$Position,$GuideLength+$PAMLength);
				if(substr($Guide,-2) eq 'GG') {
					my $GuideSequence = substr($Guide,0,$GuideLength);
					if(substr($GuideSequence,0,1) eq 'G' || $GuideLength==20) {
						my $ValidGuide={};
						#Store position relative to TSS
						my $GuidePosition=-$TSSUpstream+$Position-($MaximalGuideLength-$GuideLength);
						$ValidGuide->{'Position'}=$GuidePosition;
						$ValidGuide->{'GuideSequence'}=$GuideSequence;
						$ValidGuide->{'Sense'}='Sense';
						$ValidGuide->{'GuideLength'}=$GuideLength;
						$ValidGuide->{'Extra'}=0;
						$ValidGuide->{'Extra'}=1 if(substr($GuideSequence,0,1) ne 'G');
						$ValidGuide->{'Selected'}=1;
						push(@ValidGuides,$ValidGuide);
					}
				}
			}
						
			#Then search in the anstisense direction.
			for(my $Position=($MaximalGuideLength+$PAMLength-1);$Position<length($SearchSequence)-($MaximalGuideLength+$PAMLength-1);$Position++) {
				my $Guide=substr($SearchSequence,$Position,$GuideLength+$PAMLength);
				if(substr($Guide,0,2) eq 'CC') {
					
					my $GuideSequence=substr($Guide,-$GuideLength);
					$GuideSequence =~ tr/ACTG/TGAC/;
					$GuideSequence = reverse($GuideSequence);
					
					if(substr($GuideSequence,0,1) eq 'G' || $GuideLength==20) {
						my $ValidGuide={};
						#Store position relative to TSS
						my $GuidePosition=-$TSSUpstream+$Position-($MaximalGuideLength+$PAMLength-1);
						$ValidGuide->{'Position'}=$GuidePosition;
						$ValidGuide->{'GuideSequence'}=$GuideSequence;
						$ValidGuide->{'Sense'}='Antisense';
						$ValidGuide->{'GuideLength'}=$GuideLength;
						$ValidGuide->{'Extra'}=0;
						$ValidGuide->{'Extra'}=1 if(substr($GuideSequence,0,1) ne 'G');
						$ValidGuide->{'Selected'}=1;
						push(@ValidGuides,$ValidGuide);
					}
				}
			}				
		}
		$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}=\@ValidGuides;
	}
}

#Perform guide filtering and convert data structure into guide sequence centered output hash
my %GuidesOfGene;
print "Protospacer filtering\n";
foreach my $Gene (keys %TargetSites) {
#	print "Working on gene $Gene\n";
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		my @TargetSitesOfRefSeq=@{$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}};
		foreach my $TargetSiteOfRefSeq (@TargetSitesOfRefSeq) {
			#Filter out guides that have 4 Ts in a row as this may result in premature transcriptional termination
			$TargetSiteOfRefSeq->{'Selected'} = 0 if($TargetSiteOfRefSeq->{'GuideSequence'} =~ /TTTT/);
			$GuidesOfGene{$TargetSiteOfRefSeq}->{$Gene}=1;						
		}
	}
}
foreach my $Guide (keys %GuidesOfGene) {
	if(keys $GuidesOfGene{$Guide} > 1) {
		#If a guide is present in more than one gene, find all the guides of this sequence and deselect them
		foreach my $Gene (keys $GuidesOfGene{$Guide}) {
			foreach my $RefSeq (keys $TargetSites{$Gene}) {
				my @TargetSitesOfRefSeq=@{$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}};
				foreach my $TargetSiteOfRefSeq (@TargetSitesOfRefSeq) {
					$TargetSiteOfRefSeq->{'Selected'}=0 if($Guide eq $TargetSiteOfRefSeq->{'GuideSequence'});
				}
			} 
		}
	}
}

#Output to file
print "Writing output file\n";
my $OutputFile=$InputFile . '.protospacers';
open (OUT, ">", $OutputFile) or die "Cannot open outpufile\n";
print OUT "Gene\tRefSeq\tChromosome\tTSS\tOrientation\tGuide Position\tGuide Sequence\tGuide sense\tGuide length\tGuide Extra\n";
			
foreach my $Gene (keys %TargetSites) {
	foreach my $RefSeq (keys $TargetSites{$Gene}) {
		my @TargetSitesOfRefSeq=@{$TargetSites{$Gene}->{$RefSeq}->{'TargetSites'}};
		foreach my $TargetSiteOfRefSeq (sort {$a->{'Position'} <=> $b->{'Position'}} @TargetSitesOfRefSeq) {
			if ($TargetSiteOfRefSeq->{'Selected'}) {
				print OUT $Gene . "\t" . $RefSeq . "\t";
				print OUT $TargetSites{$Gene}->{$RefSeq}->{'Chromosome'} . "\t";
				print OUT $TargetSites{$Gene}->{$RefSeq}->{'TSS'} . "\t";		
				if($TargetSites{$Gene}->{$RefSeq}->{'Orientation'}) {
					print OUT "+\t";
				}
				else {
					print OUT "-\t";
				} 		
				print OUT $TargetSiteOfRefSeq->{'Position'} . "\t";
				print OUT $TargetSiteOfRefSeq->{'GuideSequence'} . "\t";
				print OUT $TargetSiteOfRefSeq->{'Sense'} . "\t";
				print OUT $TargetSiteOfRefSeq->{'GuideLength'} . "\t";
				print OUT $TargetSiteOfRefSeq->{'Extra'} . "\n";
			}
		}
	}
}
close (OUT) or die "ERROR in $0: Cannot close outputfile\n";