use Getopt::Std;
use warnings;
use strict;

sub IsPositionInGene($$);

my $ScriptName="QualityScoreSites.pl";
my %opts;
my %RelativesFound;
my %ExonicRelativesFound;
my %OutputText;
our %Genes;
my $ProcessRelatives = 0;
my $AdditionalIdenticalTargetsFound;
my $Margin = 100;
my $RefSeqFile = '../refseq/hg19.txt';

getopt( 'oids', \%opts );
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
die "ERROR in $ScriptName: No Inputfile given.\n" unless my $InputFile = $opts{'i'};
my $Depth = $opts{'d'};
die "ERROR in $ScriptName: No seed length given.\n" unless my $SeedLength = $opts{'s'};
getopts('l',\%opts);
my $SkipLookingAtLast = $opts{'l'};

#Load in all exonic position and make a hash off that, extending the exons with a certain margin
open (IN, $RefSeqFile) or die "ERROR in $0: Cannot open refseqfile $RefSeqFile\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @RefSeqValues = split( /\t/, $Line );
	my $Chromosome = substr( $RefSeqValues[2], 3 );
	my $GeneStart = $RefSeqValues[4];
	my $GeneEnd = $RefSeqValues[5];
	push (@{$Genes{$Chromosome}}, [$GeneStart, $GeneEnd]); 
}
close (IN) or die "ERROR in $0: Cannot close refseqfile $RefSeqFile\n";

#Create a file with all relatives of the unique target sites with a certain distance
#Read in the sam file, skip records with additional identical targets, skip records for which relatives were found in the previous run, unless l flag is set 
open (IN, $InputFile) or die "ERROR in $ScriptName: Cannot open inputfile $InputFile\n";
my $PotentialRelativesFile = $InputFile . ".PotentialRelatives";
my $RelativesMatchFile = $InputFile . ".PotentialRelatives.matched";
while (defined(my $Line=<IN>)) {
	chomp($Line);
	my @SAMValues = split( /\t/, $Line );
	my $NumberOfColumns = (scalar @SAMValues);
	my $ProtospacerSequence = $SAMValues[10];
	if($SAMValues[2] & 16) {
		$ProtospacerSequence =~ tr/ACTG/TGAC/;
		$ProtospacerSequence = reverse($ProtospacerSequence);
	}
	my $TargetSequence=substr($ProtospacerSequence,0,20);
	$OutputText{$TargetSequence}=$Line;
	
	$ExonicRelativesFound{$TargetSequence}=0;
	my $FoundInPreviousRun = 0;
	if (!$SkipLookingAtLast) {
		$FoundInPreviousRun = $SAMValues[$NumberOfColumns - 2];
	}
	$RelativesFound{$TargetSequence} = -1;
	$ExonicRelativesFound{$TargetSequence} = -1;
	if(!$FoundInPreviousRun || ($SkipLookingAtLast && $FoundInPreviousRun > -1) || ($FoundInPreviousRun == 1 && $Depth == 1)) { 
		$ProcessRelatives = 1;
		$RelativesFound{$TargetSequence} = 0;
		$ExonicRelativesFound{$TargetSequence} = 0;
		`perl FindRelativesOfOligo.pl -i $TargetSequence -o $PotentialRelativesFile -d $Depth -s $SeedLength`;
	}
}
close (IN) or die "ERROR in $ScriptName: Cannot close inputfile $InputFile\n";

#Map all the relatives to the genome. Write everything that maps to a separate file
if ($ProcessRelatives) {
	`/media/Data/iKRUNC/bowtie2-2.1.0/bowtie2 /media/Data/iKRUNC/hg19-index/hg19 -f $PotentialRelativesFile -t --no-hd --score-min L,-5,0 -a -S $RelativesMatchFile`;	
	
	#Read in everything that was mapped
	open (IN, $RelativesMatchFile) or die "ERROR in $ScriptName: Cannot open relatives match file $RelativesMatchFile\n";
	open (OUT, ">", ($OutputFile . ".relatives")) or die "ERROR in $ScriptName: Cannot open outputfile " . ($OutputFile . ".relatives") . "\n";
	while (defined(my $Line=<IN>)) {
		chomp($Line);
		my @SAMValues = split( /\t/, $Line );
		if(!($SAMValues[1] & 4)) {
			my $ProtospacerSequence = $SAMValues[0];
			my $TargetSequence=substr($ProtospacerSequence,13,20);
			$RelativesFound{$TargetSequence}++;
			my $IsInGene=0;
			if (IsPositionInGene(substr($SAMValues[2], 3),($SAMValues[3] + ($SAMValues[1] & 16 ? 6 : 17)))) {
				$ExonicRelativesFound{$TargetSequence}++;
				$IsInGene=1;
			}
			print OUT $TargetSequence . "\t" . $Line . "\t" . $IsInGene . "\n";
		}
	}
	close (IN) or die "ERROR in $ScriptName: Cannot close relatives match file $RelativesMatchFile\n";
	close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile " . ($OutputFile . ".relatives") . "\n";
#	unlink $RelativesMatchFile;
#	unlink $PotentialRelativesFile;
}

#Write the output file, with an extra column indicating how many relatives were found and one indicating how many of those were found near exons
open (OUT, ">", $OutputFile) or die "ERROR in $ScriptName: Cannot open outputfile $OutputFile\n";
my $TargetsFound;
my $TargetsNotFound;
foreach my $TargetSequence (sort {$a cmp $b} keys %OutputText) {
	print OUT $OutputText{$TargetSequence} . "\t". $RelativesFound{$TargetSequence} . "\t" . $ExonicRelativesFound{$TargetSequence} . "\n";
}	
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";

sub IsPositionInGene($$) {
	my($Chromosome, $Position) = @_;
	my $IsPositionInGene = 0;
	foreach my $Gene (@{$Genes{$Chromosome}}) {
		my $StartSite = ${$Gene}[0];
		my $EndSite = ${$Gene}[1];
		
		if ($Position > $StartSite && $Position <= $EndSite) {
			$IsPositionInGene = 1;
			last;		
		}
	}
	return $IsPositionInGene;
}