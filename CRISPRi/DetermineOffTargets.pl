use Getopt::Std;
use warnings;
use strict;
print "Usage:perl $0 -i InputFile [-o OutputFile (default:InputFile.{degree})] -d Degree of off-targets to look for [-l flag:Skip looking at the results of the previous degree]\n";

my %opts;
my %RelativesFound;
my %NearGeneRelativesFound;
my %OutputText;
my %GenePositions;
my $ProcessRelatives = 0;
my $AdditionalIdenticalTargetsFound;
my $Margin = 300;
my $RefSeqFile = '../refseq/hg19.txt';

#Parse command parameters
getopt( 'oidl', \%opts );
die "ERROR in $0: No Inputfile given.\n" unless my $InputFile = $opts{'i'};
die "ERROR in $0: No depth given" unless my $Depth = $opts{'d'};
my $OutputFile;
$OutputFile=$InputFile . "." . $Depth unless $OutputFile = $opts{'o'};
my $SkipLookingAtLast = $opts{'l'};

#Load in all gene positions and make an array of that
open (IN, $RefSeqFile) or die "ERROR in $0: Cannot open refseqfile $RefSeqFile\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @RefSeqValues = split( /\t/, $Line );
	my $Chromosome = $RefSeqValues[2];
	my $GeneStart = $RefSeqValues[4];
	my $GeneEnd = $RefSeqValues[4];
	push(@{$GenePositions{$Chromosome}},[$GeneStart,$GeneEnd]);
}
close (IN) or die "ERROR in $0: Cannot close refseqfile $RefSeqFile\n";

#Create a file with all relatives of the unique target sites with a certain distance
#Read in the sam file, skip records with additional identical targets, skip records for which relatives were found in the previous run, unless l flag is set 
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
#Skip header
my $Header=<IN>;
chomp($Header);
my $PotentialRelativesFile = $InputFile . ".PotentialRelatives";
my $RelativesMatchFile = $InputFile . ".PotentialRelatives.matched";
while (defined(my $Line=<IN>)) {
	chomp($Line);
	my @LineValues = split( /\t/, $Line );
	my $NumberOfColumns = (scalar @LineValues);
	my $TargetSequence = $LineValues[6];
	$OutputText{$TargetSequence}=$Line;
	my $FoundInPreviousRun = 0;
	if ($NumberOfColumns>=11) {
		$FoundInPreviousRun = $LineValues[$NumberOfColumns - 2];
	}
	$RelativesFound{$TargetSequence} = -1;
	$NearGeneRelativesFound{$TargetSequence} = -1;
	if(!$FoundInPreviousRun || ($SkipLookingAtLast && $FoundInPreviousRun > -1) || ($FoundInPreviousRun == 1 && $Depth == 1)) { 
		$ProcessRelatives = 1;
		$RelativesFound{$TargetSequence} = 0;
		$NearGeneRelativesFound{$TargetSequence} = 0;
		`perl FindRelativesOfOligo.pl -i $TargetSequence -o $PotentialRelativesFile -d $Depth`;
	}
}
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";

#Map all the relatives to the genome
if ($ProcessRelatives) {
	`/media/Data/iKRUNC/bowtie2-2.1.0/bowtie2 /media/Data/iKRUNC/hg19-index/hg19 -f $PotentialRelativesFile -t --no-hd --score-min L,-5,0 -a -S $RelativesMatchFile -p 3`;
	#Read in everything that was mapped
	open (IN, $RelativesMatchFile) or die "ERROR in $0: Cannot open relatives match file $RelativesMatchFile\n";
	while (defined(my $Line=<IN>)) {
		chomp($Line);
		my @SAMValues = split( /\t/, $Line );
		if(!($SAMValues[1] & 4)) {
			my $ProtospacerSequence = $SAMValues[0];
			my $TargetSequence=substr($ProtospacerSequence,13);
			$RelativesFound{$TargetSequence}++;
			my $Chromosome=$SAMValues[2];
			my $Position=$SAMValues[3];
			if(defined(@{$GenePositions{$Chromosome}})) {
				my @GenePositionsOnChromosome=@{$GenePositions{$Chromosome}};
				foreach my $GenePositionOnChromosome (@GenePositionsOnChromosome) {
					my $GeneStart=@{$GenePositionOnChromosome}[0]-$Margin;
					my $GeneEnd=@{$GenePositionOnChromosome}[1]+$Margin;
					if($Position>=$GeneStart && $Position<=$GeneEnd) {
						$NearGeneRelativesFound{$TargetSequence}++;
					}
				}
			}
		}
	}
	close (IN) or die "ERROR in $0: Cannot close relatives match file $RelativesMatchFile\n";
}
unlink $RelativesMatchFile;
unlink $PotentialRelativesFile;

#While parsing the input file, write the output file, with an extra column indicating how many relatives were found and one indicating how many of those were found near genes
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
open (OUT, ">", $OutputFile) or die "ERROR in $0: Cannot open outpufile $OutputFile\n";
<IN>; #Skip header line
print OUT $Header . "\t$Depth MM off-targets\t$Depth MM off-targets near genes\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @LineValues = split( /\t/, $Line );
	print OUT $Line . "\t";
	print OUT $RelativesFound{$LineValues[6]} . "\t";
	print OUT $NearGeneRelativesFound{$LineValues[6]} . "\n";
}
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";
