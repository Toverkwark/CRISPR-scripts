use Getopt::Std;
use File::Basename;
use warnings;
use strict;

my $ScriptName="FilterAndSelect.pl";
my $SelectNumberOfProtospacers=10;
my %opts;
my %Protospacers;
my %ProtospacerSequences;
getopt( 'io', \%opts );
die "ERROR in $ScriptName: No inputfile given.\n" unless my $InputFile = $opts{'i'};
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
open (OUT, ">", $OutputFile) or die "ERROR in $ScriptName: Cannot open outputfile $OutputFile\n";
open (OUTSEQ, ">", $OutputFile . ".seq") or die "ERROR in $ScriptName: Cannot open outputfile with sequences " . $OutputFile . ".seq\n";

#Find position of the gene to make sure browser opens in right position
my $RefSeq = substr(basename($InputFile),0,index(basename($InputFile),'.'));
my $RefSeqInfo = `grep -P "$RefSeq\t" ~/RefSeq/refGene.txt`;
die "ERROR in $ScriptName: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
my @RefSeqValues = split( /\t/, $RefSeqInfo );
my $Chromosome = substr( $RefSeqValues[2], 3 );
my $GeneStart    = $RefSeqValues[4];
my $GeneEnd      = $RefSeqValues[5];

#Write browser and track information to the output file
print OUT "browser position chr" .$Chromosome . ":" . $GeneStart . "-" . $GeneEnd . "\n";
print OUT "browser full refGene\n";
print OUT "track name='CRISPRs' description='Valid CRISPR protospacers' visibility='pack' itemRgb='On' db='hg19'\n";

#Loop over the input file
open (IN, $InputFile) or die "ERROR in $ScriptName: Cannot open inputfile $InputFile\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	#For every line, extract the necessary information
	my @SAMValues = split( /\t/, $Line );
	my $ProtospacerSequence = $SAMValues[9];
	my $Orientation = '+';
	if($SAMValues[1] & 16) {
		$Orientation = '-';
		$ProtospacerSequence =~ tr/ACTG/TGAC/;
		$ProtospacerSequence = reverse($ProtospacerSequence);
	}
	$Chromosome = $SAMValues[2];
	my $Start = $SAMValues[3];
	my $End = $Start+22;
	my $NumberOfIdentical3PrimeTargets=$SAMValues[19];
	my $NumberOfIdentical3PrimeTargetsNearExons=$SAMValues[20];
	my $Degree=1;
	my $ClosestRelatives=$SAMValues[21];
	my $ClosestRelativesNearExons=$SAMValues[22];
	for (my $i=23;$i<=27;$i=$i+2) {
		unless ($ClosestRelatives > 0) {
			$ClosestRelatives=$SAMValues[$i];
			$ClosestRelativesNearExons=$SAMValues[$i+1];
			$Degree=$Degree+1;
		}
	}
	my $Name = $NumberOfIdentical3PrimeTargets . "(" . $NumberOfIdentical3PrimeTargetsNearExons . "):" . $Degree . "MM:" . $ClosestRelatives . "(" . $ClosestRelativesNearExons . ")";
	
	#Determine the color and the score of any protospacer
	my $Score =0;
	my $ColorRed=0;
	my $ColorGreen=0;
	my $ColorBlue=1;
	
	#The score of any protospacer will be determined as follows
	#There will be one point added for every Degree, meaning the degree below that did not have relatives
	#3 Points will be added for protospacers that do not have OTHER identical 3' targets
	#2 Points will be added for protospacers that have 1 OTHER identical 3' target, which is not near Exons
	#1 Points will be added for protospacers that have more than 1 OTHER identical 3' target unless any of them near Exons
	#0.0001 points will be subtracted for the number of relatives in the closest degree
	#0.01 points will be substracted for the number of relatives in the closest degree that are near Exons
	$Score=$Degree;
	$Score=$Score + 3 if ($NumberOfIdentical3PrimeTargets == 1);
	$Score=$Score + 2 if ($NumberOfIdentical3PrimeTargets == 2 && $NumberOfIdentical3PrimeTargetsNearExons == 0);
	$Score=$Score + 1 if ($NumberOfIdentical3PrimeTargets > 2 && $NumberOfIdentical3PrimeTargetsNearExons == 0);
	$Score=$Score - (0.0001 * $ClosestRelatives);
	$Score=$Score - (0.01 * $ClosestRelativesNearExons);
	
	#Now. filter the protospacers
	#Let's try first to not have protospacers that have 1st or 2nd degree relatives
#	if ($Degree>1) {
		$Protospacers{"$Chromosome\t " . ($Start -1) . "\t$End\t$Name\t$Score\t$Orientation\t " . ($Start-1) . "\t$End\t$ColorRed,$ColorGreen,$ColorBlue\n"} = $Score;
		$ProtospacerSequences{$ProtospacerSequence}=$Score;
#	}
}

my $NumberOfProtospacers = 0;
foreach my $Protospacer (sort {$Protospacers{$b} <=> $Protospacers{$a}} keys %Protospacers) {
	print OUT $Protospacer;
	$NumberOfProtospacers++;
#	last if ($NumberOfProtospacers >= $SelectNumberOfProtospacers);
}

$NumberOfProtospacers = 0;
foreach my $ProtospacerSequence (sort {$ProtospacerSequences{$b} <=> $ProtospacerSequences{$a}} keys %ProtospacerSequences) {
	print OUTSEQ $ProtospacerSequence . "\t" . $ProtospacerSequences{$ProtospacerSequence} . "\t$RefSeq\n";
	$NumberOfProtospacers++;
#	last if ($NumberOfProtospacers >= $SelectNumberOfProtospacers);
}

close (IN) or die "ERROR in $ScriptName: Cannot close inputfile $InputFile\n";
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";
close (OUTSEQ) or die "ERROR in $ScriptName: Cannot close outputfile with sequences " . $OutputFile . ".seq\n";
