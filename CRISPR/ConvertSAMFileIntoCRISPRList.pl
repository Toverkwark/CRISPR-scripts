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
my $RefSeqInfo = `grep $RefSeq /media/Data/iKRUNCH/refGene.txt`;
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
	my $Name = $ProtospacerSequence;
	
	#Determine the color and the score of any protospacer
	my $Score =0;
	my $ColorRed=0;
	my $ColorGreen=0;
	my $ColorBlue=0;
	
	$Protospacers{"$Chromosome\t " . ($Start -1) . "\t$End\t$Name\t$Score\t$Orientation\t " . ($Start-1) . "\t$End\t$ColorRed,$ColorGreen,$ColorBlue\n"} = $Score;
	$ProtospacerSequences{$ProtospacerSequence}=$Score;
}

my $NumberOfProtospacers = 0;
foreach my $Protospacer (sort {$Protospacers{$b} <=> $Protospacers{$a}} keys %Protospacers) {
	print OUT $Protospacer;
	$NumberOfProtospacers++;
	#last if ($NumberOfProtospacers >= $SelectNumberOfProtospacers);
}

$NumberOfProtospacers = 0;
foreach my $ProtospacerSequence (sort {$ProtospacerSequences{$b} <=> $ProtospacerSequences{$a}} keys %ProtospacerSequences) {
	print OUTSEQ $ProtospacerSequence . "\t" . $ProtospacerSequences{$ProtospacerSequence} . "\t$RefSeq\n";
	$NumberOfProtospacers++;
	#last if ($NumberOfProtospacers >= $SelectNumberOfProtospacers);
}

close (IN) or die "ERROR in $ScriptName: Cannot close inputfile $InputFile\n";
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";
close (OUTSEQ) or die "ERROR in $ScriptName: Cannot close outputfile with sequences " . $OutputFile . ".seq\n";
