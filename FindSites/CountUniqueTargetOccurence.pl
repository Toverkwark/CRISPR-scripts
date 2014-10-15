use Getopt::Long;
use warnings;
use strict;

sub IsPositionInGene($$);

our %Genes;
my %Protospacers;
my $InputFile;
GetOptions(
	"input=s"  => \$InputFile,
);

#Read in RefSeq information
my $RefSeqFile="../refseq/hg19.txt";
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


my $CountFile = $InputFile . ".counts";
my $OutputFile = $InputFile . ".out";
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
open (COUNT, ">", $CountFile) or die "ERROR in $0: Cannot open outpufile $CountFile\n";
open (OUT, ">", $OutputFile) or die "ERROR in $0: Cannot open outpufile $OutputFile\n";
my $RecordsAnalyzed=0;
while (defined(my $Line = <IN>)) {
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 1000000 ) ) {
		print "Analyzing record $RecordsAnalyzed of inputfile $InputFile\n";
	}
	chomp($Line);
	my @SAMValues = split( /\t/, $Line );
	if(!($SAMValues[1] & 4)) {
		my $ProtospacerSequence = $SAMValues[9];
		if($SAMValues[1] & 16) {
			$ProtospacerSequence =~ tr/ACTG/TGAC/;
			$ProtospacerSequence = reverse($ProtospacerSequence);
		}
		my $TargetSequence=substr($ProtospacerSequence,0,20);
		$Protospacers{$TargetSequence}->[0]++;
		if (IsPositionInGene(substr($SAMValues[2], 3),($SAMValues[3] + ($SAMValues[1] & 16 ? 6 : 17)))) {
			print OUT $TargetSequence . "\t" . $Line . "\t1\n";
			$Protospacers{$TargetSequence}->[1]++;
		}
		else {
			print OUT $TargetSequence . "\t" . $Line . "\t0\n";	
		}
	}
}
foreach my $TargetSequence (keys %Protospacers) {
	#Output protospacers with counts
	print COUNT $TargetSequence . "\t" . $Protospacers{$TargetSequence}->[0];
	if($Protospacers{$TargetSequence}->[1]) {
		print COUNT "\t" . $Protospacers{$TargetSequence}->[1] . "\n";
	}
	else {
		print COUNT "\t0\n";
	}
		 
}
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";
close (COUNT) or die "ERROR in $0: Cannot close outputfile $CountFile\n";
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";

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