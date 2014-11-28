use Getopt::Std;
use warnings;
use strict;

my %opts;
my %OutputStructure;
sub PrintArray($@);
getopt( 'i', \%opts );
open(IN, my $InputFile=$opts{'i'}) or die "ERROR in $0:Could not open inputfile\n";

#Read header
my $Header=<IN>;

#Loop over all guides in inputfile
while (defined (my $Line=<IN>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
	my $Gene=$LineValues[0];
	my $RefSeq=$LineValues[1];
	my $Chromosome=$LineValues[2];
	my $TSS=$LineValues[3];
	my $Orientation=$LineValues[4];
	my $Position=$LineValues[5];
	my $GuideSequence = $LineValues[6];
	my $Sense=$LineValues[7];
	my $GuideLength=$LineValues[8];
	my $IsExtra=$LineValues[9];
	my $IsUnique=$LineValues[10];
	
	push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'RefSeqs'}}, $RefSeq);
	push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'TSSs'}}, $TSS);
	push(@{$OutputStructure{$Gene}->{$GuideSequence}->{'Positions'}}, $Position);
	$OutputStructure{$Gene}->{$GuideSequence}->{'Chromosome'}=$Chromosome;
	$OutputStructure{$Gene}->{$GuideSequence}->{'Orientation'}=$Orientation;
	$OutputStructure{$Gene}->{$GuideSequence}->{'Sense'}=$Sense;
	$OutputStructure{$Gene}->{$GuideSequence}->{'GuideLength'}=$GuideLength;
	$OutputStructure{$Gene}->{$GuideSequence}->{'IsExtra'}=$IsExtra;
	$OutputStructure{$Gene}->{$GuideSequence}->{'IsUnique'}=$IsUnique;
}

my $OutputHandle;
open($OutputHandle,">",$InputFile . ".restructured") or die "ERROR in $0:Could not open outputfile for writing\n";
print {$OutputHandle} "Gene\tGuide Sequence\tSense\tGuide Length\tChromosome\tOrientation\tRefSeqs\tTSSs\tPositions\tIsExtra\tIsUnique\n";
foreach my $Gene (sort {$a cmp $b} keys %OutputStructure) {
	foreach my $GuideSequence (keys $OutputStructure{$Gene}) {
		print {$OutputHandle} $Gene . "\t";
		print {$OutputHandle} $GuideSequence . "\t";
		print {$OutputHandle} $OutputStructure{$Gene}->{$GuideSequence}->{'Sense'} . "\t";
		print {$OutputHandle} $OutputStructure{$Gene}->{$GuideSequence}->{'GuideLength'} . "\t";
		print {$OutputHandle} $OutputStructure{$Gene}->{$GuideSequence}->{'Chromosome'} . "\t";
		print {$OutputHandle} $OutputStructure{$Gene}->{$GuideSequence}->{'Orientation'} . "\t";
		PrintArray ($OutputHandle,@{$OutputStructure{$Gene}->{$GuideSequence}->{'RefSeqs'}});
		PrintArray ($OutputHandle,@{$OutputStructure{$Gene}->{$GuideSequence}->{'TSSs'}});
		PrintArray ($OutputHandle,@{$OutputStructure{$Gene}->{$GuideSequence}->{'Positions'}});
		print {$OutputHandle}  $OutputStructure{$Gene}->{$GuideSequence}->{'IsExtra'} . "\t";
		print {$OutputHandle}  $OutputStructure{$Gene}->{$GuideSequence}->{'IsUnique'};
		print {$OutputHandle}  "\n";
	}
}

close(IN) or die "ERROR in $0:Could not close inputfile\n";
close($OutputHandle) or die "ERROR in $0:Could not close outputfile\n";

sub PrintArray($@) {
	my $FileHandle=shift(@_);
	my @ArrayToPrint=@_;
	my $ArraySize=scalar @ArrayToPrint;
	my $ItemNumber=0;
	foreach my $Item(@ArrayToPrint) {
		$ItemNumber++;
		print {$FileHandle} $Item;
		unless($ItemNumber>=$ArraySize) {
			print {$FileHandle} ",";
		}
	}
	print {$FileHandle} "\t";
}