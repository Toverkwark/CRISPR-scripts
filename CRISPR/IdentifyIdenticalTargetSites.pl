use Getopt::Std;
use warnings;
use strict;

my$ScriptName="IdentifyIdenticalTargetSites.pl";
my %opts;
my %Protospacers;
getopt( 'io', \%opts );
die "ERROR in $ScriptName: No inputfile given.\n" unless my $InputFile = $opts{'i'};
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
open (IN, $InputFile) or die "ERROR in $ScriptName: Cannot open inputfile $InputFile\n";
open (OUT, ">", $OutputFile) or die "ERROR in $ScriptName: Cannot open outpufile $OutputFile\n";
while (defined(my $Line = <IN>)) {
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
		$Protospacers{$TargetSequence}->[1]=$Line;
	}
}
foreach my $TargetSequence (keys %Protospacers) {
	#Only output protospacers with unique target sites 
	if($Protospacers{$TargetSequence}->[0]==1) {
		print OUT $Protospacers{$TargetSequence}->[1] . "\n";
	}
}
close (IN) or die "ERROR in $ScriptName: Cannot close inputfile $InputFile\n";
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";