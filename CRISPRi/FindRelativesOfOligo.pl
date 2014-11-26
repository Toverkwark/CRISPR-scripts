use Getopt::Std;
use warnings;
use strict;
my %opts;
my %Relatives;

getopt( 'oid', \%opts );
die "ERROR in $0: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
die "ERROR in $0: No Input sequence given.\n" unless my $OriginalSeedSequence = $opts{'i'};
#Cap the search for off-target effects to the 20 most PAM proximal nt
my $SeedSequence=substr($OriginalSeedSequence,-20);
my $MaxDepth = $opts{'d'};
open (OUT, ">>", $OutputFile) or die "Cannot open outputfile $OutputFile\n";


#Create all possible distance relatives
$Relatives{0}->{$SeedSequence}++;
for (my $Depth=1;$Depth<=$MaxDepth;$Depth++)
{
	foreach my $MutateSequence (keys $Relatives{$Depth-1}) {
		for (my $i=0;$i<length($MutateSequence);$i++) {
			my $LeftSequence = substr($MutateSequence,0,$i);
			my $RightSequence = substr($MutateSequence,$i+1);
			for (my $j=0;$j<=3;$j++) {
				my $NewSequence = $LeftSequence . substr("ACTG",$j,1) . $RightSequence;
				my $SequenceFoundAlready = 0;
				for (my $InvestigatingDepth=$Depth-1;$InvestigatingDepth>=0;$InvestigatingDepth--)
				{
					if($Relatives{$InvestigatingDepth}->{$NewSequence}) {
						$SequenceFoundAlready = 1;
						last;
					}
				}
				$Relatives{$Depth}->{$NewSequence}++ unless $SequenceFoundAlready;
			}
		}
	}  	
}
my $Depth=$MaxDepth;
foreach my $Relative (keys $Relatives{$Depth}) {
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "GGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "AGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "TGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "CGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "GAG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "AAG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "TAG\n";
	print OUT ">MM"  . $Depth . "_Relative_$OriginalSeedSequence\n$Relative" . "CAG\n";
}
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";
