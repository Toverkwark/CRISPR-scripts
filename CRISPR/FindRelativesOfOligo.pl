use Getopt::Std;
use warnings;
use strict;
my $ScriptName="FindRelativesOfOligo.pl";
my %opts;
my %Relatives;

getopt( 'oids', \%opts );
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
die "ERROR in $ScriptName: No Input sequence given.\n" unless my $InputSequence = $opts{'i'};
my $MaxDepth = $opts{'d'};
die "ERROR in $ScriptName: No 3' seed length given.\n" unless my $SeedLength = $opts{'s'};
open (OUT, ">>", $OutputFile) or die "Cannot open outputfile $OutputFile\n";

#Create all possible distance relatives
my $SeedSequence=substr($InputSequence,20-$SeedLength,$SeedLength);
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
	print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "GGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "AGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "TGG\n";
	print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "CGG\n";
	#print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "GAG\n";
	#print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "AAG\n";
	#print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "TAG\n";
	#print OUT ">MM"  . $Depth . "_Relative_$InputSequence\n$Relative" . "CAG\n";
}
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";
