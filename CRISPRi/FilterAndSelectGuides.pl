use Getopt::Std;
use warnings;
use strict;

my %opts;
my %PotentialGuides;
my %GuidesPerGene;
my %SelectedGuides;
my $MaxGuidesPerTranscript=5;
getopt( 'i', \%opts );
open(IN, my $InputFile=$opts{'i'}) or die "ERROR in $0:Could not open inputfile\n";

#Read header
my $Header=<IN>;
chomp($Header);

#Loop over all guides in inputfile and place all potential guides in a hash
while (defined (my $Line=<IN>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
	my @Transcripts=split(/,/,$LineValues[6]);
	my @TSSs=split(/,/,$LineValues[7]);
	my @Positions=split(/,/,$LineValues[8]);
	$GuidesPerGene{$LineValues[0]}->{$LineValues[1]}->{'LineValues'}=\@LineValues;
	$GuidesPerGene{$LineValues[0]}->{$LineValues[1]}->{'Score'}=0;
	for (my $i=0;$i<scalar @Transcripts;$i++) {
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}=0;
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'TSS'}=$TSSs[$i];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Position'}=$Positions[$i];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Sense'}=$LineValues[2];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'GuideLength'}=$LineValues[3];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Chromosome'}=$LineValues[4];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Orientation'}=$LineValues[5];
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'IsExtra'}=$LineValues[9];
		#Build score:
		#Add 100 points if this guide targets this transcript close to the TSS
		#Add 10 points for every other transcript it hits close to the TSS
		#Add 1 point for every other transcript it hits
		foreach (@Positions) {
			if ($_ <=100 && $_ >=-25) {
				$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}=$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}+10;
			}
			else {
				$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}++;
			}
		}
		$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}=$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'Score'}+90 if ($Positions[$i] <=100 && $Positions[$i]>=-25);
		
		for (my $j=11;$j<=18;$j++) {
			push(@{$PotentialGuides{$LineValues[0]}->{$Transcripts[$i]}->{$LineValues[1]}->{'OffTargets'}},$LineValues[$j]);
		}
	}	
}

#Perform selection procedure
#For every transcript, sort guide list according to the scores built up in the previous loop. 
#Select max 10 and move on to next transcript
foreach my $Gene (keys %PotentialGuides) {
	foreach my $Transcript (keys $PotentialGuides{$Gene}) {
		my $IncludePromiscuousGuides=0;
		#Run the loop twice, one to search for guides without 1MM mismatches and then if necessary again to include those
		for(my $i=0;$i<2;$i++) {
			foreach my $GuideSequence (sort {$PotentialGuides{$Gene}->{$Transcript}->{$b}->{'Score'} <=> $PotentialGuides{$Gene}->{$Transcript}->{$a}->{'Score'}} keys $PotentialGuides{$Gene}->{$Transcript}) {
				if ($PotentialGuides{$Gene}->{$Transcript}->{$GuideSequence}->{'GuideLength'}==20) {
					if ($IncludePromiscuousGuides || $PotentialGuides{$Gene}->{$Transcript}->{$GuideSequence}->{'OffTargets'}->[0]==0) {
						#Push this particular guide into the selected guides of all transcripts it targets close, unless it's already there
						my @TargetedTranscripts=split(/,/,$GuidesPerGene{$Gene}->{$GuideSequence}->{'LineValues'}->[6]);
						my @TargetedPositions=split(/,/,$GuidesPerGene{$Gene}->{$GuideSequence}->{'LineValues'}->[8]);
						foreach (@TargetedTranscripts) {
							#This line adds the guide to the selected guides of the current transcript or to the selected guides of any other transcripts it targets close (i.e. of which the score is >=100)
							#it also only adds it if it's not already there
							if ($_ eq $Transcript || $PotentialGuides{$Gene}->{$_}->{$GuideSequence}->{'Score'} >=100) {
								push(@{$SelectedGuides{$Gene}->{$_}->{'Guides'}},$GuideSequence) unless ($GuideSequence ~~ @{$SelectedGuides{$Gene}->{$_}->{'Guides'}});
							}	
						}
						#Stop selecting guides if max is reached
						last if (@{$SelectedGuides{$Gene}->{$Transcript}->{'Guides'}}>=$MaxGuidesPerTranscript);
					}
				}
			}
			if (@{$SelectedGuides{$Gene}->{$Transcript}->{'Guides'}}<$MaxGuidesPerTranscript) {
				$IncludePromiscuousGuides=1;
			} 
			else {
				#In case we already reached max guides for this transcript, make sure it doesn't loop into finding the same ones again
				$i++;
			}
		}
	}
}

#TODO: Subselect guides in case any gene has more than a certain number of guides selected
my %OutputGuides;
foreach my $Gene (keys %SelectedGuides) {
	foreach my $Transcript (keys $SelectedGuides{$Gene}) {
		foreach (@{$SelectedGuides{$Gene}->{$Transcript}->{'Guides'}}) {
			$OutputGuides{$Gene}->{$_}++;
		}
	}
}

open (OUT,">",$InputFile . ".selected") or die "ERROR in $0:Could not open for output file " . $InputFile . ".selected\n";
print OUT $Header . "\n";
foreach my $Gene (sort {$a cmp $b} keys %OutputGuides) {
	foreach my $Guide (keys $OutputGuides{$Gene}) {
		foreach (@{$GuidesPerGene{$Gene}->{$Guide}->{'LineValues'}}) {
			print OUT $_ . "\t";
		} 
		print OUT "\n";
	}
}
foreach my $Gene (sort {$a cmp $b} keys %OutputGuides) {
	open(FASTQ,">","fastq/" . $Gene . ".fastq");
	my $Counter=0;
        foreach my $Guide (keys $OutputGuides{$Gene}) {
		$Counter++;
                my @values=@{$GuidesPerGene{$Gene}->{$Guide}->{'LineValues'}};
		print FASTQ ">" . $Gene . "_" . $Counter . "\n";
		print FASTQ $values[1] . "NGG\n";
        }
	close(FASTQ);
}

close (OUT) or die "ERROR in $0:Could not close output file\n";
close(IN) or die "ERROR in $0:Could not close inputfile\n";
