use Getopt::Std;
use warnings;
use strict;
use LocalSettings;
my %LocalSettings=getconfig();
my $Bowtie=$LocalSettings{'Bowtie'};
my $IndexedHumanGenome=$LocalSettings{'IndexedHumanGenome'};
my $NumberOfCoresToUse=$LocalSettings{'NumberOfCoresToUse'};
my %opts;
my %InputGuides;
my %Guides;
getopt( 'io', \%opts );
print ("Usage:perl $0 -i InputFile [-o OutputFile (default=InputFile.unique])\n");
die "ERROR in $0: No inputfile given.\n" unless my $InputFile = $opts{'i'};
my $OutputFile;
$OutputFile=$InputFile . ".unique" unless $OutputFile = $opts{'o'};

#Parse the input file to read all possible Guide sequences
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
my $Header=<IN>; #Skip header line
chomp($Header);
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @LineValues = split( /\t/, $Line );
	#Add any guide sequence to a hash. In this way, if there are multiple guideRNAs with the same sequence, targeting perhaps multiple RefSeqs of the same gene, only one sequence will be evaluated for uniqueness because in the next paragraph, only one line per unique guide will be written1
	$InputGuides{$LineValues[6]}++;
}
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";

#Write all possible target sequences to evaluate to a temporary output file. 
open (TEMPOUTPUT, ">", ($InputFile . ".tmp")) or die "ERROR in $0: Cannot open tmp file $InputFile.tmp\n";
foreach my $InputGuide (keys %InputGuides) {
	print TEMPOUTPUT $InputGuide . "AGG\n";
	print TEMPOUTPUT $InputGuide . "TGG\n";
	print TEMPOUTPUT $InputGuide . "CGG\n";
	print TEMPOUTPUT $InputGuide . "GGG\n";
}
close (TEMPOUTPUT) or die "ERROR in $0: Cannot close tmp file $InputFile.tmp\n";


#Search the genome for occurrences of these guides
`$Bowtie $IndexedHumanGenome/hg19 -r $InputFile.tmp -t -S $InputFile.tmp.matched --no-hd --score-min L,-5,0 -k 2 -p $NumberOfCoresToUse`;

#Parse the resulting file and count occurrences for every guides
open (IN, $InputFile . ".tmp.matched") or die "ERROR in $0: Cannot open temp file for guide parsing\n";
while (defined(my $Line=<IN>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
	my $TargetSequence = $LineValues[9];
	if($LineValues[1] & 16) {
		$TargetSequence =~ tr/ACTG/TGAC/;
		$TargetSequence = reverse($TargetSequence);
	}
	my $GuideSequence=substr($TargetSequence,0,length($TargetSequence)-3);
	if(!($LineValues[1] & 4)) {
		$Guides{$GuideSequence}->[0]++;
	}
}
close (IN) or die "ERROR in $0: Cannot close temp file\n";

#Parse the input file and add a column indicating the # of genomic occurrences
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
open (OUT, ">", $OutputFile) or die "ERROR in $0: Cannot open outpufile $OutputFile\n";
<IN>; #Skip header line
print OUT $Header . "\tUnique\n";
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @LineValues = split( /\t/, $Line );
	#Only output guides for which targets are unique in the genome
	if($Guides{$LineValues[6]}->[0] == 1) {
		print OUT $Line . "\t1\n";
	}
}
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";
unlink($InputFile . ".tmp");
unlink($InputFile . ".tmp.matched");

