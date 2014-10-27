use Getopt::Long;
use List::Util qw(sum min max);
use Statistics::Basic qw(:all);
use Term::ANSIColor;
use strict;

print "Usage:perl $0 -input -output\n-input\tName of input file\n-output\tName of output file. Default is inputfile.svg\n";

GetOptions(
	"input=s"  => \my $InputFile,
	"output=s" => \my $OutputFile,
);

if ( !$OutputFile ) {
	$OutputFile = $InputFile . ".svg";
}

#Write header
`cat header >$OutputFile`;


#Read in information from inputfile
print "Reading information from inputfile\n";
open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
my $Line=<INPUT>;
my @Screen;
while (defined (my $Line=<INPUT>)) {
	chomp($Line);
	my @LineValues=split(/\t/,$Line);
	push(@{$Screen[0]},$LineValues[0]);
	push(@{$Screen[1]},substr($LineValues[0],0,index($LineValues[0],'-')));
	for (my $i=2;$i<=13;$i++) {
		push(@{$Screen[$i]},$LineValues[$i]);
	}
}
close(INPUT) or die "ERROR in $0:Could not close inputfile $InputFile\n";

#Add 1 to all read counts to prevent problems with taking the log
print "Adding 1 to all read counts to preven problems with taking the log\n";
for my $Tag (0..@{$Screen[0]}-1) {
	for my $Column (2..13) {
		$Screen[$Column][$Tag]++;	
	}
}
	
#Normalize data
print "Normalizing data\n";
my @SizeFactorData;
my @SizeFactors;
my @NormalizedData;
#First, gather sizefactor data
for my $Tag (0..@{$Screen[0]}-1) {
	my @LogReads;
	for my $Column (2..13) {
		push(@LogReads,log($Screen[$Column][$Tag])/log(2));
	}
	my $MeanLogReads=mean(@LogReads);
	for my $Column (2..13) {
		$SizeFactorData[$Column-2][$Tag]=$LogReads[$Column-2]-$MeanLogReads;
	}
}

#Calculate size factors
for my $Column (0..11) {
	push(@SizeFactors,median(@{$SizeFactorData[$Column]}));
}

#Apply the normalization
for my $Column (2..13) {
	for my $Tag (0..@{$Screen[0]}-1) {
		$NormalizedData[$Column-2][$Tag]=2**((log($Screen[$Column][$Tag])/log(2))-$SizeFactors[$Column-2]);
	}	
}

my $ControlColumn=5;
my $TreatedColumn=9;
my $HorizontalPixels=333;
my $VerticalPixels=250;
my $PixelSize=1;
my $MinimalX=0;
my $MaximalX=10;
my $MinimalY=-5;
my $MaximalY=5;
my @XPositions;
my @YPositions;
#Create plot
open( OUTPUT, ">>", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";
for my $Tag (0..@{$NormalizedData[0]}-1) {
	my $XPosition=0.5*(log($NormalizedData[$ControlColumn][$Tag])/log(2)+log($NormalizedData[$TreatedColumn][$Tag])/log(2));
	my $YPosition=(log($NormalizedData[$TreatedColumn][$Tag])/log(2))-(log($NormalizedData[$ControlColumn][$Tag])/log(2));
	
	$XPosition=$XPosition-$MinimalX;
	$XPosition=$XPosition*$HorizontalPixels/($MaximalX-$MinimalX);
	$YPosition=$YPosition-$MinimalY;
	$YPosition=$YPosition*$VerticalPixels/($MaximalY-$MinimalY);
	push(@XPositions,$XPosition);
	push(@YPositions,$YPosition);
	my $Line='<circle class="' . $Screen[1][$Tag] . '" selected="no" cx="' . $XPosition . '" cy="' . $YPosition . '" r="' . $PixelSize . '" stroke="black" stroke-width=".1" fill="red" onclick="ClickDetected(evt)" onmouseover="MouseOverDetected(evt)" onmouseout="MouseOutDetected(evt)"/>';
	print OUTPUT $Line . "\n";
}

print "Minimum X:" . (min @XPositions) . "\tMaximum X:" . (max @XPositions) . "\tMinimum Y:" . (min @YPositions) . "\tMaximum Y:" . (max @YPositions) . "\n";


#Write footer
print OUTPUT "</svg>\n";


close(OUTPUT) or die "ERROR in $0:Could not close outputfile $OutputFile\n";

sub mean {
    return sum(@_)/@_;
}