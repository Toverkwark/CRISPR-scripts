use Getopt::Long;
use List::Util qw(sum min max);
use Statistics::Basic qw(:all);
use Term::ANSIColor;
use strict;

my $ControlColumn=1;
my $TreatedColumn=4;
my $HorizontalPixels=666;
my $VerticalPixels=500;
my $PixelSize=1;
my $MinimalX=0;
my $MaximalX=20;
my $MinimalY=-10;
my $MaximalY=10;

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
my $Row=-1;
while (defined (my $Line=<INPUT>)) {
	$Row++;
	chomp($Line);
	my @LineValues=split(/\t/,$Line);	
	push(@{$Screen[$Row]},substr($LineValues[0],0,index($LineValues[0],'-')));
	push(@{$Screen[$Row]},@LineValues);
}
close(INPUT) or die "ERROR in $0:Could not close inputfile $InputFile\n";

#Add 1 to all read counts to prevent problems with taking the log
print "Adding 1 to all read counts to preven problems with taking the log\n";
for my $Row (0..@Screen-1) {
	for my $Column (3..@{$Screen[0]}-1) {
		$Screen[$Row][$Column]++;	
	}
}

#First, gather sizefactor data
print "Gathering size factor data\n";
my @SizeFactorData;
for my $Row (0..@Screen-1) {
	my @LogReads;
	for my $Column (3..@{$Screen[0]}-1) {
		push(@LogReads,log($Screen[$Row][$Column])/log(2));
	}
	my $MeanLogReads=mean(@LogReads);
	for my $Column (0..@LogReads-1) {
		$SizeFactorData[$Row][$Column]=$LogReads[$Column]-$MeanLogReads;
	}
}

#Calculate size factors
print "Calculating size factors\n";
my @SizeFactors;
for my $Column (0..@{$SizeFactorData[0]}-1) {
	push(@SizeFactors,median(map {$$_[$Column]} @SizeFactorData));
}

#Apply the normalization
print "Applying Normalization\n";
my @NormalizedData;
for my $Row (0..@Screen-1) {
	for my $Column (0..@SizeFactors-1) {
		$NormalizedData[$Row][$Column]=2**((log($Screen[$Row][$Column+3])/log(2))-$SizeFactors[$Column]);
	}	
}

#Create plot, first create axes
print "Creating plot\n";
open( OUTPUT, ">>", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";
print OUTPUT '<line x1="0" y1="0" x2="' . $HorizontalPixels . '" y2="0" stroke-width="0.5" stroke="black" />ń';
print OUTPUT '<line x1="0" y1="0" x2="0" y2="' . $VerticalPixels . '" stroke-width="0.5" stroke="black" />ń';
print OUTPUT '<line x1="' . $HorizontalPixels . '" y1="0" x2="' . $HorizontalPixels . '" y2="' . $VerticalPixels . '" stroke-width="0.5" stroke="black" />ń';
print OUTPUT '<line x1="0" y1="' . $VerticalPixels . '" x2="' . $HorizontalPixels . '" y2="' . $VerticalPixels . '" stroke-width="0.5" stroke="black" />ń';
for (my $i=$MinimalY;$i<=$MaximalY;$i++) {
	my $VerticalLevel=(-$i-$MinimalY)*$VerticalPixels/($MaximalY-$MinimalY);
	print OUTPUT '<line x1="0" y1="' . $VerticalLevel . '" x2="' . $HorizontalPixels . '" y2="' . $VerticalLevel .'" stroke-width=".2" stroke="black" />ń';
	print OUTPUT '<text font-family="Arial" font-size="10" x="' . 0.01*$HorizontalPixels . '" y="' . $VerticalLevel . '">' . $i . '</text>\n"';
}

#Create plot, write all the datapoints
my @jet=@{ReadColorMap("jet")};
for my $Row (0..@NormalizedData-1) {
	#Determine MA Values
	my $A=0.5*(log($NormalizedData[$Row][$ControlColumn])/log(2)+log($NormalizedData[$Row][$TreatedColumn])/log(2));
	my $M=-((log($NormalizedData[$Row][$TreatedColumn])/log(2))-(log($NormalizedData[$Row][$ControlColumn])/log(2)));
	#Normalize to current min/max settings
	my $XPosition=($A-$MinimalX)*$HorizontalPixels/($MaximalX-$MinimalX);
	my $YPosition=($M-$MinimalY)*$VerticalPixels/($MaximalY-$MinimalY);
	$YPosition=0 if $YPosition<0;
	$YPosition=$VerticalPixels if $YPosition > $VerticalPixels;
	$XPosition=0 if $XPosition<0;
	$XPosition=$HorizontalPixels if $XPosition > $HorizontalPixels;
	#Map color
	my $ColorPosition=($M-0.25*$MinimalY)*$VerticalPixels/(0.25*$MaximalY-0.25*$MinimalY);
	$ColorPosition=0 if $ColorPosition<0;
	$ColorPosition=$VerticalPixels if $ColorPosition > $VerticalPixels;
	my $ColorNumber=int(63*$ColorPosition/$VerticalPixels);
	my @FillColor=@{$jet[$ColorNumber]};
	my $Line='<circle class="' . ($Screen[$Row][0]) . '" selected="no" cx="' . $XPosition . '" cy="' . $YPosition . '" r="' . $PixelSize . '" stroke="black" stroke-width=".1" fill="rgb(' . $FillColor[0] . ',' .$FillColor[1] . ',' . $FillColor[2] . ')" onclick="ClickDetected(evt)" onmouseover="MouseOverDetected(evt)" onmouseout="MouseOutDetected(evt)"/>';
	print OUTPUT $Line . "\n";
}

#Write footer
print OUTPUT "</svg>\n";
close(OUTPUT) or die "ERROR in $0:Could not close outputfile $OutputFile\n";

sub mean {
    return sum(@_)/@_;
}

sub MapColor ($$$) {
	my ($YPosition, $MinimalY, $MaximalY) = @_;
	$YPosition
}

sub ReadColorMap($) {
	my @Colormap;
	my $FileName=shift(@_);
	my $Row=-1;
	open(IN,$FileName) or die "Could not open colormap file $FileName\n";
	while(defined(my $Line=<IN>)) {
		$Row++;
		chomp($Line);
		my @Color=split('\t',$Line);
		$Color[0]=int($Color[0]*255);
		$Color[1]=int($Color[1]*255);
		$Color[2]=int($Color[2]*255);	
		push(@{$Colormap[$Row]},@Color);
	}
	return \@Colormap;
}
