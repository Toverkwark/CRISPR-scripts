use Getopt::Std;
use warnings;
use strict;
use Cwd 'abs_path';
use File::Basename;

#This script creates an svg graphical representation of the coding sequence of a gene with CRISPR location mapped onto it
#As input arguments, it takes the refseq id (r) and the species (s) for which to find the info
#The whole image is created such that it is 1500 pixels wide. This can best be adjusted by using the viewBox argument in the svg header in the header file
my $DirName = dirname(abs_path($0));
my $HeaderFile = $DirName . "/header";
my $FooterFile = $DirName . "/footer";
my ($Genome, $RefSeqFile, $QualitiesFileLocation);

#First, set some parameters that determine how the svg file will look
my $YOffset=0;
my $AnimationDistance = 30;
my $TriangleHeight = 10;
my $TriangleWidth = 4;
my $CollisionWidth = 8;
my $LevelOffset = 10;

#Get options given to the script
my %ScriptOptions;
getopt('irs', \%ScriptOptions);
print "Usage:$0 -irs\n-i File containing the sites to be displayed\n-r RefSeq ID\n-s Species. Default human\n" if !%ScriptOptions;
die "ERROR in $0: No input file given.\n" unless my $InputFile = $ScriptOptions{'i'};
die "ERROR in $0: No RefSeq ID given.\n" unless my $RefSeqID = $ScriptOptions{'r'};
my $Species;
$Species = 'human' unless $Species = $ScriptOptions{'s'};

#Assign proper genome id
$Species = lc($Species);
if ($Species eq 'mouse') {
	$Genome = "mm10";
	$RefSeqFile = $DirName . "../refseq/mm10.txt";
}
else {
	if ($Species eq 'human') {
		$Genome = "hg19";
		$RefSeqFile = "../refseq/hg19.txt";
	}
	else {
		die "ERROR in $0: Species $Species is not currently handled by the iKRUNC scripts\n";
	}
}

#Open inputfile and svg file
open (IN, $InputFile) or die "ERROR in $0: Cannot open inputfile $InputFile\n";
open (OUT, ">", $RefSeqID . ".svg") or die "ERROR in $0: Cannot create svg file\n";
open (OUTHTML, ">", $RefSeqID . ".svg.html") or die "ERROR in $0: Cannot create svg html file\n";


#First read in all intron/exon information of the gene we're plotting
my $RefSeqInfo = `grep -P "$RefSeqID\t" $RefSeqFile`;
die "ERROR in $0: RefSeq $RefSeqID cannot be found in the database.\n" if !$RefSeqInfo;
my @RefSeqValues = split( /\t/, $RefSeqInfo );
my $Chromosome = $RefSeqValues[2];
my $GeneStart    = $RefSeqValues[4];
my $GeneEnd      = $RefSeqValues[5];
my $ProteinStart = $RefSeqValues[6];
my $ProteinEnd   = $RefSeqValues[7];
my $GeneOrientation = $RefSeqValues[3];
my $NumberOfExons  = $RefSeqValues[8];
my @ExonStartSites = split( /,/, $RefSeqValues[9] );
my @ExonEndSites   = split( /,/, $RefSeqValues[10] );

#Determine the size of the mRNA
my $mRNASize = 0;
for (my $i=0;$i<$NumberOfExons;$i++) {
	$mRNASize =  $mRNASize + $ExonEndSites[$i] - $ExonStartSites[$i];
}

#Read in all CRISPRS to be mapped
my %DisplayObjects;
while (defined(my $Line = <IN>)) {
	chomp($Line);
	my @TargetSites = split( /\t/, $Line );
	my $TargetChromosome = $TargetSites[2];
	my $TargetCutSite = $TargetSites[3];
	my $TargetOrientation = $TargetSites[4];
	
	#Verify that the target is in the gene, or within 250nt of the TSS
	if ($TargetChromosome eq $Chromosome && ($TargetCutSite >= $GeneStart-250 && $TargetCutSite <= $GeneEnd + 250)) {
		#Verify that the target is targeting the intended RefSeqID
		my @ListOfTargetRefSeqIDs=split(/,/,$TargetSites[10]);
		foreach my $TargetRefSeqID (@ListOfTargetRefSeqIDs) {
			print "Investigating $TargetRefSeqID\t";
			if($RefSeqID eq $TargetRefSeqID) {
				@{$DisplayObjects{$TargetOrientation}->{$TargetCutSite}}=@TargetSites;
			}
		}
	}
}
close (IN);

#Determine coding sequence position of cutsites
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
		my $CodingSequencePosition = 0;
		if($GeneOrientation eq '+') {
			for (my $i=0;$i<$NumberOfExons;$i++) {
				if ($TargetCutSite <= $ExonEndSites[$i]) {
					if((($TargetCutSite-$ExonStartSites[$i])>0) || $i==0) {
						$CodingSequencePosition = $CodingSequencePosition + ($TargetCutSite-$ExonStartSites[$i]);
						last;
					}
				} 
				else
				{
					$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i] - $ExonStartSites[$i]);
				}
			}
		}
		else {
			for (my $i=$NumberOfExons-1;$i>=0;$i--) {
				if ($TargetCutSite >= $ExonStartSites[$i]) {
					if((($ExonEndSites[$i]-$TargetCutSite)>0) || $i==($NumberOfExons-1)) {
						$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i]-$TargetCutSite);
						last;
					}
				} 
				else
				{
					$CodingSequencePosition = $CodingSequencePosition + ($ExonEndSites[$i] - $ExonStartSites[$i]);
				}
			}
		}
		my $RelativeMarkerPosition=1400*$CodingSequencePosition/$mRNASize;
		$RelativeMarkerPosition=$RelativeMarkerPosition+100;	
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[11]=$RelativeMarkerPosition;
	}
}

#Determine and assign colors
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) { 
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[12]="0,0,255";
	}
}
		
#Do collision detection and resolution
my $MaxCollisionLevelSenseStrand=1;
my $MaxCollisionLevelAntisenseStrand=1;
foreach my $Orientation (keys %DisplayObjects) {
	my $CollisionsDetected=1;
	my $CollisionLevel=1;
	my %AllPositions;
	while ($CollisionsDetected) {
		$CollisionsDetected = 0;
		foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
			my $CollisionDetectedForThisTarget = 0;
			my $RelativeMarkerPosition=$DisplayObjects{$Orientation}->{$TargetCutSite}->[11];
			if(!($DisplayObjects{$Orientation}->{$TargetCutSite}->[13])) {
				if($AllPositions{$CollisionLevel}) {
					foreach my $Position (keys $AllPositions{$CollisionLevel}) {
						if ($RelativeMarkerPosition <= ($Position+$CollisionWidth) && $RelativeMarkerPosition >= ($Position-$CollisionWidth)) {
						#	print "Detected a collision with $RelativeMarkerPosition at position $Position at level $CollisionLevel\n";
							$CollisionsDetected++;
							$CollisionDetectedForThisTarget++;
							last;
						}
					}
				}
				if (!$CollisionDetectedForThisTarget) {
					$AllPositions{$CollisionLevel}->{$RelativeMarkerPosition}++;
					$DisplayObjects{$Orientation}->{$TargetCutSite}->[13]=$CollisionLevel;
				}
			}
		}
		$CollisionLevel++;
	}
	if($CollisionLevel > $MaxCollisionLevelSenseStrand && $Orientation eq '+') {
		$MaxCollisionLevelSenseStrand=$CollisionLevel-1;
	}
	if($CollisionLevel > $MaxCollisionLevelAntisenseStrand && $Orientation eq '-') {
		$MaxCollisionLevelAntisenseStrand=$CollisionLevel-1;
	}
}
$YOffset=$YOffset+$AnimationDistance+$TriangleHeight+($MaxCollisionLevelSenseStrand-1)*$LevelOffset;

#Now print rectangles for all exons
my $CurrentExon=0;
#Start 100 away from the beginning to accomodate for CRISPRs that are 5' of the TSS
my $CurrentPosition=0;
my $BoxSize;
$CurrentPosition=1500 if ($GeneOrientation eq '-');
#Start by drawing a piece of DNA upstream from the TSS
$BoxSize =  100;
my $SVGFile = "  <rect x=\"" . 0 . "\" y=\"" . (11+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"8\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" />\'\n";
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');

#Until you meet the proteinstart, plot full exons as UTR
while(!($ProteinStart >= $ExonStartSites[$CurrentExon] && $ProteinStart <= $ExonEndSites[$CurrentExon])) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	$SVGFile = $SVGFile . "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (7+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
}
#Print the piece of exon until the proteinstart as UTR
$BoxSize =  1400*($ProteinStart - $ExonStartSites[$CurrentExon])/$mRNASize;
$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
$SVGFile = $SVGFile .  "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (7+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
my $UseAlternativeStart=1;
#Until you meet the proteinend, plot exons as exon
while(!($ProteinEnd >= $ExonStartSites[$CurrentExon] && $ProteinEnd <= $ExonEndSites[$CurrentExon])) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ProteinStart)/$mRNASize if $UseAlternativeStart;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	$SVGFile = $SVGFile . "  <rect x=\"" . $CurrentPosition . "\" y=\"" . ($YOffset) . "\" width=\"" . $BoxSize . "\" height=\"30\" fill=\"url(#grad1)\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
	$UseAlternativeStart=0;
}
#Print the piece of exon until the proteinend as exon
$BoxSize =  1400*($ProteinEnd - $ExonStartSites[$CurrentExon])/$mRNASize;
$BoxSize =  1400*($ProteinEnd - $ProteinStart)/$mRNASize if $UseAlternativeStart;
$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
$SVGFile = $SVGFile . "  <rect x=\"" . $CurrentPosition . "\" y=\"" . ($YOffset) . "\" width=\"" . $BoxSize . "\" height=\"30\" fill=\"url(#grad1)\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
$UseAlternativeStart=1;	
#Print the rest as UTR
while($CurrentExon<$NumberOfExons) {
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ExonStartSites[$CurrentExon])/$mRNASize;
	$BoxSize =  1400*($ExonEndSites[$CurrentExon] - $ProteinEnd)/$mRNASize if $UseAlternativeStart;
	$CurrentPosition=$CurrentPosition-$BoxSize if ($GeneOrientation eq '-');
	$SVGFile = $SVGFile . "  <rect x=\"" . $CurrentPosition . "\" y=\"" . (7+$YOffset) . "\" width=\"" . $BoxSize . "\" height=\"16\" fill=\"black\" stroke=\"black\" stroke-width=\"1\" onmousemove=\"ShowTooltip(evt, \'Exon " . ($GeneOrientation eq '+' ? ($CurrentExon+1) : $NumberOfExons-$CurrentExon) . "\')\" onmouseout=\"HideTooltip(evt)\"/>\'\n";	
	$CurrentPosition = $CurrentPosition + $BoxSize if ($GeneOrientation eq '+');
	$CurrentExon++;
	$UseAlternativeStart=0;
}

#Perform printing of cut sites to output file
my $ExonOffset=30;
my $TargetID=0;
my %TableObjects;
foreach my $Orientation (keys %DisplayObjects) {
	#Loop through all cut sites sorted by score
	foreach my $TargetCutSite (sort {$DisplayObjects{$Orientation}->{$a}->[0] cmp $DisplayObjects{$Orientation}->{$b}->[0]} keys $DisplayObjects{$Orientation}) {		
		my @TargetSites = @{$DisplayObjects{$Orientation}->{$TargetCutSite}};		
		my $TargetChromosome = $TargetSites[2];
		my $TargetCutSite = $TargetSites[3];
		my $TargetOrientation = $TargetSites[4];
		my $TargetSequence=$TargetSites[1];
		my $TargetNumberOfIdentical3PrimeTargets=$TargetSites[5];
		my $TargetNumberOfIdentical3PrimeTargetsNearExons=$TargetSites[6];
		my $TargetDegree=$TargetSites[7];
		my $TargetClosestRelatives=$TargetSites[8];
		my $TargetClosestRelativesNearExons=$TargetSites[9];
		$TargetID=$TargetID+1;		
		my $RelativeMarkerPosition=$TargetSites[11];
		my $TriangleColor=$TargetSites[12];
		my $LayerOffset = $LevelOffset*(($TargetSites[13]) - 1);
		my $TargetLabel = $TargetSequence;
		if($Orientation eq '+') {
			$SVGFile = $SVGFile . "<polygon id=\"" . $TargetID . "\" class=\"Triangle\" points=\"" . ($RelativeMarkerPosition - $TriangleWidth) . "\," . ($YOffset-$LayerOffset-$TriangleHeight) . " " . $RelativeMarkerPosition . "\," . ($YOffset-$LayerOffset) . " " . ($RelativeMarkerPosition+$TriangleWidth) . "\," . ($YOffset-$LayerOffset-$TriangleHeight) . "\" style=\"fill:rgb(" . $TriangleColor . ");stroke:black;stroke-width:1\" onmousemove=\"ShowTooltip(evt, \'" . $TargetLabel . "\')\" onmouseout=\"HideTooltip(evt)\" onclick=\"ClickTriangle(" . $TargetID . ")\">";
			my $RandomTime=1*rand();
			$SVGFile = $SVGFile . "<animateTransform attributeName=\"transform\" attributeType=\"XML\" type=\"translate\" from=\"0 -" . $AnimationDistance . "\" to=\"0 0\" dur=\"" . $RandomTime . "s\"/>";
			$SVGFile = $SVGFile . "</polygon>\n";
		}
		else {
			$SVGFile = $SVGFile . "<polygon id=\"" . $TargetID . "\"  class=\"Triangle\" points=\"" . ($RelativeMarkerPosition - $TriangleWidth) . "\," . ($YOffset+$ExonOffset+$LayerOffset+$TriangleHeight) . " " . $RelativeMarkerPosition . "\," . ($YOffset+$ExonOffset+$LayerOffset) . " " . ($RelativeMarkerPosition+$TriangleWidth) . "\," . ($YOffset+$ExonOffset+$TriangleHeight+$LayerOffset) . "\" style=\"fill:rgb(" . $TriangleColor . ");stroke:black;stroke-width:1\" onmousemove=\"ShowTooltip(evt, \'" . $TargetLabel . "\')\" onmouseout=\"HideTooltip(evt)\" onclick=\"ClickTriangle(" . $TargetID . ")\">";
			my $RandomTime=1*rand();
			$SVGFile = $SVGFile . "<animateTransform attributeName=\"transform\" attributeType=\"XML\" type=\"translate\" from=\"0 " . $AnimationDistance . "\" to=\"0 0\" dur=\"" . $RandomTime . "s\"/>";
			$SVGFile = $SVGFile . "</polygon>\n";
		}	
		$TableObjects{$TargetID}->[1]=$TargetOrientation;
		$TableObjects{$TargetID}->[2]=$TargetCutSite;
		$TableObjects{$TargetID}->[3]=substr($TargetSequence,0,20) . "<b>" . substr($TargetSequence,20,3) . "</b>";
		$TableObjects{$TargetID}->[4]=$TargetNumberOfIdentical3PrimeTargets;
		$TableObjects{$TargetID}->[5]=$TargetNumberOfIdentical3PrimeTargetsNearExons;
		$TableObjects{$TargetID}->[6]=$TargetDegree;
		$TableObjects{$TargetID}->[7]=$TargetClosestRelatives;
		$TableObjects{$TargetID}->[8]=$TargetClosestRelativesNearExons;
	}
} 

#Calculate the total height of the figure:
my $TotalHeight=$YOffset+$ExonOffset+($LevelOffset*($MaxCollisionLevelAntisenseStrand-1))+$TriangleHeight+$AnimationDistance;

#Start writing the svg file
#First, write the header
print OUT '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" onload="init(evt)" width="100%" height="100%" viewBox="0 0 1500 ' . $TotalHeight . '" id="svgimage">\n';
print OUT `cat "$HeaderFile"`;

#Write the image file
print OUT $SVGFile;

#Print the footer
print OUT `cat "$FooterFile"`;

#Print the table file
print OUTHTML "<!DOCTYPE html>\n<html lang='en'>\n\t<head>\n\t\t<meta charset='utf-8'>\n\t\t<link rel='stylesheet' href='../style.css' type='text/css'></link>\n\t</head>\n";
print OUTHTML "<body>\n<table id='svgtable' class='svgtable' width='100%'>\n";
print OUTHTML "\t<th align='left'>Chromosome</th>\n";
print OUTHTML "\t<th align='left'>Orientation</th>\n";
print OUTHTML "\t<th align='left'>Position</th>\n";
print OUTHTML "\t<th align='left'>Sequence</th>\n";
print OUTHTML "\t<th align='left'>Identical 3'12nt sites</th>\n";
print OUTHTML "\t<th align='left'>Identical 3'12nt sites near exons</th>\n";
print OUTHTML "\t<th align='left'>Off-target #mismatches</th>\n";
print OUTHTML "\t<th align='left'># sites</th>\n";
print OUTHTML "\t<th align='left'># sites near exons</th>\n";
my $Rank=0;
foreach my $TableRow (sort {$TableObjects{$a}->[2] <=> $TableObjects{$b}->[2]} keys %TableObjects) {
	$Rank=$Rank+1;
	print OUTHTML "\t<tr id='" . $TableRow . ".table' onclick=parent.ClickTableRow('" . $TableRow . "')>\n";
	print OUTHTML "\t\t<td>" . $Chromosome . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[1] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[2] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[3] . "</td>\n";
	print OUTHTML "\t\t<td>" . ($TableObjects{$TableRow}->[4] -1) . "</td>\n";
	print OUTHTML "\t\t<td>" . ($TableObjects{$TableRow}->[5] -1) . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[6] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[7] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[8] . "</td>\n";
}
print OUTHTML "</table>\n</body>\n";
print OUTHTML "</html>";
close (OUT);
close (OUTHTML);

