use Time::HiRes qw/gettimeofday/;
use Getopt::Std;
use warnings;
use strict;

#This script creates an svg graphical representation of the coding sequence of a gene with CRISPR location mapped onto it
#As input arguments, it takes the refseq id (r) and the species (s) for which to find the info
#As output, it prints the name of the svg file created.
#An html file containing a table to all the CRISPR sites in the svg file is also created
#The whole image is created such that it is 1500 pixels wide. This can best be adjusted by using the viewBox argument in the svg header in the header file

my $HomepagePrefix = "/home/NKI/b.evers/public_html/ikrunc/";
my $HeaderFile = $HomepagePrefix . "scripts/header";
my $FooterFile = $HomepagePrefix . "scripts/footer";
my $ScriptName="ReportTargetSites.pl";
my ($Genome, $RefSeqFile, $QualitiesFileLocation);

#First, set some parameters that determine how the svg file will look
my $YOffset=0;
my $AnimationDistance = 30;
my $TriangleHeight = 10;
my $TriangleWidth = 4;
my $CollisionWidth = 8;
my $LevelOffset = 10;

#Get options given to the script
my %opts;
getopt( 'rsnfa', \%opts );
die "ERROR in $ScriptName: No RefSeq ID given.\n" unless my $RefSeq = $opts{'r'};
die "ERROR in $ScriptName: No species given.\n" unless my $Species = $opts{'s'};
my $ReportNumberOfNucleotides;
$ReportNumberOfNucleotides=-1 unless $ReportNumberOfNucleotides = $opts{'n'};
my $ReportFractionOfTranscript;
$ReportFractionOfTranscript=1 unless $ReportFractionOfTranscript = $opts{'f'};
my $ReportNucleotidesAroundStartSite;
$ReportNucleotidesAroundStartSite=0 unless $ReportNucleotidesAroundStartSite = $opts{'a'};

#Assign proper genome id
$Species = lc($Species);
if ($Species eq 'mouse') {
	$Genome = "mm10";
	$RefSeqFile = "/home/NKI/b.evers/mm10/scripts/RefSeq/refGene.txt";
	$QualitiesFileLocation = "/home/NKI/b.evers/mm10/output";
}
else {
	if ($Species eq 'human') {
		$Genome = "hg19";
		$RefSeqFile = "/home/NKI/b.evers/hg19/scripts/RefSeq/refGene.txt";
		$QualitiesFileLocation = "/home/NKI/b.evers/hg19/output";
	}
	else {
		die "ERROR in $ScriptName: Species $Species is not currently handled by the iKRUNC scripts\n";
	}
}

#Delete all the temporary files that were created in the last 5 minutes, to keep the temp directory small
my $DeleteTempFilesCommand = "find " . $HomepagePrefix . "temp/* -mmin +5 -exec rm {} \\;";
`$DeleteTempFilesCommand`;

#Filter the qualities.4 file for relevant target sites. As output, use a temporary file
my $timestamp = gettimeofday;
print "temp/" . $timestamp . ".svg";
`perl $HomepagePrefix/scripts/FilterAndSelect.pl -i $QualitiesFileLocation/$RefSeq.qualities.4 -o $HomepagePrefix/temp/$timestamp -c $ReportFractionOfTranscript -e $ReportNucleotidesAroundStartSite -n $ReportNumberOfNucleotides -r $RefSeq -s $Species`;

#Open this temporary file and for output open a file with a similar name, but with extension .svg for the image and one with extension .svg.html for the linked table
my $InputFile=$HomepagePrefix . "temp/" . $timestamp . ".seq";
open (IN, $InputFile) or die "ERROR in $ScriptName: Cannot open inputfile $InputFile\n";
open (OUT, ">", $HomepagePrefix . "temp/" . $timestamp . ".svg") or die "ERROR in $ScriptName: Cannot create svg file\n";
open (OUTHTML, ">", $HomepagePrefix . "temp/" . $timestamp . ".svg.html") or die "ERROR in $ScriptName: Cannot create svg html file\n";

#First read in all intron/exon information of the gene we're plotting
my $RefSeqInfo = `grep -P "$RefSeq\t" $RefSeqFile`;
die "ERROR in $ScriptName: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
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
	my $TargetChromosome = $TargetSites[1];
	my $TargetOrientation = $TargetSites[2];
	my $TargetCutSite = $TargetSites[3];
		
	#Verify that the target is in the gene, or within 250nt of the TSS
	if ($TargetChromosome eq $Chromosome && ($TargetCutSite >= $GeneStart-250 && $TargetCutSite <= $GeneEnd + 250)) {
		@{$DisplayObjects{$TargetOrientation}->{$TargetCutSite}}=@TargetSites;
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
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[13]=$RelativeMarkerPosition;
	}
}

#Determine and assign colors
foreach my $Orientation (keys %DisplayObjects) {
	foreach my $TargetCutSite (keys $DisplayObjects{$Orientation}) {
		my $TriangleValue=$DisplayObjects{$Orientation}->{$TargetCutSite}->[12];
		$TriangleValue = int (511*($TriangleValue));
		my $Blue = 255;
		my $Red = $TriangleValue;
		my $Green = $TriangleValue;
		if ($TriangleValue >= 256) {
			$Blue = 510-$TriangleValue;
			$Red=255;
			$Green=255;
		}
		my $TriangleColor = "$Red,$Green,$Blue"; 
		$DisplayObjects{$Orientation}->{$TargetCutSite}->[14]=$TriangleColor;
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
			my $RelativeMarkerPosition=$DisplayObjects{$Orientation}->{$TargetCutSite}->[13];
			if(!($DisplayObjects{$Orientation}->{$TargetCutSite}->[15])) {
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
					$DisplayObjects{$Orientation}->{$TargetCutSite}->[15]=$CollisionLevel;
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
	foreach my $TargetCutSite (sort {$DisplayObjects{$Orientation}->{$a}->[4] <=> $DisplayObjects{$Orientation}->{$b}->[4]} keys $DisplayObjects{$Orientation}) {		
		my @TargetSites = @{$DisplayObjects{$Orientation}->{$TargetCutSite}};		
		my $TargetChromosome = $TargetSites[1];
		my $TargetOrientation = $TargetSites[2];
		my $TargetCutSite = $TargetSites[3];
		my $TargetScore = $TargetSites[4];
		my $TargetSequence=$TargetSites[5];
		my $TargetNumberOfIdentical3PrimeTargets=$TargetSites[6];
		my $TargetNumberOfIdentical3PrimeTargetsNearExons=$TargetSites[7];
		my $TargetDegree=$TargetSites[8];
		my $TargetClosestRelatives=$TargetSites[9];
		my $TargetClosestRelativesNearExons=$TargetSites[10];
		my $TargetLabel = $TargetSites[11];
		$TargetID=$TargetID+1;		
		my $RelativeMarkerPosition=$TargetSites[13];
		my $TriangleColor=$TargetSites[14];
		my $LayerOffset = $LevelOffset*(($TargetSites[15]) - 1);
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
		$TableObjects{$TargetID}->[9]=$TargetScore;
	}
} 

#Calculate the total height of the figure:
my $TotalHeight=$YOffset+$ExonOffset+($LevelOffset*($MaxCollisionLevelAntisenseStrand-1))+$TriangleHeight+$AnimationDistance;

#Start writing the svg file
#First, write the header
print OUT '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" onload="init(evt)" width="100%" height="100%" viewBox="0 0 1500 ' . $TotalHeight . '" id="svgimage">\n';
print OUT `cat $HeaderFile`;

#Write the image file
print OUT $SVGFile;

#Print the footer
print OUT `cat $FooterFile`;

#Print the table file
print OUTHTML "<!DOCTYPE html>\n<html lang='en'>\n\t<head>\n\t\t<meta charset='utf-8'>\n\t\t<link rel='stylesheet' href='../style.css' type='text/css'></link>\n\t</head>\n";
print OUTHTML "<body>\n<table id='svgtable' class='svgtable' width='100%'>\n";
print OUTHTML "\t<th>Rank</th>\n";
print OUTHTML "\t<th>Chromosome</th>\n";
print OUTHTML "\t<th>Orientation</th>\n";
print OUTHTML "\t<th>Position</th>\n";
print OUTHTML "\t<th>Sequence</th>\n";
print OUTHTML "\t<th>Identical 3'12nt sites</th>\n";
print OUTHTML "\t<th>Identical 3'12nt sites near exons</th>\n";
print OUTHTML "\t<th>Off-target #mismatches</th>\n";
print OUTHTML "\t<th># sites</th>\n";
print OUTHTML "\t<th># sites near exons</th>\n";
my $Rank=0;
foreach my $TableRow (sort {$TableObjects{$b}->[9] <=> $TableObjects{$a}->[9]} keys %TableObjects) {
	$Rank=$Rank+1;
	print OUTHTML "\t<tr id='" . $TableRow . ".table' onclick=parent.ClickTableRow('" . $TableRow . "')>\n";
	print OUTHTML "\t\t<td>" . $Rank . "</td>\n";
	print OUTHTML "\t\t<td>" . $Chromosome . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[1] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[2] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[3] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[4] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[5] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[6] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[7] . "</td>\n";
	print OUTHTML "\t\t<td>" . $TableObjects{$TableRow}->[8] . "</td>\n";
}
print OUTHTML "</table>\n</body>\n";
print OUTHTML "</html>";
close (OUT);
close (OUTHTML);

