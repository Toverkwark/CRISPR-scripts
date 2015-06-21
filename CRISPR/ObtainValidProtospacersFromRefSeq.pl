use Getopt::Std;
use warnings;
use strict;
require '../SharedTools/FetchGenomicSequence.pl';

sub FetchGenomicSequence($$$);
my $ScriptName="ObtainValidProtospacersFromRefSeq.pl";

#Cas9 cuts here:
#XXXXXXXXXXXXXXXXX|XXXNGG

#This value determines how far the 3'G of the PAM site should be minimally removed from the exon boundary
#If CRISPRs are designed too close to the boundary, indels may result in no phenotype
my $MinimumDistanceToExonBoundary = 0;

#Determine in what part of the protein coding sequence protospacers should be found
my $ProtospacerSearchFraction=1;

#Determine how many nucleotides around the transcription start site, protospacers should be found
my $ProtospacerSearchAroundTSS=100;

my %opts;
my %TargetSites;

getopt( 'ro', \%opts );
die "ERROR in $ScriptName: No RefSeq ID given.\n" unless my $RefSeq = $opts{'r'};
my $RefSeqInfo = `grep -P "$RefSeq\t" ../GenomeInfo/hg19.txt`;
die "ERROR in $ScriptName: RefSeq $RefSeq cannot be found in the database.\n" if !$RefSeqInfo;
die "ERROR in $ScriptName: No Outputfile given.\n" unless my $OutputFile = $opts{'o'};
open (OUT, ">", $OutputFile) or die "Cannot open outpufile $OutputFile\n";
my @RefSeqValues = split( /\t/, $RefSeqInfo );

my $Chromosome = substr( $RefSeqValues[2], 3 );

#The refGene.txt file start sites are always one nt 5' of the actual nt, so compensate for that
my $GeneStart    = $RefSeqValues[4]+1;
my $GeneEnd      = $RefSeqValues[5];
my $ProteinStart = $RefSeqValues[6]+1;
my $ProteinEnd   = $RefSeqValues[7];
my $Orientation = 0;
$Orientation = 1 if($RefSeqValues[3] eq '+');
my $NumberOfExons  = $RefSeqValues[8];
my @ExonStartSites = split( /,/, $RefSeqValues[9]);
foreach (@ExonStartSites) {
		$_++;
}
my @ExonEndSites   = split( /,/, $RefSeqValues[10]);

my $Pos;
my $Neg;

#First, determine the size of the protein coding region
my $CDSSize=0;
for (my $Exon = 0;$Exon < $NumberOfExons; $Exon++) {
	my $ExonStart= ($ExonStartSites[$Exon]);
	my $ExonEnd = ($ExonEndSites[$Exon]); 
	
	#Check if the current exon has coding sequence. If it has, set StartSite and EndSite according to coordinates that fall within protein coding sequence
	if($ProteinStart<=$ExonEnd && $ProteinEnd>=$ExonStart) {
		my $StartSite = ($ProteinStart > $ExonStart) ? $ProteinStart : $ExonStart;
		my $EndSite = ($ProteinEnd < $ExonEnd) ? $ProteinEnd : $ExonEnd;
		$CDSSize = $CDSSize + ($EndSite - $StartSite) + 1; 
	}
}

#Set the limit of how far to look into the CDS
my $CDSLimit = int(($CDSSize * $ProtospacerSearchFraction) + 0.5);
$CDSSize = 0;

#Cycle over the exons and scan for valid target sites
if($Orientation == 0) {
	@ExonStartSites = reverse @ExonStartSites;
	@ExonEndSites = reverse @ExonEndSites;
}

for ( my $Exon = 0 ; $Exon < $NumberOfExons ; $Exon++ ) {
	my $ExonStart= ($ExonStartSites[$Exon]);
	my $ExonEnd = ($ExonEndSites[$Exon]); 
	
	#Check if the current exon has coding sequence. If it has, set StartSite and EndSite according to coordinates that fall within protein coding sequence
	if($ProteinStart<=$ExonEnd && $ProteinEnd>=$ExonStart) {
		my $StartSite = ($ProteinStart > $ExonStart) ? $ProteinStart : $ExonStart;
		my $EndSite = ($ProteinEnd < $ExonEnd) ? $ProteinEnd : $ExonEnd;
		$CDSSize = $CDSSize + ($EndSite - $StartSite) + 1;
		#Check if adding all or part of the current exon falls within the CDS Size limit. 
		#This is the case when either the current CDSSize is smaller than the limit, or the CDSSize without the currently investaged exon investigated is.
		if (($CDSSize<$CDSLimit) || (($CDSSize - (($EndSite - $StartSite)+1)) <$CDSLimit)) {
			#In case adding the current exon exceeds the size limit within the current exon, readjust the search space to fit the limit
			if($CDSSize>$CDSLimit) {
				if($Orientation==1) {
					$EndSite=$StartSite + ($CDSLimit - ($CDSSize - (($EndSite - $StartSite) + 1)));	
				}
				else {
					$StartSite=$EndSite - ($CDSLimit - ($CDSSize - (($EndSite - $StartSite) + 1)));
				}
			}
			
			#Clip off Minimum distance to exon boundary	
			$StartSite=($StartSite > ($ExonStart + $MinimumDistanceToExonBoundary)) ? $StartSite : $ExonStart + $MinimumDistanceToExonBoundary;
			$EndSite=($EndSite < ($ExonEnd - $MinimumDistanceToExonBoundary)) ? $EndSite : $ExonEnd - $MinimumDistanceToExonBoundary;
		
			#To search for forward hits, add 17 nt at the beginning and 6 nt to the end
			#To search for reverse hits, add 17 nt at the end and 6nt at the beginning
			my $ExonSequenceFW = FetchGenomicSequence($Chromosome, $StartSite-17, $EndSite+6);
			my $ExonSequenceRV = FetchGenomicSequence($Chromosome, $StartSite-6, $EndSite+17);
	
			#Search for valid protospacers in the sense strand
			#while ( $ExonSequenceFW =~ /(?<=(.{21}(G|A)))(G)/g ) { //UNCOMMENT THIS LINE FOR INCLUDING NAG PAMs AS VALID PROTOSPACERS
			while ( $ExonSequenceFW =~ m/(?=(.{21}GG))/g ) {
				$TargetSites{substr( $1, 0, 20 )}++;
				$Pos++;
			}
	
			#Search for valid protospacers in the antisense strand
			#while ( $ExonSequenceRV =~ /(?=(.(C|T).{21}))(C)/g ) { //UNCOMMENT THIS LINE FOR INCLUDING NAG PAMs AS VALID PROTOSPACERS
			#while ( $ExonSequenceRV =~ m/(?=(.(C).{21}))(C)/g ) {
			while ( $ExonSequenceRV =~ m/(?=(CC.{21}))/g ) {	
				my $TargetSequence = substr($1,3,20);				
				$TargetSequence =~ tr/ACTG/TGAC/;
				$TargetSequence = reverse($TargetSequence);
				$TargetSites{$TargetSequence}++;
				$Neg++;
			}
		} 
	}
}

#Now, add the potential target sites around the TSS
my $ExonSequenceFW;
my $ExonSequenceRV;
if($Orientation==1) {
	$ExonSequenceFW = FetchGenomicSequence($Chromosome, $ProteinStart-$ProtospacerSearchAroundTSS-17, $ProteinStart+$ProtospacerSearchAroundTSS+6);
	$ExonSequenceRV = FetchGenomicSequence($Chromosome, $ProteinStart-$ProtospacerSearchAroundTSS-6, $ProteinStart+$ProtospacerSearchAroundTSS+17);
}
else {
	$ExonSequenceFW = FetchGenomicSequence($Chromosome, $ProteinEnd-$ProtospacerSearchAroundTSS-17, $ProteinEnd+$ProtospacerSearchAroundTSS+6);
	$ExonSequenceRV = FetchGenomicSequence($Chromosome, $ProteinEnd-$ProtospacerSearchAroundTSS-6, $ProteinEnd+$ProtospacerSearchAroundTSS+17);
}
#Search in sense strand
while ( $ExonSequenceFW =~ m/(?=(.{21}GG))/g ) {
	$TargetSites{substr( $1, 0, 20 )}++;
	$Pos++;
}
#And in antisense strand
while ( $ExonSequenceRV =~ m/(?=(CC.{21}))/g ) {	
	my $TargetSequence = substr($1,3,20);
	$TargetSequence =~ tr/ACTG/TGAC/;
	$TargetSequence = reverse($TargetSequence);
	$TargetSites{$TargetSequence}++;
	$Neg++;
}
			
foreach my $TargetSite ( keys %TargetSites ) {
	#Test for >=4xT in the sequence
	unless ($TargetSite =~ /TTTT/) {
		print OUT "$TargetSite" . "CGG\n";
		print OUT "$TargetSite" . "TGG\n";
		print OUT "$TargetSite" . "AGG\n";
		print OUT "$TargetSite" . "GGG\n";
		#print OUT "$TargetSite" . "CAG\n";
		#print OUT "$TargetSite" . "TAG\n";
		#print OUT "$TargetSite" . "AAG\n";
		#print OUT "$TargetSite" . "GAG\n";
	}
}

print "$Pos Positive and $Neg Negative target sites found.\n";
close (OUT) or die "ERROR in $ScriptName: Cannot close outputfile $OutputFile\n";
