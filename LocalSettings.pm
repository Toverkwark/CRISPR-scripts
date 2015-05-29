package LocalSettings;
require Exporter;
use strict;
our @ISA=qw(Exporter);
our @EXPORT=qw(getconfig);
my %confighash;
my $ExpectedLeadingSequence;
my $ExpectedTrailingSequence;

########################################################################################################################################################################################
########################################################################################################################################################################################
#General settings
my $Location='Cluster';
my $ScreenType='GIN_A';
my $ExpectedTrailingNucleotides=0; #Set at 0 to include all 
my $BowtieLocation;
my $HumanGenomeLocation;
my $IndexedHumanGenomeLocation;
my $NumberOfCoresToUse = 64;
my $BarcodeOffset=0;
my $BarcodeLength=6;
my $ExpectedInsertLength=20;
my $ErrorThresholdLeading = 10; #This number of mutations or indels can be present in the leading  sequences
my $ErrorThresholdTrailing = 10; #This number of mutations or indels can be present in the trailing sequences
########################################################################################################################################################################################
########################################################################################################################################################################################

#Set relevant expected sequences
if($ScreenType eq 'GIN_A') {
        #For CRISPRi Libraries clone into pGIN_A:
        $ExpectedLeadingSequence = "CCCTTGGAGAAAAGCCTTGTTTG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
        $ExpectedTrailingSequence = "GTTTAAGAGCTAGAAA"; #Sequence that is expected to come after the gRNA/shRNA sequence
}
if($ScreenType eq 'Geckov2') {
	#For GECKO v2 Libraries:
	$ExpectedLeadingSequence = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
	$ExpectedTrailingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"; #Sequence that is expected to come after the gRNA/shRNA sequence
}
if($ScreenType eq 'iKRUNCv1') {
	#For iKRUNC v1 Libraries:
	$ExpectedLeadingSequence = "CCCTATCAGTGATAGAGACTCGAG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
	$ExpectedTrailingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT"; #Sequence that is expected to come after the gRNA/shRNA sequence
}
if($ScreenType eq 'iKRUNCv2Short') {
	#For iKRUNC v2 short Libraries:
	$ExpectedLeadingSequence = "CCCTATCAGTGATAGAGACTCGAG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
	$ExpectedTrailingSequence = "GTTTAAGAGCTAGAAATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT"; #Sequence that is expected to come after the gRNA/shRNA sequence
}
if($ScreenType eq 'iKRUNCv2Long') {
	#For iKRUNC v2 long Libraries:
	$ExpectedLeadingSequence = "CCCTATCAGTGATAGAGACTCGAG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
	$ExpectedTrailingSequence = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT"; #Sequence that is expected to come after the gRNA/shRNA sequence
}
if($ScreenType eq 'TRC') {
	#For TRC Libraries (ALSO ADJUST EXPECTED INSERT LENGTH):
	$ExpectedLeadingSequence = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCGG"; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
	$ExpectedTrailingSequence = "TTTTT"; #Sequence that is expected to come after the gRNA/shRNA sequence
	$ExpectedInsertLength=21+6+21;
}
if($ExpectedTrailingNucleotides > 0) {
	$ExpectedTrailingSequence=substr($ExpectedTrailingSequence,0,$ExpectedTrailingNucleotides);
}	

if($Location eq 'Home') {
	$BowtieLocation="/home/bastiaan/data/genomestuff/bowtie2-2.1.0/bowtie2";
	$HumanGenomeLocation="/home/bastiaan/data/genomestuff/hg19";
	$IndexedHumanGenomeLocation="/home/bastiaan/data/genomestuff/hg19-index";
}
else {
	if($Location eq 'Work') {
	        $BowtieLocation="/media/Data/iKRUNC/bowtie2-2.1.0/bowtie2";
        	$HumanGenomeLocation="/media/Data/iKRUNC/hg19";
	        $IndexedHumanGenomeLocation="/media/Data/iKRUNC/hg19-index";
	}
	else {
	        if($Location eq 'Cluster') {
			use lib "/home/NKI/b.evers/perl5/lib/perl5";
	                $BowtieLocation="bowtie2";
        	        $HumanGenomeLocation="/home/NKI/b.evers/Genomes/hg19/genome-fasta";
                	$IndexedHumanGenomeLocation="/home/NKI/b.evers/Genomes/hg19/genome-index";
		}
        }

}

%confighash=(
	"Bowtie" => $BowtieLocation, 
	"HumanGenome" => $HumanGenomeLocation,
	"IndexedHumanGenome" => $IndexedHumanGenomeLocation,
	"NumberOfCoresToUse" => $NumberOfCoresToUse,
	"BarcodeOffset" => $BarcodeOffset,
	"BarcodeLength" => $BarcodeLength,	
	"ExpectedInsertLength" => $ExpectedInsertLength,
	"ErrorThresholdLeading" => $ErrorThresholdLeading,
	"ErrorThresholdTrailing" => $ErrorThresholdTrailing,
	"ExpectedLeadingSequence" => $ExpectedLeadingSequence,
	"ExpectedTrailingSequence" => $ExpectedTrailingSequence	
);

sub getconfig {
	return %confighash;
}
