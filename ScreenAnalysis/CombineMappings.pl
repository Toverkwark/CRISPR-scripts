use Getopt::Long;
use Term::ANSIColor;
use strict;
use lib '..';
use LocalSettings;

print "Usage:perl $0 -left -right -leftlib -rightlib -output\n";
my $LeftFile;
my $RightFile;
my $LeftLibraryFile;
my $RightLibraryFile;
my $OutputFile;
my %InsertCounts;
my @Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);

GetOptions(
	"left=s"  => \$LeftFile,
	"right=s" => \$RightFile,
	"leftlib=s" => \$LeftLibraryFile,
	"rightlib=s" => \$RightLibraryFile,
	"output=s" => \$OutputFile
);

open( LEFT, $LeftFile ) or die "ERROR in $0:Left input file $LeftFile is not accessible.\n";
open( RIGHT, $RightFile ) or die "ERROR in $0:Right input file $RightFile is not accessible.\n";
open( LEFTLIB, $LeftLibraryFile ) or die "ERROR in $0:Left library file $LeftLibraryFile is not accessible.\n";
open( RIGHTLIB, $RightLibraryFile ) or die "ERROR in $0:Right library file $RightLibraryFile is not accessible.\n";
open( OUTPUT,">", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";

#Start by reading in the library files
#The format of this file should be [ID],[GENE],[SEQUENCE]
print "Reading left library file\n";
my %LeftLibrary;
my %LeftGenes;
my $LeftInsertsFound;
while ( defined( my $Line = <LEFTLIB> ) ) {
	$LeftInsertsFound++;
	chomp($Line);
	my @values = split( /\,/, $Line );

	#Include this line because excel generated csv files have CRLF line endings
	#$values[2]=substr($values[2],0,length($values[2])-1);
	
	$LeftGenes{$values[0]}=$values[1];
}
close(LEFTLIB) or die "Could not close file $LeftLibraryFile\n";
print "$LeftInsertsFound inserts found in library file $LeftLibraryFile\n";

print "Reading right library file\n";
my %RightLibrary;
my %RightGenes;
my $RightInsertsFound;
while ( defined( my $Line = <RIGHTLIB> ) ) {
	$RightInsertsFound++;
	chomp($Line);
	my @values = split( /\,/, $Line );

	#Include this line because excel generated csv files have CRLF line endings
	#$values[2]=substr($values[2],0,length($values[2])-1);
	
	$RightGenes{$values[0]}=$values[1];
}
close(RIGHTLIB) or die "Could not close file $RightLibraryFile\n";
print "$RightInsertsFound inserts found in library file $RightLibraryFile\n";

#Loop over input files to combine info
my $TotalInsertCounts=0;
my $LeftInsertCounts=0;
while(defined(my $LeftLine=<LEFT>)) {
	my $RightLine=<RIGHT>;
	chomp($LeftLine);
	chomp($RightLine);
	my @values=split(/\t/,$LeftLine);
	my $LeftTag=$values[0];
	my $LeftBarcode=$values[1];
	@values=split(/\t/,$RightLine);
	my $RightTag=$values[0];
	my $RightBarcode=$values[1];
	if($LeftBarcode && $RightBarcode) {
		$InsertCounts{$LeftBarcode}->{$LeftTag}->{$RightTag}++;
		$TotalInsertCounts++;	
	}
	$LeftInsertCounts++;
}
print "$TotalInsertCounts of combined reads mapped out of $LeftInsertCounts for each individual file\n";

#Output the read counts of library hits
print "Writing combined library insert counts\n";
print OUTPUT "Left Tag\tLeft Gene\tRight Tag\tRight Gene";
foreach my $Barcode (@Barcodes) {
	print OUTPUT "\t$Barcode";
}
print OUTPUT "\n";
foreach my $LeftTag (sort {$a cmp $b} keys %LeftGenes) {
	foreach my $RightTag (sort {$a cmp $b} keys %RightGenes) {
		print OUTPUT $LeftTag . "\t" . $LeftGenes{$LeftTag} . "\t" . $RightTag . "\t" . $RightGenes{$RightTag};
		foreach my $Barcode (@Barcodes) {
			if($InsertCounts{$Barcode}->{$LeftTag}->{$RightTag}) {
				print OUTPUT "\t" . $InsertCounts{$Barcode}->{$LeftTag}->{$RightTag};	
			} else {
				print OUTPUT "\t0";
			}
		}
		print OUTPUT "\n";
	}
}

close(OUTPUT) or die "Could not close output file $OutputFile.\n";
close(LEFT) or die "Could not close left input file $LeftFile.\n";
close(RIGHT) or die "Could not close right input file $RightFile.\n";
