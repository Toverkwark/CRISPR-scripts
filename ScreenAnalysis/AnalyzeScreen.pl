use Getopt::Long;
use Term::ANSIColor;
use strict;
#use threads;
#use threads::shared;
use lib '..';
use LocalSettings;
use Parallel::ForkManager;
require "Damerau.pl";
require "ProcessReads.pl";
sub MatchBarcode($@);
sub ScoreTwoStrings($$);

print "Usage:perl $0 -input -output -report -library -keepmapfiles\n-input\tName of input file\n-output\tName of output file. Default is inputfile.stripped\n-report\tName of report file. Default is inputfile.report\n-library\tName of library file to which inserts are mapped\n-keepmapfiles\tGive Y as argument to make sure map files are kept for paired end matching\n";

my $StartTime=time;

#Define screen specific settings
my %LocalSettings=getconfig();
my $NumberOfThreads = $LocalSettings{'NumberOfCoresToUse'};
print "Number of cores to use:" . $NumberOfThreads . "\n";
my $BarcodeOffset = $LocalSettings{'BarcodeOffset'}; #Position of start of barcode
my $BarcodeLength = $LocalSettings{'BarcodeLength'}; #Number of nucleotides that the barcode is long
my $ExpectedInsertLength = $LocalSettings{'ExpectedInsertLength'}; #Number of nucleotides of the insert between leading and trailing sequence
my $ExpectedLeadingSequence = $LocalSettings{'ExpectedLeadingSequence'}; #Sequence that is expected to come between the barcode and the start of the gRNA/shRNA sequence
my $ExpectedTrailingSequence = $LocalSettings{'ExpectedTrailingSequence'}; #Sequence that is expected to come after the gRNA/shRNA sequence
my $ErrorThresholdLeading = $LocalSettings{'ErrorThresholdLeading'}; #This number of mutations or indels can be present in the leading  sequences
my $ErrorThresholdTrailing = $LocalSettings{'ErrorThresholdTrailing'}; #This number of mutations or indels can be present in the trailing sequences
my @Barcodes;
@Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
my %Results;
my %LeadingErrors;
my %TrailingErrors;
my %InsertLengths;
my %QualitiesByBarcode;
my %InsertCounts;
my %PerfectInsertCounts;
my $InputFile;
my $OutputFile;
my $ReportFile;
my $LibraryFile;
my $KeepMapFiles;
my $RecordsAnalyzed=0;
my $NotAnalyzed="";

GetOptions(
	"input=s"  => \$InputFile,
	"output=s" => \$OutputFile,
	"report=s" => \$ReportFile,
	"library=s" => \$LibraryFile,
	"keepmapfiles=s" => \$KeepMapFiles
);

if ( !$OutputFile ) {
	$OutputFile = $InputFile . ".output";
}

if ( !$ReportFile ) {
	$ReportFile = $InputFile . ".report";
}

open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
open( OUTPUT, ">", $OutputFile ) or die "ERROR in $0:Output file $OutputFile is not accessible.\n";
open( PERFECTOUTPUT, ">", ($OutputFile . ".perfect")) or die "ERROR in $0:Output file $OutputFile.perfect is not accessible.\n";
open( REPORT, ">", $ReportFile ) or die "ERROR in $0:Report file $ReportFile is not accessible.\n";
open( LIBRARY, $LibraryFile ) or die "ERROR in $0:Library file $LibraryFile is not accessible.\n";
if($KeepMapFiles eq 'Y') {
	open( MAPPEDOUTPUT, ">", $OutputFile . ".mapped") or die "ERROR in $0:Output file $OutputFile.mapped is not accessible.\n";
	open( MAPPEDPERFECTOUTPUT, ">", ($OutputFile . ".mapped.perfect")) or die "ERROR in $0:Output file $OutputFile.mapped.perfect is not accessible.\n";	
}

#Start by reading in the library file
#The format of this file should be [ID],[GENE],[SEQUENCE]
print "Reading library file\n";
my %Library;
my %Genes;
my $InsertsFound;
while ( defined( my $Line = <LIBRARY> ) ) {
	$InsertsFound++;
	chomp($Line);
	my @values = split( /\,/, $Line );

	#Include this line because excel generated csv files have CRLF line endings
	#$values[2]=substr($values[2],0,length($values[2])-1);
	
	$Library{$values[2]} = $values[0];
	$Genes{$values[0]}=$values[1];
}
close(LIBRARY) or die "Could not close file $LibraryFile\n";
print "$InsertsFound inserts found in library file $LibraryFile\n";

#Split up file in the same number of subfiles as we're gonna use threads to analyze
my %FileHandlers;
for (my $i=1;$i<=$NumberOfThreads;$i++) {
	open ($FileHandlers{$i}, ">", ($InputFile . "." . $i));
}
my $Thread=1;
while (defined (my $Line=<INPUT> )) {
	$Thread=1 if($Thread>$NumberOfThreads);
	print {$FileHandlers{$Thread}} $Line;
	$Line=<INPUT>;
	print {$FileHandlers{$Thread}} $Line;
	$Line=<INPUT>;
	print {$FileHandlers{$Thread}} $Line;
	$Line=<INPUT>;
	print {$FileHandlers{$Thread}} $Line;
	$Thread++;
}
close (INPUT) or die "Could not close input file $InputFile";

for (my $i=1;$i<=$NumberOfThreads;$i++) {
	close($FileHandlers{$i});
}
my $ProcessManager=new Parallel::ForkManager($NumberOfThreads);
for ($Thread=1;$Thread<=$NumberOfThreads;$Thread++) {
	$ProcessManager->start and next;
	ProcessReads($InputFile . "." . $Thread,$BarcodeLength,$BarcodeOffset,$ExpectedInsertLength,$ExpectedLeadingSequence,$ExpectedTrailingSequence,$ErrorThresholdLeading,$ErrorThresholdTrailing,\%Library,@Barcodes);
	$ProcessManager->finish;
}
$ProcessManager->wait_all_children;

#Loop through all temp files to obtain data
open( NOTANALYZED, ">", $InputFile . ".notanalyzed") or die "ERROR in $0:File " . $InputFile . ".analyzed is not accessible.\n";
for ($Thread=1;$Thread<=$NumberOfThreads;$Thread++) {
	print "Reading results from tempfile $Thread\n";
	open(INPUT,($InputFile . "." .$Thread . ".tmp")) or die "Could not open temporary result file $InputFile.$Thread.tmp\n";
	#Read Results
	my $Line=<INPUT>;
	while(($Line=<INPUT>) ne "***LEADING ERRORS***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$Results{$Values[0]}->[$Values[1]]+=$Values[2];
	}
	#Read Leading Errors
	while(($Line=<INPUT>) ne "***TRAILING ERRORS***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$LeadingErrors{$Values[0]}->{$Values[1]}+=$Values[2];
	}
	#Read Trailing Errors
	while(($Line=<INPUT>) ne "***INSERT LENGTHS***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$TrailingErrors{$Values[0]}->{$Values[1]}+=$Values[2];
	}
	#Read Insert Lengths
	while(($Line=<INPUT>) ne "***INSERT COUNTS***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$InsertLengths{$Values[0]}+=$Values[1];
	}
	#Read Insert Counts
	while(($Line=<INPUT>) ne "***PERFECT INSERT COUNTS***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$InsertCounts{$Values[0]}->{$Values[1]}+=$Values[2];
	}
	#Read Perfect Insert Counts
	while(($Line=<INPUT>) ne "***QUALITIES BY BARCODE***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$PerfectInsertCounts{$Values[0]}->{$Values[1]}+=$Values[2];
	}
	#Read Qualities By Barcode
	while(($Line=<INPUT>) ne "***RECORDS ANALYZED***\n") {
		chomp($Line);
		my @Values=split(/\t/,$Line);
		$QualitiesByBarcode{$Values[0]}->{$Values[1]}+=$Values[2];
	}
	#Read Records Analyzed
	$Line=<INPUT>;
	chomp($Line);
	$RecordsAnalyzed=$RecordsAnalyzed+$Line;
	#Read Records not analyzed
	$Line=<INPUT>;
	while (defined ($Line=<INPUT>)) {
		print NOTANALYZED $Line;
	}
	close(INPUT) or die "Could not open temporary result file $InputFile$Thread.tmp\n";
	if($KeepMapFiles eq 'Y') {
		#Insert code here to combine intermediate mapped files into one
	}
}
close(NOTANALYZED) or die "Could not close file $InputFile.notanalyzed\n";

#Delete intermediate files
for ($Thread=1;$Thread<=$NumberOfThreads;$Thread++) {
	unlink($InputFile . "." . $Thread);
	unlink($InputFile . "." . $Thread . ".tmp");
	unlink($InputFile . "." . $Thread . ".mapped");
	unlink($InputFile . "." . $Thread . ".perfect.mapped");
}

#Print report of read filtering per barcode
print "Writing filtering steps per barcode\n";
my @Totals;
print REPORT "Filtering results\n";
print REPORT "Barcode\tBarcodes Found\tPerfect Barcodes Found\tLeader found\tPerfect leaders found\tTrailers found\tPerfect trailers found\tLeader and trailes found\tCorrect insert length\tInserts in library\tPerfect Leader and trailes found\tPerfect Correct insert length\tPerfect Inserts in library\n";
foreach my $Barcode ( @Barcodes ) {
	print REPORT "$Barcode\t";
	for ( my $counter = 0 ; $counter <= 11 ; $counter++ ) {
		print REPORT sprintf( "%i", $Results{$Barcode}->[$counter] ) . "\t";
		$Totals[$counter] += $Results{$Barcode}->[$counter];
	}
	print REPORT "\n";
}
print REPORT sprintf( "TOTAL\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", $Totals[0], $Totals[1], $Totals[2], $Totals[3], $Totals[4], $Totals[5], $Totals[6], $Totals[7], $Totals[8], $Totals[9], $Totals[10], $Totals[11]);

print REPORT sprintf(
	"\n\t\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%\t%4.2f%%",
	100 * $Totals[1] / $Totals[0],
	100 * $Totals[2] / $Totals[0],
	100 * $Totals[3] / $Totals[0],
	100 * $Totals[4] / $Totals[0],
	100 * $Totals[5] / $Totals[0],
	100 * $Totals[6] / $Totals[0],
	100 * $Totals[7] / $Totals[0],
	100 * $Totals[8] / $Totals[0],
	100 * $Totals[9] / $Totals[0],
	100 * $Totals[10] / $Totals[0],
	100 * $Totals[11] / $Totals[0],
);
print REPORT "\n\nAnalyzed reads\t" . sprintf( "%i", $RecordsAnalyzed ) . "\n\n";

#Print quality information
print "Writing quality information to report\n";
print REPORT "Quality Information\n";
foreach my $Barcode ( @Barcodes ) {
	print REPORT "$Barcode\t";
	foreach my $CharNumber (sort {$a<=>$b} keys %{$QualitiesByBarcode{$Barcode}}) {
		print REPORT sprintf("%4.2f", $QualitiesByBarcode{$Barcode}->{$CharNumber} / $Results{$Barcode}->[0]) . "\t";
	}
	print REPORT "\n";
}

#Output the insert size distribution
print "Writing insert size distribution to report\n";
print REPORT "\nInsert sizes\n";
foreach my $InsertLength (sort {$a <=> $b} keys %InsertLengths) {
	print REPORT "Size\t" . ($InsertLength) . "\t" . sprintf("%4.3f%%",(100*$InsertLengths{$InsertLength} / $Totals[6])) . "\n";
}
print REPORT "\n";

#Output the leading error information
print "Writing leading error information to report\n";
print REPORT "Leading errors\nPosition\tInsertions\tDeletions\tMutations\n";
foreach my $Position (sort {$a <=> $b} keys %LeadingErrors) {
	my $TotalErrors=0;
	print REPORT "$Position\t";
	my @ErrorTypes = ('Insertion', 'Deletion', 'Mutation');
	foreach my $ErrorType (@ErrorTypes) {
		$TotalErrors=$TotalErrors+$LeadingErrors{$Position}->{$ErrorType};
		print REPORT sprintf("%4.3f%%",(100*$LeadingErrors{$Position}->{$ErrorType} / $Totals[2])) . "\t";
	}
	print REPORT sprintf("%4.3f%%",(100*$TotalErrors / $Totals[2])) . "\n";
}
print REPORT "\n";

#Output the trailing error information
print "Writing trailing error information to report\n";
print REPORT "Trailing errors\nPosition\tInsertions\tDeletions\tMutations\n";
foreach my $Position (sort {$a <=> $b} keys %TrailingErrors) {
	my $TotalErrors=0;
	print REPORT "$Position\t";
	my @ErrorTypes = ('Insertion', 'Deletion', 'Mutation');
	foreach my $ErrorType (@ErrorTypes) {
		$TotalErrors=$TotalErrors+$TrailingErrors{$Position}->{$ErrorType};
		print REPORT sprintf("%4.3f%%",(100*$TrailingErrors{$Position}->{$ErrorType} / $Totals[4])) . "\t";
	}
	print REPORT sprintf("%4.3f%%",(100*$TotalErrors / $Totals[4])) . "\n";
}
print REPORT "\n";


#Output the read counts of library hits
print "Writing filtered library insert counts\n";
print OUTPUT "LibraryID\tGene";
foreach my $Barcode (@Barcodes) {
	print OUTPUT "\t$Barcode";
}
print OUTPUT "\n";
foreach my $InsertSequence (sort {$Library{$a} cmp $Library{$b}} keys %Library) {
	my $LibraryGene=substr($Library{$InsertSequence},0,index($Library{$InsertSequence},'-'));
	print OUTPUT $Library{$InsertSequence} . "\t" . $Genes{$Library{$InsertSequence}};
	foreach my $Barcode (@Barcodes) {
		if($InsertCounts{$InsertSequence}->{$Barcode}) {
			print OUTPUT "\t" . $InsertCounts{$InsertSequence}->{$Barcode};	
		}
		else {
			print OUTPUT "\t0";
		}
	}
	print OUTPUT "\n";
}

#Output the read counts of library hits that are 100% perfect
print "Writing filtered library insert counts that are 100% perfect\n";
print PERFECTOUTPUT "LibraryID\tSequence";
foreach my $Barcode (@Barcodes) {
	print PERFECTOUTPUT "\t$Barcode";
}
print PERFECTOUTPUT "\n";
foreach my $InsertSequence (sort {$Library{$a} cmp $Library{$b}} keys %Library) {
	my $LibraryGene=substr($Library{$InsertSequence},0,index($Library{$InsertSequence},'-'));
	print PERFECTOUTPUT $Library{$InsertSequence} . "\t" . $Genes{$Library{$InsertSequence}};
	foreach my $Barcode (@Barcodes) {
		if($PerfectInsertCounts{$InsertSequence}->{$Barcode}) {
			print PERFECTOUTPUT "\t" . $PerfectInsertCounts{$InsertSequence}->{$Barcode};	
		}
		else {
			print PERFECTOUTPUT "\t0";
		}
	}
	print PERFECTOUTPUT "\n";
}

close(OUTPUT) or die "Could not close output file $OutputFile.\n";
close(PERFECTOUTPUT) or die "Could not close output file $OutputFile.perfect.\n";
close(REPORT) or die "Could not close report file $ReportFile.\n";

my $EndTime=time;
print "Time spent analyzing:" . ($EndTime-$StartTime) . " seconds. Analyzed $RecordsAnalyzed sequences, or " . ($RecordsAnalyzed/($EndTime-$StartTime)) . " sequences per second\n";




