$InputFile="../3110-1_S1_L001_R1_001.fastq";
#$InputFile="../small";
$OutputFile="barcode";
@Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);
open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
my %FileHandlers;
for (my $i=1;$i<=4;$i++) {
	open ($FileHandlers{$i}, ">", ($OutputFile . "." . $i));
}
my $RecordsAnalyzed=0;
while (defined (my $Line1=<INPUT> )) {
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 100000 ) ) {
		print "Analyzing record $RecordsAnalyzed of inputfile $InputFile\n";
	}
	my $Line2=<INPUT>;
	my $Line3=<INPUT>;
	my $Line4=<INPUT>;
	my $Barcode=substr($Line2,0,6);
	if ( grep( /$Barcode/, @Barcodes ) ) {
		my $PrintTo=0;
		$PrintTo=1 if($Barcode eq $Barcodes[0] || $Barcode eq $Barcodes[1] || $Barcode eq $Barcodes[2]);
		$PrintTo=2 if($Barcode eq $Barcodes[3] || $Barcode eq $Barcodes[4] || $Barcode eq $Barcodes[5]);
		$PrintTo=3 if($Barcode eq $Barcodes[6] || $Barcode eq $Barcodes[7] || $Barcode eq $Barcodes[8]);
		$PrintTo=4 if($Barcode eq $Barcodes[9] || $Barcode eq $Barcodes[10] || $Barcode eq $Barcodes[11]);
		if($PrintTo>0) {
			print {$FileHandlers{$PrintTo}} $Line1;
			print {$FileHandlers{$PrintTo}} $Line2;
			print {$FileHandlers{$PrintTo}} $Line3;
			print {$FileHandlers{$PrintTo}} $Line4;
		}
	}	
}
for (my $i=1;$i<=4;$i++) {
	close($FileHandlers{$i});
}
close (INPUT) or die "Could not close input file $InputFile\n";