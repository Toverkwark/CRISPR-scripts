use Getopt::Long;

GetOptions(
	"input=s"  => \$InputFile,
);
open( INPUT, $InputFile ) or die "ERROR in $0:Input file $InputFile is not accessible.\n";
$OutputFile=($InputFile . ".barcode");
@Barcodes = qw(CGTGAT ACATCG GCCTAA TGGTCA CACTGT ATTGGC GATCTG TCAAGT CTGATC AAGCTA GTAGCC TACAAG);

my %FileHandlers;
for (my $i=0;$i<12;$i++) {
	open ($FileHandlers{$Barcodes[$i]}, ">", ($OutputFile . "." . ($i+1)));
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
		my $PrintTo=$Barcode;
		print {$FileHandlers{$PrintTo}} $Line1;
		print {$FileHandlers{$PrintTo}} $Line2;
		print {$FileHandlers{$PrintTo}} $Line3;
		print {$FileHandlers{$PrintTo}} $Line4;
	}	
}
for (my $i=1;$i<=4;$i++) {
	close($FileHandlers{$i});
}
close (INPUT) or die "Could not close input file $InputFile\n";
