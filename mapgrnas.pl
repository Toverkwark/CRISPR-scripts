$LibraryFile="../library.csv";
open( LIBRARY, $LibraryFile ) or die "ERROR in $0:Library file $LibraryFile is not accessible.\n";

#Start by reading in the library file
#The format of this file should be [ID] tab [SEQUENCE]
print "Reading library file\n";
my %Library;
my $InsertsFound;
my $Line=<LIBRARY>;
while ( defined($Line = <LIBRARY> ) ) {
	$InsertsFound++;
	chomp($Line);
	my @values = split( /\t/, $Line );
	$Library{$values[7]} = $values[3];
}
close(LIBRARY) or die "Could not close file $LibraryFile\n";
print "$InsertsFound inserts found in library file $LibraryFile\n";
open (INPUT, "../Stefno.stripped");
open(OUTPUT,">","../Stefano.filtered");
$Line=<INPUT>;
while(defined($Line=<INPUT>)) {
	chomp($Line);
	my @values=split(/\t/,$Line);
	$gRNA=substr($values[0],1);
	if($Library{$gRNA} && substr($values[0],0,1) eq 'G') {
		print OUTPUT $Library{$gRNA};
		for (my $i=0;$i<12;$i++) {
			print OUTPUT "\t" . $values[$i];
		}
		print OUTPUT "\n";
	}
}
close(OUTPUT);
close(INPUT);