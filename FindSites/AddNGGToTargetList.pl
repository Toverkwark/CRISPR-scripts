use Getopt::Long;
use warnings;
use strict;
my $InputFile;
GetOptions(
	"input=s"  => \$InputFile,
);
open(IN,$InputFile) or die "ERROR in $0:Cannot open inputfile $InputFile\n";
my $OutputFile=$InputFile . ".potentialtargets";
open(OUT,">", $OutputFile) or die "ERROR in $0:Cannot open inputfile $OutputFile\n";
my $RecordsAnalyzed=0;
while(defined(my $Line=<IN>)) {
	$RecordsAnalyzed++;
	if ( !( $RecordsAnalyzed % 10000 ) ) {
		print "Analyzing record $RecordsAnalyzed of inputfile $InputFile\n";
	}
	chomp($Line);
	print OUT "$Line" . "CGG\n";
	print OUT "$Line" . "TGG\n";
	print OUT "$Line" . "AGG\n";
	print OUT "$Line" . "GGG\n";
}
close(IN) or die "ERROR in $0:Cannot close inputfile $InputFile\n";
close(OUT) or die "ERROR in $0:Cannot close outputfile $OutputFile\n";
