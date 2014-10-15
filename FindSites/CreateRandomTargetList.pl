use Getopt::Long;
use warnings;
use strict;
my $OutputFile;
GetOptions(
	"output=s"  => \$OutputFile,
);
open(OUT,">", $OutputFile) or die "ERROR in $0:Cannot open inputfile $OutputFile\n";
for(my $i=1;$i<=100000000;$i++) {
	my $Sequence;
	for (my $j=1;$j<=20;$j++) {
		my $Getal=int(rand(4));
		$Sequence=$Sequence . "A" if($Getal==0);
		$Sequence=$Sequence . "C" if($Getal==1);
		$Sequence=$Sequence . "T" if($Getal==2);
		$Sequence=$Sequence . "G" if($Getal==3);
	}
	print OUT $Sequence . "\n";
}
close(OUT) or die "ERROR in $0:Cannot close outputfile $OutputFile\n";
