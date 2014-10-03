use Getopt::Std;
use warnings;
use strict;

my %opts;
my %TargetSites;
my $QuerySequence;
my $InputFile;
my $OutputFile;
my $Pos;
my $Neg;

print("Usage:perl $0 -i inputfile [-o outputfile]\n");
getopt( 'io', \%opts );
die "ERROR in $0: No inputfile given.\n" unless $InputFile = $opts{'i'};
$OutputFile=($InputFile . ".protospacers") unless $OutputFile = $opts{'o'};
open (IN,$InputFile) or die "Error in $0:Cannot open inputfile $InputFile\n";
open (OUT, ">", $OutputFile) or die "Error in $0:Cannot open outpufile $OutputFile\n";

#Build the query sequence
while (defined(my $Line=<IN>)) {
	chomp($Line);
	$QuerySequence = $QuerySequence . $Line;
}

#Search for valid protospacers in the sense strand
while ( $QuerySequence =~ /(?<=(.{22}))(G)/g ) {
	$TargetSites{substr( $1, 0, 20 )}++;
	$Pos++;
}
	
#Search for valid protospacers in the antisense strand
while ( $QuerySequence =~ /(?=(.{23}))(C)/g ) {
	my $TargetSequence = substr($1,3,20);
	$TargetSequence =~ tr/ACTG/TGAC/;
	$TargetSequence = reverse($TargetSequence);
	$TargetSites{$TargetSequence}++;
	$Neg++;
}

			
foreach my $TargetSite ( sort {$a cmp $b} keys %TargetSites ) {
	#Test for >5xT in the sequence
	unless ($TargetSite =~ /TTTTT/) {
		print OUT "$TargetSite" . "AAG\n";
		print OUT "$TargetSite" . "CAG\n";
		print OUT "$TargetSite" . "GAG\n";
		print OUT "$TargetSite" . "TAG\n";
		print OUT "$TargetSite" . "ACG\n";
		print OUT "$TargetSite" . "CCG\n";
		print OUT "$TargetSite" . "GCG\n";
		print OUT "$TargetSite" . "TCG\n";
		print OUT "$TargetSite" . "AGG\n";
		print OUT "$TargetSite" . "CGG\n";
		print OUT "$TargetSite" . "GGG\n";
		print OUT "$TargetSite" . "TGG\n";
		print OUT "$TargetSite" . "ATG\n";
		print OUT "$TargetSite" . "CTG\n";
		print OUT "$TargetSite" . "GTG\n";
		print OUT "$TargetSite" . "TTG\n";
	}
}

print "$Pos Positive and $Neg Negative target sites found.\n";
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";
