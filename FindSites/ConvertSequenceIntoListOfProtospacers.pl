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
while ( $QuerySequence =~ /(?<=(.{21}))(G)/g ) {
	$TargetSites{substr( $1, 0, 20 )}++;
	$Pos++;
}
	
#Search for valid protospacers in the antisense strand
while ( $QuerySequence =~ /(?=(.{22}))(C)/g ) {
	my $TargetSequence = substr($1,2,20);
	$TargetSequence =~ tr/ACTG/TGAC/;
	$TargetSequence = reverse($TargetSequence);
	$TargetSites{$TargetSequence}++;
	$Neg++;
}

			
foreach my $TargetSite ( sort {$a cmp $b} keys %TargetSites ) {
	#Test for >5xT in the sequence
	unless ($TargetSite =~ /TTTTT/) {
		print OUT "$TargetSite" . "AGA\n";
		print OUT "$TargetSite" . "CGA\n";
		print OUT "$TargetSite" . "GGA\n";
		print OUT "$TargetSite" . "TGA\n";
		print OUT "$TargetSite" . "AGC\n";
		print OUT "$TargetSite" . "CGC\n";
		print OUT "$TargetSite" . "GGC\n";
		print OUT "$TargetSite" . "TGC\n";
		print OUT "$TargetSite" . "AGG\n";
		print OUT "$TargetSite" . "CGG\n";
		print OUT "$TargetSite" . "GGG\n";
		print OUT "$TargetSite" . "TGG\n";
		print OUT "$TargetSite" . "AGT\n";
		print OUT "$TargetSite" . "CGT\n";
		print OUT "$TargetSite" . "GGT\n";
		print OUT "$TargetSite" . "TGT\n";
	}
}

print "$Pos Positive and $Neg Negative target sites found.\n";
close (OUT) or die "ERROR in $0: Cannot close outputfile $OutputFile\n";
close (IN) or die "ERROR in $0: Cannot close inputfile $InputFile\n";
