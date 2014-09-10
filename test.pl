require "Damerau.pl";
$Sequence=                 "GTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC";
$ExpectedLeadingSequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC";
my %DamerauResults;
DetermineDamerauLevenshteinDistance($ExpectedLeadingSequence,$Sequence,\%DamerauResults);
print "Distance:" . $DamerauResults{'Distance'} . "\n";
print "Accurately Determined:" . $DamerauResults{'AccuratelyDetermined'} . "\n";
foreach my $Result (keys $DamerauResults{'Changes'}) {
	print "$Result:" . $DamerauResults{'Changes'}->{$Result} . "\n";
}
