#!/bin/perl
require "Damerau.pl";
$A="TTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTCC";
$B="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTC";
my %DamerauResults;
DetermineDamerauLevenshteinDistance($B,$A,\%DamerauResults);
print "AccuratelyDetermined:" . $DamerauResults{'AccuratelyDetermined'} . "\n";
print "Distance:" . $DamerauResults{'Distance'} . "\n";

foreach my $Change (keys $DamerauResults{'Changes'}) {
	$LeadingErrors{$Change}->{$DamerauResults{'Changes'}->{$Change}}++;
	print $DamerauResults{'Changes'}->{$Change} . " @" . $Change . " - ";
}
print "\n";
