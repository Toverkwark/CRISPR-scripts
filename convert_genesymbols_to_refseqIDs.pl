use warnings;
use strict;
use Getopt::Std;
my $RefSeqFile = "refseq/hg19.txt";

#Import script options
my %ScriptOptions;
getopt( 'l', \%ScriptOptions );
die
"ERROR in $0:This script convert an input list of gene symbols into a list of RefSeq IDs and a list of non-matched genes. Usage: perl convert_genesymbols_to_refseqIDs.pl -l list-of-input-genesymbols\n"
  unless my $GeneSymbolList = $ScriptOptions{'l'};

#Open inputfile and read in all gene symbols, write them to output including the refseq IDs, or to non-matched when no refseq IDs are found
open (OUT, ">", $GeneSymbolList . '.out') or die "Could not open output file\n";
open (NONMATCHED, ">", $GeneSymbolList . '.nonmatched') or die "Could not open non-matched file\n";
open( IN, $GeneSymbolList ) or die "Could not open input file $GeneSymbolList\n";
while ( defined(my $GeneSymbol = <IN> )) {
	chomp($GeneSymbol);
	print "Found gene symbol $GeneSymbol, finding transcript IDs\n";
	open( REFSEQ, $RefSeqFile ) or die "ERROR in $0:Cannot open RefSeq ID file $RefSeqFile\n";
	my $GeneFound = 0;
	while ( defined( my $Line = <REFSEQ> ) ) {
		chomp($Line);
		my @RefSeqValues = split( /\t/, $Line );
		my $Gene = $RefSeqValues[12];
		if ( $Gene eq $GeneSymbol ) {
			$GeneFound++;
			print "Refseq ID found for gene $GeneSymbol:" . $RefSeqValues[1] . "\n";
			print OUT "$GeneSymbol\t" . $RefSeqValues[1] . "\n";
		}
	}
	if ( $GeneFound == 0 ) {
		print "No RefSeq IDs found that match gene symbol $GeneSymbol";
		print NONMATCHED "$GeneSymbol\n";
	}
}
close(IN) or die "Could not close input file $GeneSymbolList\n";
close(OUT) or die "Could not close output file\n";
close(NONMATCHED) or die "Could not close non-matched file\n";