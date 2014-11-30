use warnings;
use strict;
use Getopt::Std;
my $ApprisFile = "../refseq/appris_data.principal.txt";

#Import script options
my %ScriptOptions;
getopt( 'l', \%ScriptOptions );
die
"ERROR in $0:This script convert an input list of gene symbols into a list of ensembl IDs according to the appris principal list and a list of non-matched genes. Usage: perl convert_genesymbols_to_EnsemblIDs.pl -l list-of-input-genesymbols\n"
  unless my $GeneSymbolList = $ScriptOptions{'l'};

#Open inputfile and read in all gene symbols, write them to output including the refseq IDs, or to non-matched when no refseq IDs are found
open (OUT, ">", $GeneSymbolList . '.out') or die "Could not open output file\n";
open (NONMATCHED, ">", $GeneSymbolList . '.nonmatched') or die "Could not open non-matched file\n";
open( IN, $GeneSymbolList ) or die "Could not open input file $GeneSymbolList\n";
while ( defined(my $GeneSymbol = <IN> )) {
	chomp($GeneSymbol);
	print "Found gene symbol $GeneSymbol, finding transcript IDs\n";
	my @Transcripts=`grep -P "$GeneSymbol\t" $ApprisFile`;
		if ( @Transcripts ) {
		foreach (@Transcripts) {
			my @TranscriptIDs=split(/\t/,$_);
			print "Refseq ID found for gene $GeneSymbol:" . $TranscriptIDs[2] . "\n";
			print OUT "$GeneSymbol\t" . $TranscriptIDs[2] . "\n";
		}
	}
	else {
		print "No RefSeq IDs found that match gene symbol $GeneSymbol";
		print NONMATCHED "$GeneSymbol\n";
	}
}
close(IN) or die "Could not close input file $GeneSymbolList\n";
close(OUT) or die "Could not close output file\n";
close(NONMATCHED) or die "Could not close non-matched file\n";
