#!/bin/perl
$Search="GGCTTTATATATCTTGTGGAAAGGACGAAACACCG";
$Library="/home/NFS/users/b.evers/libraries/crispr.epigenetics.library.csv";
my $File=open(IN,"small");
my $Line=<IN>;
my $Line=<IN>;
my $Line=<IN>;
my $Line=<IN>;
while (defined($Line=<IN>)) {
	$Line=<IN>;
	$Position=index($Line,"GGCTTTATATATCTTGTGGAAAGGACGAAACACCG");
	if($Position) {
		$Sequence=substr($Line,$Position+35,20);
		$Sequences{$Sequence}++;
	}

	$Line=<IN>;
	$Line=<IN>;
}
$Total=0;
$InLibraryTotal=0;
foreach $Sequence (sort {$Sequences{$a} <=> $Sequences{$b}} keys %Sequences) {
	$InLibrary=`grep -P "$Sequence" $Library`;
	print "$Sequence:\t" . $Sequences{$Sequence} . "\t";
	if ($InLibrary) {
		print "yes\n";
		$InLibraryTotal=$InLibraryTotal+$Sequences{$Sequence};
	}
	else {
		print "no\n";
	}
	

	$Total=$Total+$Sequences{$Sequence};
}
print "In Library:$InLibraryTotal\n";
print "Total:$Total\n";

