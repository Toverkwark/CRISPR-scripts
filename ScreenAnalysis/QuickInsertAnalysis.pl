#!/bin/perl
$Search="GGCTTTATATATCTTGTGGAAAGGACGAAACACCG";
my $File=open(IN,"small");
my $Line=<IN>;
my $Line=<IN>;
my $Line=<IN>;
my $Line=<IN>;
while (defined($Line=<IN>)) {
        $Line=<IN>;
        $Position=index($Line,"GGCTTTATATATCTTGTGGAAAGGACGAAACACCG");
        if($Position) {
                $Sequence=substr($Line,$Position+36,20);
                $Sequences{$Sequence}++;
        }

        $Line=<IN>;
        $Line=<IN>;
}
foreach $Sequence (sort {$Sequences{$a} <=> $Sequences{$b}} keys %Sequences) {
        print "$Sequence:\t" . $Sequences{$Sequence} . "\n";
}
