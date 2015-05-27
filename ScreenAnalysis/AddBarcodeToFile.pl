open(IN,"3500_1_pGIN_Screen_ATCACG_L001_R2_001.fastq");
open(OUT,">","temp");
while (defined(my $Line=<IN>)) {
	print OUT $Line;
	$Line=<IN>;
	$Line = "CGTGAT" . $Line;
	print OUT $Line;
	$Line=<IN>;
	print OUT $Line;
	$Line=<IN>;
	print OUT $Line;
}

close(IN);
close(OUT);
