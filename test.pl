#!/bin/perl
use threads;
$numberofjobs=64;
for ($i=1;$i<=$numberofjobs;$i++) {
	print "i=$i\n";
	$thr[$i]=threads->create('doewat',$i);
}
for ($i=1;$i<=$numberofjobs;$i++) {
	$thr[$i]->join();
}

sub doewat($) {
	my ($k)=@_;
	print "Starting job $k\n";
	for ($i=0;$i<100000000;$i++) {
	}
	print "Done with job $k\n";
}
