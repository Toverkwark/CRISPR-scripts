use Getopt::Long;

print "Usage:perl $0 -input -output -records -offset\n-input\t\tName of inputfile\n-output\t\tName of outputfile. If none chosen, inputfile +.sub is used\n-records\tNumber of records to be printed\n-offset\t\tNumber of records to skip before printing starts\n";
GetOptions (   "input=s"    		=> \$InputFile,
		"output=s"		=> \$OutputFile,
		"records=i"		=> \$Records,
		"offset=i" 		=> \$Offset,
);

if(!$OutputFile)
{
	$OutputFile=$InputFile . ".sub";
}

open (INPUT, $InputFile) or die "Input file $InputFile is not accessible.\n";
open (OUTPUT, ">", $OutputFile) or die "Output file $OutputFile is not accessible.\n";

for ($Record=1;$Record<=$Offset;$Record++)
{
	for ($RecordLine=1;$RecordLine<=4;$RecordLine++)
	{
		$line=<INPUT>;
	}
}

for ($Record=1;$Record<=$Records;$Record++)
{
		$line=<INPUT>;
		print OUTPUT $line;		
		$line=<INPUT>;
		print OUTPUT $line;		
		$line=<INPUT>;
		print OUTPUT $line;		
		$line=<INPUT>;
		print OUTPUT $line;		
}

close (INPUT) or die "Could not close input file $InputFile.\n";
close (OUTPUT) or die "Could not close output file $OutputFile.\n";
