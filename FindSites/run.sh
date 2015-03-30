perl CreateRandomTargetList.pl -output=foo
perl AddNGGToTargetList.pl -input=foo
/media/Data/iKRUNC/bowtie2-2.1.0/bowtie2 /media/Data/iKRUNC/hg19-index/hg19 -r foo.potentialtargets -t -S foo.potentialtargets.matched --no-hd --score-min L,-5,0 -p 64
perl CountUniqueTargetOccurence.pl -input=foo.potentialtargets.matched
perl QualityScoreSites.pl -i foo.potentialtargets.matched.out -o foo.qualities.1 -d 1 -s 20 -l
perl QualityScoreSites.pl -i foo.qualities.1 -o foo.qualities.2 -d 2 -s 20 
perl QualityScoreSites.pl -i foo.qualities.2 -o foo.qualities.3 -d 3 -s 20 
perl QualityScoreSites.pl -i foo.qualities.3 -o foo.qualities.4 -d 4 -s 20 


