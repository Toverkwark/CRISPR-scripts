#!/bin/sh
#Make a list of genes (list) and run in through the refseq finder:
echo $1 >$1.list
perl ConvertGeneSymbolsToRefSeqIDs.pl -l $1.list

#Then, create a list of all possible CRISPRi guides:
perl ObtainCRISPRiGuidesForGene.pl -i $1.list.out 

#Determine which of these target unique genes in the genome:
perl IdentifyUniqueTargetSites.pl -i $1.list.out.protospacers 

#Determine off target effects:
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique -d 1 -l
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.1 -d 2
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.1.2 -d 3
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.1.2.3 -d 4

#Delete intermediate files
rm $1.list
rm $1.list.out
rm $1.list.out.protospacers
rm $1.list.out.protospacers.unique
rm $1.list.out.protospacers.unique.1
rm $1.list.out.protospacers.unique.1.2
rm $1.list.out.protospacers.unique.1.2.3

