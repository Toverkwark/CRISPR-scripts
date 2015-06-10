#!/bin/sh
#Make a list of genes (list) and run through the ensembl ID finder:
#That finder bases its information on the appris database to find those transcripts that appear most relevant for every gene
perl ConvertGeneSymbolsToEnsemblIDs.pl -l $1.list

#Then, create a list of all possible CRISPRi guides:
#This script only selects guides that start with a G that are between 19 and 22 (constraints are set at beginning of script)
#or are 20nt long and start with any other letter. It also already filters out guides that have >=4Ts or
#guides that are present in more than one gene
perl ObtainCRISPRiGuidesForGene.pl -i $1.list.out 

#Determine which of these target unique genes in the genome:
perl IdentifyUniqueTargetSites.pl -i $1.list.out.protospacers 

#Restructure file
perl RestructureGuideFile.pl -i $1.list.out.protospacers.unique

#Determine off target effects:
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.restructured -d 1 -l
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.restructured.1 -d 2
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.restructured.1.2 -d 3
perl DetermineOffTargets.pl -i $1.list.out.protospacers.unique.restructured.1.2.3 -d 4

#Filter and select
#This script prioritizes guides that target between -25 and +100 of the TSS
#It then selects max 10 guides per transcript. A Todo here is to limit the total number of guides per gene.
#Also, somehow in this script guides are selected to have exactly 20nt.
perl FilterAndSelectGuides.pl -i $1.list.out.protospacers.unique.restructured.1.2.3.4

#Delete intermediate files
rm $1.list
rm $1.list.out
rm $1.list.nonmatched
rm $1.list.out.protospacers
rm $1.list.out.protospacers.unique
rm $1.list.out.protospacers.unique.restructured.1
rm $1.list.out.protospacers.unique.restructured.1.2
rm $1.list.out.protospacers.unique.restructured.1.2.3
rm $1.list.out.protospacers.unique.restructured
