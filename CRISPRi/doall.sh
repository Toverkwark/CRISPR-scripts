#!/bin/sh
for i in `cat /home/NKI/b.evers/bastiaan/gRNAs/epigenetic/additional`;
do nice -n 20 /home/NKI/b.evers/CRISPR-scripts/CRISPRi/RunOneGene.sh $i
done
