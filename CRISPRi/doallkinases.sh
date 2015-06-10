#!/bin/sh
for i in `cat /home/NKI/b.evers/bastiaan/gRNAs/kinases/additional`;
do 
echo $i >/home/NKI/b.evers/bastiaan/gRNAs/kinases/$i.list
/home/NKI/b.evers/CRISPR-scripts/CRISPRi/RunOneKinase.sh /home/NKI/b.evers/bastiaan/gRNAs/kinases/$i
done
