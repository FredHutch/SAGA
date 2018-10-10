#!/bin/bash

gunzip *.gz
# get sample names
for s in `ls *fastq`; do echo ${s:0:-16}; done | sort | uniq  > sample_names.txt
cp /fh/scratch/delete90/hahn_s/aerijman/genomes/run_rmRep.sh .
this_dir=$(pwd)
bash run_rmRep.sh ${this_dir} ${this_dir}

for i in *align.sh
do
	bash ${i}
done

#rm *fastq *sh *list
