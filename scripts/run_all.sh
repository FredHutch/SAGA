#!/bin/bash

folders=$(cat *txt)
#$(ls -d /shared/ngs/illumina/lwarfiel/*/Unaligned/Project_lwarfiel/*Taf*/)
ct=0
for f in ${folders}
do 
	ct=$(echo ${ct} + 1 | bc)
	echo "${ct} --> ${f}" >> files.txt
	sbatch -p campus -c 4 --job-name="${ct}" --output=${ct}.output --wrap="bash run_one.sh ${f} ${ct}"
done
