#!/bin/bash

#folders=$(ls -d /shared/ngs/illumina/lwarfiel/*/Unaligned/Project_lwarfiel/*Taf*/)
# folders should be splited on "/" and the 9 element is the actual folder.

# should run in parallel, for each of the folders 
# 1. make folder
# 2. copy the gz files 
# 3. extract
# 4. map to genome
# 5. delete fastq and keep sam
# 6. make wig and metgen

folder_name=${2}_$(echo $1 | awk 'BEGIN{FS="/"}{print $9}')
mkdir ${folder_name}
cp ${1}/*fastq.gz ./${folder_name}
cp ./second.sh ./${folder_name}

cd ${folder_name}
bash second.sh

