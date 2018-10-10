# These are the directories to use.
cerevisiae=/fh/scratch/delete90/hahn_s/aerijman/genomes/cerevisiae/cerevisiae
pombe=/fh/scratch/delete90/hahn_s/aerijman/genomes/pombe/pombe
fastq_dir=$1
sam_dir=$2

# Save a sample_names.txt file with the prefix names of the samples. 
# Then make executables for each data file (files.sh).
for i in `cat ${fastq_dir}/sample_names.txt`
do
	cat  > ${i}_sample_align.sh << eof1
#!/bin/bash
ls ${fastq_dir}/${i}*R1*fastq > ${i}_fastq_R1.list
ct=0
for j in \`cat ${i}_fastq_R1.list\`
do
	k=\`echo \$j | sed 's/_R1_/_R2_/'\`
	if [ \$ct -lt 1 ];then
		bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12 -x $cerevisiae -1 \$j -2 \$k > ${sam_dir}/${i}_Sc.sam
	else
		bowtie2 --no-head --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12 -x $cerevisiae -1 \$j -2 \$k >> ${sam_dir}/${i}_Sc.sam
	fi
	ct=\`echo "\$ct + 1" | bc\`
done
eof1

	cat > ${i}_spike_align.sh << eof2
#!/bin/bash
ls ${fastq_dir}/${i}*R1*fastq > ${i}_fastq_R1.list
ct=0
for j in \`cat ${i}_fastq_R1.list\`
do
	k=\`echo \$j | sed 's/_R1_/_R2_/'\`
	if [ \$ct -lt 1 ];then
		bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12 -x $pombe -1 \$j -2 \$k > ${sam_dir}/${i}_pombe.sam
	else
		bowtie2 --no-head --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12 -x $pombe -1 \$j -2 \$k >> ${sam_dir}/${i}_pombe.sam
	fi
	ct=\`echo "\$ct + 1" | bc\`
done
eof2

	#echo "sbatch --mem=21780 -p restart -c 12 --job-name=${i}_sample --output=err_out/${pre}_sample.OUT.%J ${i}_sample_align.sh --requeue" >> jobs
	#echo "sbatch --mem=21780 -p restart -c 12 --job-name=${i}_spike  --output=err_out/${pre}_spike.OUT.%J  ${i}_spike_align.sh  --requeue" >> jobs
done


# sam to pairs.  I am using the read_samB.pl here
#mkdir $fastq_dir/pairs
#for i in `ls -lt $fastq_dir*sam | awk '{print $9}'`; do complete=$i; echo processing... $complete; short=`echo $i | awk -F "/" '{print $8}'`; perl read_samB.pl $complete | sort -k1,1 -k3n,3n -k4n,4n | awk '{if ($2 ~/^[0-9]+/ && $3 ~/^[0-9]+/ && $4 ~/^[0-9]+/ && $5 ~/^[0-9]+/ && $6 ~/^[0-9]+/) print $0}' > $fastaq_dir/pairs/${short:0:${#short}-4}.pairs; done

# Calculate the scaling factor
#for i in `cat sample_names.txt`; do scaling_factor=`cat ${fastaq_dir}/pairs/${i}_pombe.pairs | wc -l`; echo "$i $scaling_factor" >> scaling.dat; done

#mkdir wiggles

# I have changed the pairs_2_wig so it normalizes through the total counts of pombe
#a=( $(cat scaling.dat ) );for ((i=0; i<${#a[@]}; i+=2)); do filename=$fastaq_dir/pairs/${a[$i]}_Sc.pairs; scale=${a[$i+1]}; python pairs_2_wiggle.py $filename $scale; done

#for i in `ls *pombe.pairs`; do scale=( $( wc -l $i )); echo ${scale[0]} ${scale[1]:0:-12}; done > scale.dat
