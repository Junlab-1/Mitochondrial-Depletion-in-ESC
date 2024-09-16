#!/bin/bash
# conda activate R4
mkdir temp
for i in `ls *_1.fq.gz`
do
	echo $i
	file=$(echo "$i" | awk -F "_1.fq" '{print $1}')
	file_1=$"$file"_1.fq.gz
	file_2=$"$file"_2.fq.gz
	file_1_trim=$"$file"_1_trim.fq.gz
	file_2_trim=$"$file"_2_trim.fq.gz
	trim_galore --quality 30 --stringency 3 --length 20 --paired --phred33 --gzip ../RNAdata/$file_1 ../RNAdata/$file_2 --output_dir temp
	mv temp/"$file"_R1*.gz ../RNAdata/$file_1_trim
	mv temp/"$file"_R2*.gz ../RNAdata/$file_2_trim
	file_sam=$file.sam
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/human/human_index -1 ../RNAdata/$file_1_trim -2 ../RNAdata/$file_2_trim  -S $file_sam  -t
	file_bam=$file.bam
	#tranfer from sam to bam
	samtools sort -@ 8 -o $file_bam $file_sam
	mkdir result/$file
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/human/Homo_sapiens.GRCh38.110.gtf -o result/$file/output_merge.gtf $file_bam
done

~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i result/
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i result/
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i result/
