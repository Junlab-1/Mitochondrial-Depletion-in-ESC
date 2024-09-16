#!/bin/bash
# conda activate R4
mkdir temp result
listH=("HO_H7" "HO_H8" "HO_H10" "HO_O6" "HO_O12" "HO_O14")
for i in "${listH[@]}"
do
	echo $i
	#$RNA_singularity fastq-dump --split-3 --gzip $i
	file_1=$"$i"_1.fq.gz
	file_2=$"$i"_2.fq.gz
	file_1_trim=$"$i"_1_trim.fq.gz
	file_2_trim=$"$i"_2_trim.fq.gz
	trim_galore --quality 30 --stringency 3 --length 20 --paired --phred33 --gzip ~/PosdocProject/Denial/orangutan-human/clean/$i/$file_1 ~/PosdocProject/Denial/orangutan-human/clean/$i/$file_2 --output_dir temp
	mv temp/"$i"_1*.gz ./$file_1_trim
	mv temp/"$i"_2*.gz ./$file_2_trim
	hfile_sam=$"$i"_h.sam
	ofile_sam=$"$i"_o.sam
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/human/human_index -1 ./$file_1_trim -2 ./$file_2_trim  -S $hfile_sam  -t
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/orangutan/orangutan_index -1 ./$file_1_trim -2 ./$file_2_trim  -S $ofile_sam  -t
	hfile_bam=$"$i"_h.bam
	ofile_bam=$"$i"_o.bam
	#tranfer from sam to bam
	samtools sort -@ 8 -o $hfile_bam $hfile_sam
	samtools sort -@ 8 -o $ofile_bam $ofile_sam
done
