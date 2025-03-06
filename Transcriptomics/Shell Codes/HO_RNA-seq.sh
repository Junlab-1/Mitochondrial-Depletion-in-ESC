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
        hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/human_rmNUMT/human_nonumt_index -1 ../ho/$file_1_trim -2 ../ho/$file_2_trim  -S $hfile_sam  -t
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/orangutan_rmNUMT/orangutan_nonumt_index -1 ../ho/$file_1_trim -2 ../ho/$file_2_trim  -S $ofile_sam  -t \
	hfile_bam=$"$i"_h.bam
	ofile_bam=$"$i"_o.bam
	#tranfer from sam to bam
	samtools sort -@ 8 -o $hfile_bam $hfile_sam
	samtools sort -@ 8 -o $ofile_bam $ofile_sam
	human_folder=$"h""$i"_all
	orangutan_folder=$"o""$i"_all
	mkdir human_all_gtf/$human_folder orangutan_all_gtf/$orangutan_folder
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/human/rmNUMTs/rmNUMTgtf/human_filt_nomunt_MTfix.gtf ./$hfile_bam -o ./human_all_gtf/$human_folder/output_merge.gtf
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/orangutan/rmNUMTs/rmNUMTgtf/orangutan_filt_nonumt_MTfix.gtf ./$ofile_bam -o ./orangutan_all_gtf/$orangutan_folder/output_merge.gtf
done

mkdir all_genecount
cd all_genecount
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../human_all_gtf -g human_ho_gene_nonumt.csv -t human_ho_gene_transcript_nonumt.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../orangutan_all_gtf -g orangutan_ho_gene_nonumt.csv -t orangutan_ho_gene_transcript_nonumt.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../human_all_gtf -g human_ho_tpm_nonumt.csv -t human_ho_tpm_transcript_nonumt.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../orangutan_all_gtf -g orangutan_ho_tpm_nonumt.csv -t orangutan_ho_tpm_transcript_nonumt.csv


