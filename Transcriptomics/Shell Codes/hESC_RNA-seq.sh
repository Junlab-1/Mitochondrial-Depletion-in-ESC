#!/bin/bash
# conda activate R4
mkdir temp human_rmNUMT_gtf
listH=("H272_1" "H272_2" "H272_3")
for i in "${listH[@]}"
do
	echo $i
	file_1_trim=$"$i"_1_trim.fq.gz
	file_2_trim=$"$i"_2_trim.fq.gz
	hfile_sam=$"$i"_h.sam
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/human_rmNUMT/human_nonumt_index -1 ../$file_1_trim -2 ../$file_2_trim  -S $hfile_sam  -t
	hfile_bam=$"$i"_h.bam
	#tranfer from sam to bam
	samtools sort -@ 8 -o $hfile_bam $hfile_sam
	human_folder=$"h""$i"_all
	mkdir human_rmNUMT_gtf/$human_folder
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/human/rmNUMTs/rmNUMTgtf/human_filt_nomunt_MTfix.gtf ./$hfile_bam -o ./human_rmNUMT_gtf/$human_folder/output_merge.gtf
done

mkdir all_genecount
cd all_genecount
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../human_rmNUMT_gtf -g human_control_gene_nonumt.csv -t human_control_gene_transcript_nonumt.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../human_rmNUMT_gtf -g human_control_tpm_nonumt.csv -t human_control_tpm_transcript_nonumt.csv
