#!/bin/bash
# conda activate R4
mkdir temp orangutan_rmNUMT_gtf
listH=("O273_1" "O273_2" "O273_3")
for i in "${listH[@]}"
do
	echo $i
	file_1_trim=$"$i"_1_trim.fq.gz
	file_2_trim=$"$i"_2_trim.fq.gz
	ofile_sam=$"$i"_o.sam
	hisat2 -q --dta -p 8 -x ~/reference/hisat2_index/orangutan_rmNUMT/orangutan_nonumt_index -1 ../$file_1_trim -2 ../$file_2_trim  -S $ofile_sam  -t
	ofile_bam=$"$i"_o.bam
	#tranfer from sam to bam
	samtools sort -@ 8 -o $ofile_bam $ofile_sam
	orangutan_folder=$"o""$i"_all
	mkdir orangutan_rmNUMT_gtf/$orangutan_folder
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/orangutan/rmNUMTs/rmNUMTgtf/orangutan_filt_nonumt_MTfix.gtf ./$ofile_bam -o ./orangutan_rmNUMT_gtf/$orangutan_folder/output_merge.gtf
done

mkdir all_genecount
cd all_genecount
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../orangutan_rmNUMT_gtf -g orangutan_control_gene_nonumt.csv -t orangutan_control_gene_transcript_nonumt.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../orangutan_rmNUMT_gtf -g orangutan_control_tpm_nonumt.csv -t orangutan_control_tpm_transcript_nonumt.csv
