#!/bin/bash
# conda activate R4
mkdir human_all_gtf orangutan_all_gtf
listH=("HO_H7" "HO_H8" "HO_H10" "HO_O6" "HO_O12" "HO_O14")
for i in "${listH[@]}"
do
	echo $i
	hfile_sam=$"$i"_h.sam
    ofile_sam=$"$i"_o.sam
    hfile_bam=$"$i"_h.bam
    ofile_bam=$"$i"_o.bam
	human_folder=$"h""$i"_all
	orangutan_folder=$"o""$i"_all
	mkdir human_all_gtf/$human_folder orangutan_all_gtf/$orangutan_folder
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/human/Homo_sapiens.GRCh38.110.gtf ./$hfile_bam -o ./human_all_gtf/$human_folder/output_merge.gtf
	stringtie -e -B -p 8 -G ~/reference/FastaGTF/orangutan/Pongo_abelii.Susie_PABv2.111.gtf ./$ofile_bam -o ./orangutan_all_gtf/$orangutan_folder/output_merge.gtf
done
# statistic gene exp
mkdir all_genecount
cd all_genecount
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../human_all_gtf -g human_ho_gene.csv -t human_ho_gene_transcript.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../human_all_gtf -g human_ho_tpm.csv -t human_ho_tpm_transcript.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getFPKM.py -i ../human_all_gtf -g human_ho_fpkm.csv -t human_ho_fpkm_transcript.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/prepDE.py -i ../orangutan_all_gtf -g orangutan_ho_gene.csv -t orangutan_ho_gene_transcript.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getTPM.py -i ../orangutan_all_gtf -g orangutan_ho_tpm.csv -t orangutan_ho_tpm_transcript.csv
~/miniconda3/envs/py27/bin/python2.7 ~/software/newprepDE/getFPKM.py -i ../orangutan_all_gtf -g orangutan_ho_fpkm.csv -t orangutan_ho_fpkm_transcript.csv