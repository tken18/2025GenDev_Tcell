#! /bin/bash
#$ -S /bin/bash

#Identification of Super-Enhancers(SEs) by ROSE (RANK ORDERING OF SUPER-ENHANCERS) using a gff file of signals in H3K27ac ChIP-seq.


#Example
#fastq_R1_file_name=./SRR7305519_1.fastq
#sample_ID_number=SRR7305519
#sample_name_Tag=WT_DN3_H3K27Ac-ChIP
#bowtie2_mm9=/home/username/reference/index/mm9



#Quality check using FastQC (ver_0.12.1)
fastqc ${fastq_R1_file_name}


#Trimming using HOMER (v4.11)
homerTools trim -3 GATCGGAAGAGCACACGTCT -mis 1 -minMatchLength 4 -min 15 ${fastq_R1_file_name}


#Alignment using bowtie2 (v.2.3.5)
bowtie2 -p 6 -x bowtie2_mm9 ${fastq_R1_file_name}.trimmed -S ${fastq_R1_file_name}_adapter_trimmed.mm9_alignment.sam


#Make TagDirectory using HOMER (v4.11)
makeTagDirectory ${sample_name_TAG}_${sample_ID_number}_R1_Adtrimmed_mm9/ -genome mm9 -checkGC ${fastq_R1_file_name}_adapter_trimmed.mm9_alignment.sam -format sam



#Peak calling using HOMER (v4.11)
findPeaks ${sample_name_TAG}_${sample_ID_number}_R1_Adtrimmed_mm9/ -style histone -o auto -i ${Input_tag_directory}/


#Manual conversion of the peak file format to the gff file format.
#Identification of SE by ROSE (version0.1).

#Example
#gff_file=./OP9_DN3_H3K27Ac_SRR7305519_mm9_regions.gff
#taraget_sorted_bam_file=./OP9_DN3_K27ac_SRR7305519_mm9spc_sorted.bam
#output_directory_name=OP9_DN3_H3K27Ac_SRR7305519_mm9_SE_ROSE

python ROSE_main.py -g mm9 -i ${gff_file} -r ${taraget_sorted_bam_file} -o ${output_directory_name} -s 12500 -t 2500 