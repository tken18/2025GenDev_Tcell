#! /bin/bash
#$ -S /bin/bash


#Example
#fastq_R1_file_name=./SRR307056_1.fastq
#sample_ID_number=SRR307056
#sample_name_Tag=DN3_E2A_ChIP
#bowtie2_mm9=/home/username/reference/index/mm9



#Quality check using FastQC (ver_0.12.1)
fastqc ${fastq_R1_file_name}


#Trimming using HOMER (v4.11)
homerTools trim -3 GATCGGAAGAGCACACGTCT -mis 1 -minMatchLength 4 -min 15 ${fastq_R1_file_name}


#Alignment using bowtie2 (v.2.3.5)
bowtie2 -p 6 -x bowtie2_mm9 ${fastq_R1_file_name}.trimmed -S ${fastq_R1_file_name}_adapter_trimmed.mm9_alignment.sam


#Make TagDirectory using HOMER (v4.11)
makeTagDirectory ${sample_name_TAG}_${sample_ID_number}_R1_Adtrimmed_mm9/ -genome mm9 -checkGC ${fastq_R1_file_name}_adapter_trimmed.mm9_alignment.sam -format sam


#Make a bedgraph file using HOMER (v4.11)
makeUCSCfile ${sample_name_TAG}_${sample_ID_number}_R1_Adtrimmed_mm9/ -o auto


#Peak calling (for histone marks) using HOMER (v4.11)
findPeaks ${sample_name_TAG}_${sample_ID_number}_R1_Adtrimmed_mm9/ -style histone -o auto -i ${Input_tag_directory}/


#Peak calling (for transcription factors) using HOMER (v4.11)
findPeaks ${Target_tag_directory}/ -style factor -o auto -i ${Input_tag_directory}/

