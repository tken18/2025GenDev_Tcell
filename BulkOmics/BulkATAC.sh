#! /bin/bash
#$ -S /bin/bash


#Example
#fastq_R1_file_name=./KM7801_S73_R1_001.fastq
#fastq_R2_file_name=./KM7801_S73_R2_001.fastq
#sample_ID_number=KM7801
#sample_name_Tag=ATAC-B6-ETP-1
#bowtie2_mm9=/home/username/reference/index/mm9 





#Quality check using FastQC (ver_0.12.1)
fastqc ${fastq_R1_file_name}
fastqc ${fastq_R2_file_name}


#Trimming using Skewer (ver_0.2.2)
skewer -x CTGTCTCTTATA -y CTGTCTCTTATA -m pe -q 20 -l 30 -t 6 -o ${sample_ID_number} ${fastq_R1_file_name} ${fastq_R2_file_name}



#Alignment using bowtie2 (v.2.3.5)
bowtie2 -p 6 -x bowtie2_mm9 -1 ${sample_ID_number}-trimmed-pair1.fastq -2 ${sample_ID_number}-trimmed-pair2.fastq -S ${KM_sample_ID_number}_PE_skewer_mm9_alignment.sam


#Remove of mitochondrial reads from an sam file
sed '/chrM/d' ${sample_ID_number}_PE_skewer_mm9sc_alignment.sam >${sample_ID_number}_PE_skewer_mm9_alignment_wo_chrM.sam


#Make TagDirectory using HOMER (v4.11)
makeTagDirectory ${sample_name_TAG}_${sample_ID_number}_sk_wo_chrM_mm9/ -genome mm9 -checkGC ${sample_ID_number}_PE_skewer_mm9_alignment_wo_chrM.sam


#Make a bedgraph file using HOMER (v4.11)
makeUCSCfile ${TagDirectory_name}/ -o auto -style dnase


#Peak calling using HOMER (v4.11)
findPeaks ${sample_name_TAG}_${sample_ID_number}_sk_wo_chrM_mm9/ -o auto -style dnase >ATAC_peak_file.txt


#Annotate a peak file using HOMER (v4.11)
annotatePeaks.pl ATAC_peak_file.txt mm9 >output.txt


#Finding overlapping peaks using HOMER (v4.11)
mergePeaks -d given ATAC_peak_file.txt ATAC_peak_file-2.txt >ATAC_newPeakFile.txt


#Finding differentially bound peaks (4-fold more tags and a cumulative Poisson p-value less than 0.0001) using HOMER (v4.11)
getDifferentialPeaks ${ATAC_peak_file} ${target_TagDirectory} ${background_TagDirectory} >output.txt


#Counting the number of tags at each peak position from different sequencing experiments and output of rlog variance stabilized transformed data using HOMER (v4.11)
annotatePeaks.pl ${ATAC_peak_file} mm9 -size 300 -rlog -d ${TagDirectory_name1}/ ${TagDirectory_name2}/ ${TagDirectory_name3}/ ${TagDirectory_name4}/ >output_rlog.txt


