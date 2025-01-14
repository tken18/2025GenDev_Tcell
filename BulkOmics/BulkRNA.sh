
#! /bin/bash
#$ -S /bin/bash


#Example
#fastq_R1_file_name=./KM9801_S17_R1_001.fastq
#fastq_R2_file_name=./KM9801_S17_R2_001.fastq
#sample_ID_number=KM9801
#sample_name_Tag=Cont-ETP-1
#STAR_mm9=/home/username/reference/STAR_Index/STAR_mm9



#Quality check using FastQC (ver_0.12.1)
fastqc ${fastq_R1_file_name}
fastqc ${fastq_R2_file_name}


#Alignment using STAR (ver_2.7.3a)
STAR --genomeDir STAR_mm9 --runThreadN 6 --readFilesIn ${fastq_R1_file_name} ${fastq_R2_file_name} --outFileNamePrefix STAR_mm9_${sample_ID_number}_R1R2_


#Makign TagDirectory using HOMER (v4.11)
makeTagDirectory ${sample_name_TAG}_${sample_ID_number}_R1R2_STAR_mm9/ -genome mm9 -checkGC STAR_mm9_${sample_ID_number}_R1R2_Aligned.out.sam -format sam 


#Making a bedgraph file using HOMER (v4.11)
makeUCSCfile ${sample_name_TAG}_${sample_ID_number}_R1R2_STAR_mm9/ -fragLength given -o auto


#Making a Transcripts Per Million table file using HOMER 
analyzeRepeats.pl rna mm9 -strand both -count exons -condenseGenes -d Exp1r1/ Exp1r2/ Exp2r1/ Exp2r2/ Exp3r1/ Exp3r2 -tpm > tpm.txt


#Making a raw count table file using HOMER (v4.11)
analyzeRepeats.pl rna mm9 -strand both -count exons -condenseGenes -d Exp1r1/ Exp1r2 Exp2r1/ Exp2r2/ Exp3r1/ Exp3r2 -noadj > raw.txt


#Differential Expression Analysis using HOMER(v4.11)
getDiffExpression.pl raw.txt cond1 cond1 cond2 cond2 > diffExp.output.txt