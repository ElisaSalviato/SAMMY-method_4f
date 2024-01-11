# SAMMY-method 4f: bioinformatics analysis and plot generation

In this guide it is explained which are the bioinformatic steps to perform analysis of SAMMY-seq (4f protocol) as in "Biochemical properties of chromatin domains define genome compartmentalization" paper.

In this repository you are going to find three main folders: 
. `Data_processing`: where the scripts to perform the basic bioinformatics analysis starting from a fastq are stored,
. `Input_samples`: where some input files for testing are stored (in this guide we are going to analyze them)
. `Figure_making`: where the scripts to make the figures presented in the paper are stored
. `Analyzed`: where the results of the examples commands performed on `Input_samples` are stored.


## 1) Trimming

The first step of analysis is the trimming of the raw fastqs.
This operation is done through the Trimmomatic software and it has been implemented through the bash script `01-trimmer.sh` in the folder `Data_processing`. For running this script the installation of Trimmomatic is required, we suggest to install through conda the version 0.39. 

The script take as input:

1. the type of reads: it can be or `single-end` or `paired-end` according to the type of sequencing performed
2. the type of trimming to perform (according to the input fastq provided): it could be `SE` or `PE` according if we are going to trim Single-End or Paired-End reads
3. the number or cores that should be used for the analysis (an integer is required)
4. the type of phred will be used for asses the sequenced reads quality
5. the file where the log of Trimmomatic are saved
6. the file where the statistics of Trimmomatic will be saved
7. extra command used by Trimmomatic (e.g. `'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15'`)
8. the input file you want trim
9. the name of the output file

You can run the script with a command as follow:
```
sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz
```

For the sake of this tutorial we are going to analyze the data in `Input_samples` through the following bash commands:
```
sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-C004-r2_S2S_L001_R1_rand100000__raw_trimmed.fq.gz.log \
	Analyzed/01-C004-r2_S2S_L001_R1_rand100000__raw_trimmed.fq.gz.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/C004-r2_S2S_L001_R1_rand100000__raw.fq.gz \
	Analyzed/01-C004-r2_S2S_L001_R1_rand100000__raw_trimmed.fq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-C004-r2_S2L_L001_R1_rand100000__raw_trimmed.fq.gz.log \
	Analyzed/01-C004-r2_S2L_L001_R1_rand100000__raw_trimmed.fq.gz.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/C004-r2_S2L_L001_R1_rand100000__raw.fq.gz \
	Analyzed/01-C004-r2_S2L_L001_R1_rand100000__raw_trimmed.fq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-C004-r2_S3_L001_R1_rand100000__raw_trimmed.fq.gz.log \
	Analyzed/01-C004-r2_S3_L001_R1_rand100000__raw_trimmed.fq.gz.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/C004-r2_S3_L001_R1_rand100000__raw.fq.gz \
	Analyzed/01-C004-r2_S3_L001_R1_rand100000__raw_trimmed.fq.gz

sh Data_processing/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Analyzed/01-C004-r2_S4_L001_R1_rand100000__raw_trimmed.fq.gz.log \
	Analyzed/01-C004-r2_S4_L001_R1_rand100000__raw_trimmed.fq.gz.stat \
	'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_samples/C004-r2_S4_L001_R1_rand100000__raw.fq.gz \
	Analyzed/01-C004-r2_S4_L001_R1_rand100000__raw_trimmed.fq.gz
```

## 2) Alignement
Once completed the trimming, the next step is to align the the fastqs and filter the reads removing: PCR duplicates, multimapping reads, bad quality aligned reads etc.

To perform this step we provided the script `02-aligner_and_filterer.sh` in folder `Data_processing`. To use this script you need to install bwa, samtools and picard softwares.
