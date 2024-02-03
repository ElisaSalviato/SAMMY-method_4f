# SAMMY-method 4f: bioinformatics analysis and plots generation

In this guide it is explained which are the bioinformatic steps to perform analysis of SAMMY-seq (4f protocol) as in "Biochemical properties of chromatin domains define genome compartmentalization" paper.

In this repository you are going to find four main folders: 
* `Scripts`: where the scripts to perform the basic bioinformatics analysis starting from a fastq are stored
* `Input_examples`: where some input files for testing are stored (in this guide we are going to analyze them)
* `Output_examples`: where the results of the examples commands performed on `Input_examples` are stored 
* `Figure_making`: where the scripts to make the figures presented in the paper are stored.

In this tutorial we are assuming that you will store all your outputs in a folder named `Output`. If you use as input the data stored in `Input_examples` double check your results with the ones stored in `Output_examples` to be sure everything worked fine.

## Before run the scripts
Before run the script described in this tutorial you should edit the `Scripts/environment_varialbles.src` with the path of the programs you are going to use and then run the following command (assuming you are located in the same folder where this README.md is placed):
```
source Scripts/environment_variables.src
```

## 1) Pre processing
### 1.1) Trimming

**Aim.** The first step of analysis is the trimming of the raw fastqs, this is done to remove from the reads all the adaptors and ensure a good standard quality for the alignment.

**Software.** This operation is done through the Trimmomatic software and it has been implemented through the bash script `01-trimmer.sh`, which you can find in the `Data_processing` folder. For running this script the installation of Trimmomatic is required, we suggest to install, through conda software, the version 0.39. 

**Input.** The trimming script take as input:
1. the type of reads: it can be `single-end` or `paired-end`, according to the type of sequencing performed
2. the type of trimming to perform: it can be `SE` (single-end) or `PE` (paired-end), according to the type of input fastq provided
3. the number or cores that should be used for the analysis: an integer is required
4. the type of phred that will be used to asses the sequenced reads quality: it can be 'phred33' or 'phred64', depending on the Illumina pipeline used 
5. the name of the file where the log of Trimmomatic are to be saved
6. the name of the file where the statistics of Trimmomatic are to be saved
7. extra commands used by Trimmomatic (e.g. `'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15'`)
8. the path to the input file that should be trimmed
9. the name of the output trimmed file

**Output.** At the end of the step multiple output file will be created. One  containing the statistic of the trimming process (.stat), one where the log are stored (.log) and a fastq where the adapters from the original fastq have been removed plus a fastq containing the removed unpeared reads (only if two paired fastqs have been provided as input).

**Command.**

You can run the script with a command as follows:
```
sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:/path/to/TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz
```

For the sake of this tutorial we are going to analyze the data in `Input_examples` through the following bash commands:
```
sh Scripts/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.log \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20375_C2C12_2M_S2S_L003_R1_rand100000.fastq.gz \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz

sh Scripts/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.log \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20376_C2C12_2M_S2L_L002_R1_rand100000.fastq.gz \
	Input_examples/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz

sh Scripts/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.log \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20377_C2C12_2M_S3_L003_R1_rand100000.fastq.gz \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz

sh Scripts/01-trimmer.sh \
	single-end \
	SE \
	1 \
	phred33 \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.log \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.stat \
	'ILLUMINACLIP:./Input_samples/Trimmomatic/TruSeq3-SE.fa:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' \
	Input_examples/S20378_C2C12_2M_S4_L001_R1_rand100000.fastq.gz \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz
```

### 1.2) Alignment
**Aim.** Once completed the trimming, the fastqs have to be aligned and then filtered to remove PCR duplicates, multimapping reads, bad quality aligned reads etc. 

**Software.** To perform this step we provide the script `02-aligner_and_filterer.sh`, which you can find in the `Scripts` folder. To use this script, you need to install bwa (version 0.7.17-r1188), samtools (version 1.17.1) and picard (version 2.22) softwares.

**Input.** The alignment and filtering script takes as input:
1. the number of cores that should be used for the analysis: an integer is required
2. the BWA indexed reference genome that will be used for the analysis
3. the path to the fastq file that should be aligned
4. the path to the folder where the output files will be stored
5. the numeric code that will be used by samtools -F command to filter reads (e.g 1540)
6. the minimum required quality for keeping reads during the filtering process (e.g. 1)

**Output.** From this analysis several bam files will be produced as output (coupled with the correspondent .bai indexes):
1. the file ending by "_full.bam" is the bam file just aligned to the reference genome
2. the file ending by "_mrkdup.bam" is the bam file aligned to the reference genome whose duplicates have been marked
3. the file ending by "_filtered.bam" is the bam file aligned to the reference genome whose reads have been filtered according to the input parameters

**Command.**
You can run the script with a command as follows:
```
sh 02-aligner_and_filterer.sh 4 mm9_UCSC_onlycanonical.fa.gz trimmed.fq.gz ./Output 1540 1
```

In this tutorial we are going to use the previously trimmed fastqs as input and the indexed genome `./Input_examples/Reference/chr15.fa.gz`[^1] as reference for this alignment process:
```
sh Scripts/02-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/02-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/02-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1

sh Scripts/02-aligner_and_filterer.sh \
	1 \
	Input_examples/Reference/chr15.fa.gz \
	Output/01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed.fastq.gz \
	Output \
	1540 \
	1
```

## 2) Post processing
The aligned and filtered reads produced by the previous commands are then processed to perform the analysis shown in the paper.

### 2.1) Coverage/single tracks
**Aim.** One of the most relevant information we could extract from the aligment and filtering of the sequencing files is how the reads distribute over the genome. This provide an information about how much the SAMMY-seq fraction in analysis is represented in each genomic region. The signal represented by this analysis is not the same that you can directly infer from the alignment (bam) file, but it is normalized by RPKM. This information is foundamental to call compartments.

**Software.** This analysis is performed throught the bamCoverage tool of DeepTools software suite. To perform this step we provide the script `03a-genome_coverage__maker.sh` in `Scripts` folder.

**Input.** The script takes in input:
1. the bam file storing the information about the aligned and filtered reads created in the step 2) (e.g. input_filtered_aligned_trimmed.bam)
2. the output file where the reads will be stored (e.g. output.bw)
3. the output file type (e.g. bigwig)
4. the number of cores used to perform this process (e.g. 1)
5. the binsize for genome converage calculation expressed in base pairs (e.g. 5000)
6. the path to the sequences we want to exclude (e.g. blacklist.bed)
7. the genome effective size you can find on DeepTools *effective genome size* manual page (e.g. per mm9: 2620345972)
8. the method used for normalization (e.g. "RPKM")
9. the extension size of the aligned reads in base pairs (e.g. 250)

**Ouptut** The output is a bigwig file containing the reads distribution of each genomic region available in the input alignment file, normalized and without the blacklist regions.

**Command**
You can run the script as follow:
```
sh Scripts/03a-genome_coverage__maker.sh input_filtered_aligned_trimmed.bam output.bw bigwig 1 5000 blacklist.bed 2620345972 RPKM 250
```
To test the script here provided you can run one or all of the several commands:

```
sh Scripts/03a-genome_coverage__maker.sh \
	Output/02.3-01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S2S_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/03a-genome_coverage__maker.sh \
	Output/02.3-01-S20376_C2C12_2M_S2L_L002_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S2L_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/03a-genome_coverage__maker.sh \
	Output/02.3-01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S3_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250

sh Scripts/03a-genome_coverage__maker.sh \
	Output/02.3-01-S20378_C2C12_2M_S4_L001_R1_rand100000_trimmed_filtered.bam \
	Output/03a-S4_C2C12.bw \
	bigwig \
	1 \
	50 \
	Input_examples/Blacklist/mm9-blacklist.bed \
	2620345972 \
	RPKM \
	250
```

### 2.2) Comparisons/relative enrichement tracks
**Aim.** To increase the resolution of signal is often useful compare the signal coming from two different SAMMY-seq fractions, through this operation is possible better define regions of interests and check at the differences between two fractions. This operation is the same that is performed in ChIP-seq in performing a diffirence or a ratio between the signal coming from an IP and the one coming from an INPUT.

**Software.** To perform this operation we are using the Rscript `03b-comparison_maker.R` that you could find in the `Scripts` folder. For running this script `spp`, `data.table` and `Rcpp` libraries are required. We suggest to run the script using the R 3.5.2 version.

**Input.** For running this analysis/script four inputs are required:
1. the bam file of the first fraction we want compare (e.g. S2S.bam)
2. the bam file of the second fraction we want compare (e.g. S3.bam)
3. the chromosome size file, a file where for each lines there are three columns tab separated: (1) chromosome name as reported in the bams, (2) the coordinate of the first base of the chromosome (usually 1), (3) the coordinate of the last base of the chromosome. You can use as example the file `Input_examples/chr15.sizes`
4. the path and name.wig of the file where do you want to store your output (e.g. Output/S2SvsS3.wig)

**Ouput.** From running this script two (wiggle formatted) output files will results:
1. once ending with just ".wig" that is the output get.smoothed.tag.density(..., tag.shift = 0, background.density.scaling = TRUE) SPP function
2. once ending with just "_mle.wig" that is the output get.smoothed.enrichment.mle(..., tag.shift = 0, background.density.scaling = TRUE) SPP function

**Command.** 

You can run this script as follow:
```
Rscript 03b-comparison_maker.R S2S.bam S3.bam chr.sizes S2SvsS3.wig
```

If you followed the previous command of this tutorial you can test the programm running the following line:
```
Rscript Scripts/03b-comparison_maker.R Output/02.3-01-S20375_C2C12_2M_S2S_L003_R1_rand100000_trimmed_filtered.bam Output/02.3-01-S20377_C2C12_2M_S3_L003_R1_rand100000_trimmed_filtered.bam Input_examples/chr15.size Output/03b-S2SvsS3_C2C12.wig
```

## 3) Downstream analyses and data visualization
### 3.1) Compartments and sub-comparments analysis
### 3.2) Manuscript figures (?)


[^1]: This index has been produced through the command `bwa index chr15.fa.gz` performed on chromosome chr15 of mm9 genome, dowloaded from UCSC at the following [link](https://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz).
