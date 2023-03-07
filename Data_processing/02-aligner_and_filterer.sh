echo "Working example: sh 02-aligner_and_filterer.sh 4 hg38_UCSC_onlycanonical.fa.gz test.fq.gz $PWD/Alignment_test 1540 1"

# INPUT
## Arguments
CORE=$1 # Number of cores that will be used for the analysis
GENOME=$2 # BWA genome that will be used for the analysis
fastq=$3 # path to fastq that should be aligned
output_folder=$4 # Folder where the output files will be stored
filtering_code=$5 # Code that will be used by samtools -F command to filter reads (e.g 1540)
min_quality=$6 # Minimum quality to keep reads during filtering (e.g. 1)

## Variables
fastq_name=$( basename $fastq | cut -d '.' -f 1 )
date=$( date +%y-%m-%d )
sam_header="@RG\tID:${fastq_name}\tPL:unknown\tPU:${date}_${fastq_name}\tSM:${fastq_name}"

full_bam=$output_folder/01-${fastq_name}_full.bam
mrkdup_bam=$output_folder/02-${fastq_name}_mrkdup.bam
filtered_bam=$output_folder/03-${fastq_name}_filtered.bam

statistics_folder=${output_folder}/Statistics/

## Programs
bwa="path/to/bwa"
samtools="path/to/samtools"
picard="java -jar path/to/picard"


# PROCESSING
## Create folders
mkdir $output_folder $statistics_folder # Permanent
tmp_align=$( mktemp -d ${output_folder}/alignment_tempXXXXXX ) # Temporary

## 1
## Alignment
$bwa aln -k 2 -t $CORE $GENOME $fastq \
    | $bwa samse -r $sam_header $GENOME /dev/stdin $fastq  \
    | $samtools sort -@ $CORE -O bam -T ${tmp_align} /dev/stdin > $full_bam 

## Create index of full alignament
$samtools index $full_bam

## General samtools statistics
$samtools flagstat -@ $CORE $full_bam > $statistics_folder/$( basename $full_bam ).stats
echo $full_bam

## 2
## Mark duplicates using Picard
$picard MarkDuplicates I=$full_bam O=$mrkdup_bam M=$statistics_folder/$( basename $mrkdup_bam ).stats

## Create index of marked duplicated bam
$samtools index $mrkdup_bam

## 3
## Filtering
$samtools view -@ $CORE -F $filtering_code -b -q $min_quality $mrkdup_bam > $filtered_bam
echo $mrkdup_bam

### Create Index
$samtools index $filtered_bam

## General samtools statistics
$samtools flagstat -@ $CORE $filtered_bam > $statistics_folder/$( basename $filtered_bam ).stats
echo $filtered_bam

rm -r ${tmp_align}

# Working example: sh 02-aligner_and_filterer.sh 4 hg38_UCSC_onlycanonical.fa.gz test.fq.gz $PWD/Alignment_test 1540 1
