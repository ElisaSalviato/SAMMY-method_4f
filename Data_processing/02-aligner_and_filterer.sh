# INPUT
## Arguments
CORE=${1} # Number of cores that will be used for the analysis
GENOME=${2} # BWA genome that will be used for the analysis
fastq=${3} # path to fastq that should be aligned
output_folder=${4} # Folder where the output files will be stored
filtering_code=${5} # Code that will be used by samtools -F command to filter reads (e.g 1540)
min_quality=${6} # Minimum quality to keep reads during filtering (e.g. 1)

## Variables
fastq_name=$( basename ${fastq} | cut -d '.' -f 1 )
date=$( date +%y-%m-%d )
sam_header="@RG\tID:${fastq_name}\tPL:unknown\tPU:${date}_${fastq_name}\tSM:${fastq_name}"

full_bam=${output_folder}/02.1-${fastq_name}_full.bam
mrkdup_bam=${output_folder}/02.2-${fastq_name}_mrkdup.bam
filtered_bam=${output_folder}/02.3-${fastq_name}_filtered.bam

statistics_folder=${output_folder}/02-Statistics/

## Programs
bwa="${bwa}"
samtools="${samtools}"
picard="${picard}"

# PROCESSING
# Display options
echo "Cores used for this analysis: ${CORE}"
echo "Reference genome used in the analysis: ${GENOME}"
echo "Input fastq that will be aligned: ${fastq}"
echo "Folder where the output will be stored: ${output_folder}"
echo "Samtools code to filtering reads: ${filtering_code}"
echo "Minimum mapping quality to keep aligned reads: ${min_quality}"
echo 

## Create folders
mkdir -p ${output_folder} ${statistics_folder} # Permanent
tmp_align=$( mktemp -d ${output_folder}/alignment_tempXXXXXX ) # Temporary

## 1
## Alignment
echo "Alignment"
${bwa} aln -k 2 -t ${CORE} ${GENOME} ${fastq} \
    | ${bwa} samse -r ${sam_header} ${GENOME} /dev/stdin ${fastq}  \
    | ${samtools} sort -@ ${CORE} -O bam -T ${tmp_align} /dev/stdin > ${full_bam} 
echo 

## Create index of full alignament
echo "Indexing aligned file"
${samtools} index ${full_bam}
echo

## General samtools statistics
echo "Making statistics on the aligned file"
${samtools} flagstat -@ ${CORE} ${full_bam} > ${statistics_folder}/$( basename ${full_bam} ).stats
echo 

echo "Aligned bam produced: ${full_bam}"
echo
echo

## 2
## Mark duplicates using Picard
echo "Marking duplicates on aligned bam"
${picard} MarkDuplicates I=${full_bam} O=${mrkdup_bam} M=${statistics_folder}/$( basename ${mrkdup_bam} ).stats
echo

## Create index of marked duplicated bam
echo "Indexing bam file with duplicates marked"
${samtools} index ${mrkdup_bam}
echo

echo "Bam with duplicates marked: ${mrkdup_bam}"
echo
echo

## 3
## Filtering
echo "Filtering reads from marked duplicates bam"
${samtools} view -@ ${CORE} -F ${filtering_code} -b -q ${min_quality} ${mrkdup_bam} > ${filtered_bam}
echo ${mrkdup_bam}

### Create Index
echo "Indexing reads of bam after reads filtered"
${samtools} index ${filtered_bam}
echo

## General samtools statistics
echo "Statistics on bam after reads filtered"
${samtools} flagstat -@ ${CORE} ${filtered_bam} > ${statistics_folder}/$( basename ${filtered_bam} ).stats
echo

echo "Aligned bam filtered: ${filtered_bam}"

rm -r ${tmp_align}
