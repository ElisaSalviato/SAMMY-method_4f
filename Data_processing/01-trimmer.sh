echo "Working example with SE option: sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test 'ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15' input.fq.gz output.fq.gz"

# INPUT
## Common variables for analysis required both for "single-end" and "paired-end" input fastqs
reads_type=$1 # Argument that can be or "single-end" "paired-end" according to the type of sequencing performed
trimming_type=$2 # Method of trimming used by Trimmomatic ("SE" or "PE")
CORE=$3 # Number of cores that will be used for the analysis
phred=$4 # Type of phred will be used for asses the sequenced reads quality
trimlog_file=$5 # File where the log of Trimmomatic are saved
trimstat_file=$6 # File where the statistics of Trimmomatic are saved
extra_commands=$7 # String of extra command used by Trimmomatic (e.g. "ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15")

## Variables according to the type of fastq(s) in analysis
if [ $reads_type == "single-end" ]
then

    fastq=$8 # Path to single end fastq to trim
    fastq_output=$9 # Ouput file (fastq)
    
    
elif [ $reads_type == "paired-end" ]
then

    fastq_input1=$8 # Path to one of the two paired fastq to trim
    fastq_input2=$9 # Path to the other of the two paired fastq to trim
    fastq_output1=$10 # Output file after trimming $fastq_input1
    fastq_unpairedoutput1=$11 # Reads unpaired of the $fastq_input1
    fastq_output2=$12 # Output file after trimming $fastq_input2
    fastq_unpairedoutput2=$13 # Reads unpaired of the $fastq_input2
    
fi

## Program used by this script
trimmomatic="trimmomatic" # Path to trimmomatic, if you are using the java version write "java -jar path/to/trimmomatic"

# PROCESSING
if [ $reads_type == "single-end" ]
then

    $trimmomatic $trimming_type \
		 -threads $CORE \
		 -$phred \
		 -trimlog $trimlog_file \
		 -summary $trimstat_file \
		 $fastq \
		 $fastq_output \
		 $extra_commands

elif [ $reads_type == "paired-end" ]
then
    
    $trimmomatic $trimming_type \
		 -threads $CORE \
		 -$phred \
		 -trimlog $trimlog_file \
		 -summary $trimstat_file \
		 $fastq_input1 $fastq_input2 \
		 $fastq_output1 $fastq_unpairedoutput1 \
		 $fastq_output2 $fastq_unpairedoutput2 \
		 $extra_commands

else

    echo "Error: specify if are you analysing sequencing based on 'single-end' or 'paired-end'"
    echo "Warning: In case of 'paired-end' the two strand should stay in two different files"
    
fi

# Working example with SE option: sh 01-trimmer.sh single-end SE 1 phred33 log_test stat_test "ILLUMINACLIP:TruSeq3-SE:2:30:10 MINLEN:35 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15" input.fq.gz output.fq.gz
