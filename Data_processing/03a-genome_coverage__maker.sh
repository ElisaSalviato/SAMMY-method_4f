# INPUT
file=$1 # Input file (.bam)
output_file=$2 # Path to the output file
file_type=$3 # Type of the output file (e.g. "bigwig" )
CORE=$4 # Number of cores will be used for the analysis
bin=$5 # Binsize (expressed in base pairs) for genome coverage calculation (e.g. 50)
notalignable=$6 # Path to the file of the genomic regions shouldn't be included in the analysis (it need the path to a bed file)
genomesize=$7 # Estimated genome size to performing the normalization 
normalization_method=$8 # Method selecte to perform the normalizazion (e.g. "RPKM") 
extend=$9 # Number expressing (in base pairs) how much a reads should be extended (e.g. 250) 

# PROCESSING
$bamCoverage -b $file \
	    -o $output_file \
	    -of $file_type \
	    -p $CORE \
            -bs $bin \
            -bl $notalignable \
            --effectiveGenomeSize $genomesize \
            --normalizeUsing $normalization_method \
	    --extendReads $extend
