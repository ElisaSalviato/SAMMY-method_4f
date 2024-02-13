## LIBRARIES
require( "data.table" )
require( "rtracklayer" )

source( "Scripts/3.1-bw_binner.Rl" )

## FUNCTIONS

## INPUT
### External variables
args <- commandArgs( trailingOnly = TRUE )
## args <- c(
##     input_file = "/storage/data/FF/Cristiano/Pipeline2.0/Fibroblast/FlowCell_220919/FlowCell_220919_SAMMY-seq/03-Analysis/Single_tracks/COO2_untreated_10K_S2.bw",
##     condition = "untreated",
##     protocol = "4fhet",
##     biorep = "C002",
##     fraction = "S2",
##     binsize = 250000,
##     date = "20220725" 
## )

input_file <- args[ 1 ]
condition <- args[ 2 ]
protocol <- args[ 3 ]
biorep <- args[ 4 ]
fraction <- args[ 5 ]
binsize <- as.numeric( args[ 6 ] )
date <- args[ 7 ]
output_name <- paste0( paste( biorep, condition, protocol, fraction, sep = '_' ), "__date", date, "__binsize", binsize / 1000, "Kb" )

### Internal variables
genome <- args[ 8 ] # "mm9"
genome_bedfile <- args[ 9 ] # "Input/mm9_chromsizes.bed"
chrs_bedfile <- args[ 10 ] # "Input/mm9_only-canonical_chromsizes.bed"
## chrs_bedfile <- "Input/mm9_chromsizes.bed"

output_path <- args[ 11 ]

### Print imported used variables
print( paste0( "Processing file: ", input_file ) )
print( paste0( "Condition: ", condition ) )
print( paste0( "Protocol: ", protocol ) )
print( paste0( "Biological replicate: ", biorep ) )
print( paste0( "Fraction: ", fraction ) )
print( paste0( "Binsize used: ", binsize / 1000, "Kb" ) )
print( paste0( "Date: ", date ) )
print( paste0( "Operative name: ", output_name ) )
print( paste0( "Analysis based on: ", genome ) )

## PROCESSING
### Prepare analysis
#### Select chromosomes to use
gchr_selecteds <- import( chrs_bedfile, format = "bed" )
chrs <- as.character( unique( seqnames( gchr_selecteds ) ) )

#### Make output dir
output_dir <- paste0( output_path, '/', output_name )

print( paste0( "Creating ", output_dir ) )
system( command = paste0( "mkdir -p ", output_dir ) )

### Determine the bins to use for the genome
bin_list <- bins_calculator( genome_bedfile, binsize, genome )

### Import files and rebin them
bws <- import_and_rebin__bw( input_file, bin_list, genome, output_name )
bw <- bws[[ 1 ]]

## OUTPUT
for( chr in chrs ){
    
    output_file <- paste0( output_dir, '/', output_name, "__", chr, ".bw" )
    print( paste0( "Creating: ", output_file ) )
    
    chr_gr <- bw[ which( seqnames( bw ) == chr ), ]
    export( chr_gr, output_file, format = "bigWig" )
    
}
