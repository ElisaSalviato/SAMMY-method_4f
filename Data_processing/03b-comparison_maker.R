## LIBRARIES
require( "Rcpp" )
require( "data.table" )
require( "spp" ) 

## INPUT
args <- commandArgs( trailingOnly = TRUE ) 

ip_file <- args[ 1 ] # Path to the file should be used as IP for the analysis (file should be in .bam format)
input_file <- args[ 2 ] # Path to the file should be used as INPUT for the analysis (file should be in .bam format)
chromsizes_file <- args[ 3 ] # Path where the dimension of the chromosomes are specified (look for UCSC chromsize specifications)
output_file <- args[ 4 ] # Name of the output file (".wig" extension required)

mle_output_file <- gsub( ".wig", "_mle.wig", output_file )

## PROCESSING
### Import the data
ip    <- read.bam.tags( ip_file )
input <- read.bam.tags( input_file )

### Make the list of chromosomes
chromsizes <- fread( chromsizes_file, data.table = FALSE )
chrs <- as.character( chromsizes[ , 1 ] )

### Select the informative tags
ip_tags <- ip$tags
input_tags <- input$tags

### Remove the tags with anomalies
ip_rm <- remove.local.tag.anomalies( ip_tags[ chrs ] )
input_rm <- remove.local.tag.anomalies( input_tags[ chrs ] )

### Make the kernel smoothing
chip_smoothed <- get.smoothed.tag.density( ip_rm, input_rm,
                                          tag.shift = 0,
                                          background.density.scaling = TRUE )

### Make the smoothing with the loglikelihood function
chip_smoothed_mle <- get.smoothed.enrichment.mle( ip_rm, input_rm,
                                                 tag.shift = 0,
                                                 background.density.scaling = TRUE )

## OUTPUT
### Write the wig of the chip_smoothed
writewig( chip_smoothed,
         output_file,
         "Smoothed tag density difference for ChIP-seq data" )

### Write the wig of th the mle
writewig( chip_smoothed_mle,
         mle_output_file,
         "Smoothed tag density difference for ChIP-seq data" )
