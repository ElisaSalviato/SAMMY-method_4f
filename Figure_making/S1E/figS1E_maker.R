## LIBRARIES
require( "data.table" )
require( "ggplot2" )
require( "ggrepel" )

## FUNCTIONS
unalignable_bias <- function( bedfile ){

    bed_df <- fread( file = bedfile, header = FALSE, data.table = FALSE )
    names( bed_df ) <- c( "chr", "start", "end" )

    chrs <- unique( bed_df$chr )

    bias_df <- data.frame() 
    for( chr in chrs ){

        tmp_df <- bed_df[ which( bed_df$chr == chr ), ]
        bias <- sum( tmp_df$end - tmp_df$start )

        bias_row <- cbind( chr = chr, bias = bias )
        bias_df <- rbind( bias_df, bias_row )
    }
    bias_df$bias <- as.numeric( as.character( bias_df$bias ) )

    return( bias_df )
    
}

## INPUT
table_dir <- "./Input/" ## Folder where is stored a single file for each fraction of each SAMMY-seq experiment that contain the statistics of the alignement calculated with Samtools (v1.6.1)
files <- dir( table_dir, full.names = TRUE )

unalignable_file <- "/path/to/blacklist_regions.bed"

## PROCESSING
system( "mkdir -p ./Output/" )

colors <- c( S2 = "#DC5F42", S2S = "#AC0727", S2L = "#D8553E", S3 = "#89BEDA", S4 = "#5E8BC0" )

### Calculate unalignable bias
unalignable_bias_db <- unalignable_bias( unalignable_file )

### Import statistics calculated by Samtools coverage
dflist <- lapply( files, function( file ){

    df <- fread( file = file, header = TRUE, data.table = FALSE )

    id <- gsub( "\\..*", '', basename( file ) )
    sample_prot <- gsub( "_S.*", '', id )

    sample <- gsub( "_.*", '', sample_prot )
    protocol <- gsub( ".*_", '', sample_prot )
    fraction <- gsub( ".*_", '', id )

    df$id <- id
    df$sample <- sample
    df$protocol <- protocol
    df$fraction <- fraction

    return( df )

})
df <- do.call( "rbind", dflist )

df$fraction <- factor( df$fraction, levels = names( colors ) )
df$protocol <- factor( df$protocol, levels = c( "3f", "4f10K", "4f" ) )

## OUTPUT
### Percentage of coverage per chromosome
coverage_plot <- ggplot( df, aes( x = fraction, y = coverage, fill = fraction ) ) +
    geom_boxplot( width = 0.2 ) +
    geom_text_repel( aes( label = `#rname` ), fill = "white" ) + 
    facet_grid( protocol ~ sample ) +
    scale_fill_manual( values = colors, breaks = c( "S2", "S2S", "S2L", "S3", "S4" ) ) +
    ggtitle( label = "Percentage of covered bases [0..100] per chromosome" ) +
    theme_bw()
ggsave( "./Output/figS1E.pdf", width = 12, height = 9 )
