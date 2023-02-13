## LIBRARIES
require( "data.table" )
require( "ggplot2" )
require( "ggrepel" )

## FUNCTIONS
corr_per_patient <- function( fulldf, genome_bedfile ){

    genome_db <- fread( file = genome_bedfile, header = FALSE, data.table = FALSE )
    genome_db$size <- genome_db[ , 3 ] - genome_db[ , 2 ]
    rownames( genome_db ) <- genome_db[ , 1 ]
    
    cols <- unique( fulldf$col )
    rows <- unique( fulldf$row )

    dflist <- lapply( cols, function( col ){

        row_dflist <- lapply( rows, function( row ){

            tmp_df <- fulldf[ which(
                fulldf$col == col &
                fulldf$row == row ), ]

            tmp_df$chr_size <- genome_db[ tmp_df$chr, "size" ]

            genome_size <- sum( tmp_df$chr_size )

            patient_corr <- sum( tmp_df$corr_val * tmp_df$chr_size ) / genome_size

            row_df <- cbind( row = row, col = col, corr_val = patient_corr )

            return( row_df )

        })
        cols_df <- do.call( rbind, row_dflist )

    })
    df <- as.data.frame( do.call( rbind, dflist ) )

    df$corr_val <- as.numeric( as.character( df$corr_val ) )

    return( df )
    
}


## INPUT
fulldf_file <- "correlation_summary.tsv"
hg38_bedfile <- "/path/to/chromosome_sizes.bed"

colors <- c( S4vsS2 = "#DC5F42", S3vsS2S = "#AC0727" )
lit_order <- c( "H3K36me3", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "LaminA/C", "LaminB1" )

protocol_order <- c( "3f", "4f10K", "4f" )

exp_dbfile <- "exp_db.tsv"

status_db <- c( "H3K36me3" = "open", "H3K4me1" = "open", "H3K4me3" = "open", "H3K27ac" = "open", "H3K27me3" = "closed", "H3K9me3" = "closed", "LaminA/C" = "closed", "LaminB1" = "closed" )

title <- gsub( ".tsv", '', gsub( ".*_", '', fulldf_file ) )

## PROCESSING
### Import database of correlation per chromosome
fulldf <- fread( file = fulldf_file, header = TRUE, data.table = FALSE )

### Summarize correlation per patient
df <- corr_per_patient( fulldf = fulldf, genome_bedfile = hg38_bedfile )

### Marks on columns and SAMMYs on rows
small_df <- df[ which(
    grepl( "f", df$row ) &
    !grepl( "f", df$col ) ), ]

exp_db <- fread( file = exp_dbfile, header = TRUE, data.table = FALSE )
rownames( exp_db ) <- exp_db$Name

small_df$protocol <- exp_db[ as.character( small_df$row ), "Protocol" ]
small_df$sample <- exp_db[ as.character( small_df$row ), "Sample" ]
small_df$fraction <- exp_db[ as.character( small_df$row ), "Fraction" ]

small_df$macrofraction <- small_df$fraction
small_df$macrofraction[ grepl( "S2.", small_df$macrofraction ) ] <- "S2"

small_df$status <- status_db[ as.character( small_df$col ) ]
small_df$fraction <- factor( small_df$fraction, levels = names( colors ) )

small_df$protocol <- factor( small_df$protocol, levels = protocol_order )

small_df$col <- factor( small_df$col, levels = lit_order )

## OUTPUT
### Heatmap
small_df$row <- factor( small_df$row, levels = rev( levels( small_df$row ) ) )

g1 <- ggplot( small_df, aes( y = row, x = col ) ) +
    geom_tile( aes( fill = corr_val ), colour = "white" ) +
    geom_text( aes( label = round( corr_val, 2 ) ), size = 7 ) +
    facet_grid( protocol ~ ., scales = "free", shrink = FALSE ) + 
    scale_fill_gradient2( high = "#d8b365", mid = "#f5f5f5", low = "#5ab4ac" ) + 
    theme_bw() +
    theme( text = element_text( size = 18 ), axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1 ) ) +
    ylab( "Fraction" ) +
    xlab( "Chromatin marks" ) +
    ggtitle( label = paste0( "Binsize: ", title ) ) 
ggsave( filename = "figS3B.pdf", width = 18, height = 15 )

