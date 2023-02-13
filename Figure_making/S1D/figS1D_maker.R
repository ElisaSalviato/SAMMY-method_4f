## LIBRARIES
require( data.table )
require( reshape2 )
require( ggplot2 )

## INPUT
data_file <- "Input/reads_table.tsv"
conversion_table_file <- "Input/conversion-name_table.tsv"
    
out_filename <- "./Output/paper_samples_summary"
out_file <- paste0( out_filename, ".pdf" )

colors <- rev( c( "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525" ) ) ## greyscale

fraction_order <- c( "S2", "S2S", "S2L", "S3", "S4" )
exple_order <- c( "4f\nC002\nrep1", "4f\nC004\nrep1", "4f\nC004\nrep2", "4f10K\nC002\nrep1", "4f10K\nC004\nrep1", "4f10K\nC004\nrep2" )

## PROCESSING
### Import the conversion table
conversion_table <- fread( file = conversion_table_file, header = TRUE, data.table = FALSE )
names( conversion_table ) <- c( "or", "new" )

### Import the data form the file
pre_df <- fread( file = data_file, header = TRUE, data.table = FALSE )

### Change the names
or_names <- pre_df$Sample

new_names <- c()

for( name in or_names ){

    prefix <- gsub( "_L.*", '', name )
    suffix <- gsub( ".*_L", "_L", name )

    new_prefix <- unique( conversion_table[ grep( prefix, conversion_table$or ), "new" ] )
    new_names <- c( new_names, new_prefix )

    if( length( new_prefix ) > 1 ){

        print( name )
        print( unique( conversion_table[ grep( prefix, conversion_table$or ), "or" ] ) )
        print( new_prefix )
        print( '' )
    }
    
}

pre_df$Sample <- new_names

### Calculate the absolute number of reads lost in each passage ore used
trimmed <- pre_df$All - pre_df$After_trim
not_aligned <- pre_df$After_trim - pre_df$Aligned
duplicates <- pre_df$Duplicated
filtered <- pre_df$Aligned - pre_df$Duplicated - pre_df$After_filtering
used <- pre_df$After_filtering

### Create a database with the just calculated data
ng_df <- data.frame( file = pre_df$Sample, trimmed, not_aligned, duplicates, filtered, used )

### Arrange the database to make it usable by ggplot2
#### Orient correctly the dataframe
df <- melt( ng_df )
df$sample <- df$file
df$variable <- factor( df$variable, levels = rev( c( "used", "filtered", "duplicates", "not_aligned", "trimmed" ) ) )

#### Add batch specific variables
df$id <- df$sample

df$sample <- gsub( "_.*", '', df$id )
df$exp <-  gsub( ".*_", '', gsub( "_S.*", '', df$file ) )

df$fraction <- gsub( ".*_S", 'S', df$file )
df$fraction <- factor( df$fraction, fraction_order )

df$exple <- paste0( df$exp , '_', df$sample )
df$exple <- gsub( "[-,_]", '\n', df$exple )

#### Make the dataframe of relative values and join with absolute one
df$datatype <- "n_reads"

ratio_df <- df
ratio_df$datatype <- "ratio"
for( line in seq( dim( ratio_df )[ 1 ] ) ){

    file <-  as.character( ratio_df$file[ line ] )
    variable <- ratio_df[ line, "variable" ]
    tot <- sum( df[ which( as.character( df$file ) == file ), "value" ] )
    
    value <- ratio_df[ line, "value" ]
    ratio_df[ line, "value" ] <- value / tot

}
df <- rbind( df, ratio_df )

## OUTPUT
### Make plots
df$exple <- factor( df$exple, levels = exple_order )

plot <- ggplot( data = df, aes( x = fraction, y = value ) ) +
    geom_bar( aes( fill = variable  ),
             stat = "identity", position = "stack" ) +
    theme_classic() +
    scale_fill_manual( values = colors ) +
    facet_grid( datatype ~ exple, scales = "free", space = "free_x", switch = 'y' ) +
    ylab( '' ) +
    labs( fill = "Alignemnt step" ) +
    theme( strip.background = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_text( angle = 90 ) ) +
    ggtitle( "Untreated samples" )
ggsave( "./Output/suppl1d.pdf", width = 18 )

tofile_list <- lapply( unique( ng_df$file ), function( id ){

    tmp_df <- ng_df[ which( ng_df$file == id ), c( "trimmed", "not_aligned", "duplicates", "filtered", "used" ) ]

    sums_by_col <- colSums( tmp_df )

    
    return( sums_by_col )

})
tofile_df <- data.frame( do.call( rbind, tofile_list ) )
tofile_df$id <- unique( ng_df$file )
                                                      
fwrite( tofile_df, "Output/supp_table1.csv" )
