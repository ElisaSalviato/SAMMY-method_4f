## IMPORT ARGS OR EXIT
args <- commandArgs( trailingOnly = TRUE )
if( is.na( args[ 1 ] ) ){

    print( "Arguments to run the script:" )
    print( " 1) METAFILE" )
    print( " 2) PLOTNAME" )
    print( " 3) CORES" )
    print( " 4) AUX FOLDER (optional)" )
    print( "Example: Rscript automatic_gViz.R metafile.tsv myplot 12 AUX_gViz" )
    quit()

}

## LIBRARIES
suppressPackageStartupMessages({

    ### External libraries
    require( "data.table" )
    require( "parallel" )
    require( "Gviz" )
    require( "rtracklayer" )
    require( "dplyr" )
    require( "readr" )
    require( "viridis" )
    require( "GenomicRanges" )
    require( "biomaRt" )
        
    ### Personal libraries
    bin_dir <- "./Bin_autoGviz/"
    source( paste0( bin_dir, "metafile_dealer.Rl" ) )
    source( paste0( bin_dir, "track_manager.Rl" ) )
    source( paste0( bin_dir, "mygviz.Rl" ) )

})

## FUNCTIONS
get_from_header <- function( string, metafile ){

    found <- as.character(

        unlist(

            fread(

                cmd = paste0(
                    "grep \"^#",
                    string,
                    "\" ",
                    metafile ),

                header = F,
                sep = '\t'

            )

        )

    )

    
    if( length( found ) > 0 ){

        return( found[ -1 ] )
        
    } else{

        return( "NOT FOUND" )

    }

}

isColor <- function( x ){
    
    res <- try( col2rgb( x ), silent = TRUE )
    return( ! "try-error" %in% class( res ) )

}

## Works only with string and hex, not with rgb!
color2fill <- function( x, alpha = 0.1 ){
    
    color <- do.call( rgb, as.list( col2rgb( x ) / 255 ) )

    alpha_hex <- as.character( as.hexmode( round( alpha * 255 ) ) )
    
    fill <- paste0( color, alpha_hex )

    return( fill )

}

## INPUT
metafile <- args[ 1 ] # Metafile where all the information for plotting the figure are stored (i.e. "fig1b.tsv" )
plot_name <- args[ 2 ] # Name of the plot (i.e. "fig1b")
cores <- args[ 3 ] # Number of cores used by the analysis

if( length( args ) == 4 ){ aux_folder <- args [ 4 ] } else{ aux_folder <- "./AUX_Gviz/" }
    
## PROCESSING
### Check it the right metafile have been passed to the script
check_file( metafile )

### Loading explicit information from header
#### Position
pos <- get_pos( metafile ) 

#### Genome used for the analysis
genome <- get_genome( metafile )

#### Text size
cex <- get_cex( metafile )

#### Background color
axis_bg <- get_axis_background( metafile )

#### Width of the plot
plot_width <- get_width( metafile )

#### Height of the plot
plot_height <- get_height( metafile )

### Get the info about the ideotrack
ideo_sizes <- get_ideosizes( metafile )

### Import the metafile 
meta_df <- fread( cmd = paste( "grep -v '^#'", metafile ), data.table = FALSE, header = TRUE )
print( "Metafile loaded!" )

### Get the track sizes
sizes <- get_sizes( meta_df$Extra, ideo_sizes )

### If there are tracks to overlay...
to_overlay <- extract_to_overlay( meta_df )

#### Look for groups and arrange names
overlay_and_labels <- group_dealer( to_overlay )

group_lists <- overlay_and_labels[[ 1 ]]
to_overlay <- overlay_and_labels[[ 2 ]]

if( is.data.frame( to_overlay ) ){
    
#### ... make a new list of names merging the overlay ones...
    merged_names <- merge_overlay_names( meta_df )
    
#### ... update the list of sizes (some of them have to be removed)...
    sizes <- update_overlay_sizes( meta_df, sizes )

#### ... remove overlayable tracks from the normal analysis...
    meta_df <- remove_overlay( meta_df )
    
}

## Get the colors
meta_df <- set_colors( meta_df, color_path = "../colors_db.txt" )

### Make the plot object
plotlist <- plotlist_maker( meta_df, pos, genome, cores, aux_folder, cex, axis_bg )

#### If there are the tracks to overlay...
if( exists( "merged_names" ) ){
    
    #### ... overlay track groups...
    overlayed_tracks <- overlay_all(
        to_overlay = to_overlay,
        group_lists = group_lists,
        pos = pos,
        genome = genome,
        aux_folder = aux_folder,
        cex = cex,
        axis_bg = axis_bg,
        cores = cores )

    #### ... and then join them in "plotlist" (in the correct order)
    plotlist <- c( plotlist, overlayed_tracks )[ merged_names ]
    
}

## OUTPUT
### Make the plot
pdf(
    paste0( plot_name, ".pdf" ),
    height = plot_height,
    width = plot_width
)

plotTracks(
        
    trackList = unlist( plotlist ),
    from = pos$from,
    to = pos$to,
    sizes = sizes
    add = FALSE
    
)

dev.off()

print( "Plot created, the process reaches the end, congratulations!" )

