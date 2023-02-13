## LIBRARIES
require( "rtracklayer" )
require( "Gviz" )
require( "biomaRt" )

## FUNCTIONS
### Make the compartments tracks
multicol_anno_maker <- function(
                                bedfile, name, color_column, coord_gr,
                                background.title = "transparent",
                                fontcolor.title = "black"
                                ){

    comp_gr <- import.bed( bedfile, which = coord_gr )
    colors <- as.data.frame( comp_gr )[ , color_column ] 

    comp_anno <- AnnotationTrack( range = comp_gr, name = name )

    displayPars( comp_anno )$fill <- colors
    displayPars( comp_anno )$background.title <- background.title
    displayPars( comp_anno )$fontcolor.title <- fontcolor.title
    
    return( comp_anno )
    
}

### Make the polygon tracks
polygon_track_maker <- function( trackfile, coord_gr, name, range, colors, window = 500 ){

    bw <- import.bw( trackfile, which = coord_gr )

    dt <- DataTrack(
    
        bw,
        name = name,
        type = "polygon",
        fill.mountain = colors,
        ylim = range,
    
        window = window,

        col = "transparent",
        col.axis = "black",
    
        background.title = "transparent",
        fontcolor.title = "black"
    
    )

    return( dt )
    
}

### Overlay tracks [not tested]
overlay_track_maker <- function( bwfiles, sample_id, names, coord_gr, colors, window = 1000 ){

    tracks <- lapply(
        names,
        function( name ){

            bw <- import.bw( bwfiles[ name ], which = coord_gr )

            color <- ifelse(
                colors[ name ] != names[ 1 ],
                color[ name ],
                colors
            ) ## The first track plotted should have all the colors to be show them in the legend
            
            dt <- DataTrack(
                bw,
                name = sample_id,
                
                groups = factor( name, levels = names ),
                col = color,
                window = window,
                
                type = 'a',
                legend = TRUE,
                col.axis = "black",
                background.title = "transparent",
                fontcolor.title = "black"
            )

            return( dt )

        }
    )

    overlay <- OverlayTrack( rev( tracks ) )

    return( overlay )
    
}


## INPUT
### General data
coords <- list(
    
    HOXAext = c( chr = "chr7",
           start = 26950000,
           end = 27350000 ), # Aext
    
    HOXDext = c( chr = "chr2",
           start = 175810000,
           end = 176460000 ) # Dext
)

genome <- "hg38"
biomart_dataset <- "hsapiens_gene_ensembl"

### Annotrack data
compartments <- c(
    C002 = "/path/to/C002_compartments.bed",
    C004_r1 = "/path/to/C004-r1_compartments.bed",
    C004_r2 = "/path/to/C004-r2_compartments.bed",
    HiC = "/path/to/Hi-C_compartments.bed"
)

### Track data
lits <- c(
    H3K27me3 = "path/to/H3K27me3_comparison.bw",
    H3K4me3 = "path/to/H3K4me3_comparison.bw" 
)
ranges <- list(
    H3K27me3 = c( 0, 6 ),
    H3K4me3 = c( 0, 6 ),
    RNAseq_C002 = c( 0, 15 ),
    RNAseq_C004 = c( 0, 15 )
)
colors <- list(
    H3K27me3 = c( "#ABD9E9", "#ABD9E9" ),
    H3K4me3 = c( "#BC1728", "#BC1728" ),
    RNAseq_C002 = c( "#C01B28", "#C01B28" ),
    RNAseq_C004 = c( "#C01B28", "#C01B28" )
)
windows <- c(
    H3K27me3 = 200,
    H3K4me3 = 200,
    RNAseq_C002 = 1000,
    RNAseq_C004 = 1000
)

### Output
out_dir <- "./Output/"

## PROCESSING
for( coord_name in names( coords ) ){

    plot <- paste0( out_dir, coord_name, ".pdf" )
    coord <- coords[[ coord_name ]]
    
    print( coord_name )
    print( coord )
    
### Trasform coordinates in a GRanges object
    coord_gr <- GRanges(
        seqnames = coord[ "chr" ],
        IRanges(
            start = as.integer( coord[ "start" ] ),
            end = as.integer( coord[ "end" ] )
        )
    )

### Ideogram + ruler
    ideog <- IdeogramTrack( genome = genome, chromosome = coord[ "chr" ] )
    ruler <- GenomeAxisTrack()

### Compartments
    compartment_annos <- lapply(
        names( compartments ),
        function( name ){

            multicol_anno_maker(
                bedfile = compartments[ name ],
                name = name,
                coord_gr = coord_gr,
                color_column = "itemRgb"
            )
            
        }
    )

### Literature
    lit_tracks <- lapply(
        names( lits ),
        function( name ){

            polygon_track_maker(
                trackfile = lits[ name ],
                coord_gr = coord_gr,
                name = name,
                range = ranges[[ name ]],
                colors = colors[[ name ]],
                window = windows[ name ]
            )
            
        }
    )
    
### Tracks overlay
    overlay_tracks <- lapply(
        names( overlay_groups ),
        function( group_name ){
            
            overlay_track_maker(
                bwfiles = overlay_groups[[ group_name ]],
                sample_id = group_name,
                names = names( overlay_groups[[ group_name ]] ),
                coord_gr = coord_gr,
                colors = overlay_group_colors[[ group_name ]],
                window = overlay_group_windows[[ group_name ]]
            )

        }
    )

### Gene reference
    biomart <- useMart( "ensembl", dataset = biomart_dataset )

    gene_track <- BiomartGeneRegionTrack(
        biomart = biomart,
        genome = genome,
        chromosome = coord[ "chr" ],
        start = as.integer( coord[ "start" ] ),
        end = as.integer( coord[ "end" ] ),
        stacking = "squish",
        name = "Genes",
        collapseTranscripts = TRUE,
        showId = TRUE,
        shape = "arrow",
        transcriptAnnotation = "symbol"
    )


    ## OUTPUT
    pdf( plot, width = 12, height = 8 )

    plotTracks(

        chr = coord[ "chr" ],
        from = as.integer( coord[ "start" ] ),
        to = as.integer( coord[ "end" ] ),
        
        trackList = c(
            ideog,
            ruler,
            compartment_annos,
            lit_tracks,
            gene_track
        ),

        main = coord_name,
        sizes = c( 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3 )
    )

    dev.off()

}
