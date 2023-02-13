require( "Gviz" )
require( "rtracklayer" )
require( "dplyr" )
require( "readr" )
require( "viridis" )
require( "GenomicRanges" )

## Make the cytoband and genomeaxis
make_cytoband <- function( genome, chr, aux_folder, cex, fontcolor = "black" ){

    url_prefix <- "http://hgdownload.cse.ucsc.edu/goldenpath/"
    url_suffix <- "/database/cytoBand.txt.gz"
    url <- paste0( url_prefix, genome, url_suffix )

    cytoband_file <- paste0( aux_folder, "/cytoband.txt.gz" )
    
    download <- paste0( "wget -nv -O ", cytoband_file, ' ', url )
    less <- paste0( "less ", cytoband_file )

    system( command = download, ignore.stdout = TRUE, ignore.stderr = TRUE, intern = TRUE )
    cytoband_data <- fread( cmd = less, data.table = F )

    ideotrack <- IdeogramTrack( genome = genome,
                               chromosome = chr,
                               cex = cex,
                               fontcolor = fontcolor )

    return( ideotrack )
    
}

make_genomeaxis <- function( cex, color = "black", fontcolor = "black" ){

    axis_track <- GenomeAxisTrack( cex = cex,
                                  col = color,
                                  fontcolor = fontcolor )

    return( axis_track ) 
    
}

## My default track maker
track_maker <- function( file, format, name, genome, chromosome, type, ylim, line_color = "black", line_width = 0.01, fill_colors, window, cex.title = 2, axis_bg = "white", log_scale = FALSE, less_mean = FALSE ){
    
    range <- import( file, format = format )

    if( length( ylim ) != 2 ){

        ticks <- ylim
        ylim <- c( min( ylim ), max( ylim ) )

    }
    
    obj <- DataTrack( range = range,
                     name = name,
                     genome = genome,
                     chromosome = chromosome,
                     type = type,
                     ylim = ylim,
                     col = line_color,
                     lwd.mountain = line_width,
                     lwd = line_width,
                     fill.mountain = fill_colors,
                     window = window,
                     col.axis = "black", col.title = "black", background.title = axis_bg,
                     cex.title = cex.title, cex.axis = cex.title * 0.7 )

### Consider the case in which we want to specify the ticks
        if( exists( 'ticks' ) ){

            displayPars( obj ) <- list( yTicksAt = ticks )

        }

    
    if( log_scale == TRUE ){

        obj@data <- log( obj@data + 1 )

    }

    if( less_mean == TRUE ){

        mean <- mean( obj@data )
        obj@data <- obj@data - mean

    }

    return( obj )

}

## My default peak maker
peak_maker <- function( file, name, genome, chromosome, shape = "box", line_color, line_width = 1, fill_color, bg.title = "lightgrey", frame_col = "white", connector_col = "white" ){     

    obj <- AnnotationTrack( range = NULL,
                           name = "",
                           col = line_color,
                           lwd = line_width,
                           fill = fill_color,
                           background.title = bg.title,
                           col.frame = frame_col, col.line = connector_col )

    
    if( length( file ) > 0 ){

        if( file.exists( file ) ){
        
            obj <- AnnotationTrack( range = file,
                                   name = name,
                                   genome = genome,
                                   chromosome = chromosome,
                                   shape = shape,
                                   col = line_color,
                                   lwd = line_width,
                                   fill = fill_color,
                                   background.title = bg.title,
                                   col.frame = frame_col, col.line = connector_col )
        }

    }

    return( obj )

}

## Plot track and peak
track_and_peak <- function( track_path, sample_name, gen = "hg38", chr, ylim, track_colors, window, peak_path, peak_color, cex_title = 1.2, axis_bg, log_scale = FALSE, less_mean = FALSE, track_format, line_color = "black", type = "polygon" ){
    
    track <- track_maker( file = track_path,
                         format = track_format,
                         name = sample_name,
                         genome = gen,
                         chromosome = chr,
                         type = type,
                         ylim = ylim,
                         line_color = line_color,
                         line_width = 0.001,
                         fill_colors = track_colors,
                         window = window,
                         cex.title = cex_title,
			 axis_bg = axis_bg,
                         
                         log_scale = log_scale,
                         less_mean = less_mean )

    peak <- peak_maker( file = peak_path,
                       name = "",
                       genome = gen,
                       chromosome = chr,
                       shape = "box",
                       line_color = "black",
                       line_width = 1,
                       fill_color = peak_color,
                       bg.title = axis_bg,
                       frame_col = "white",
                       connector_col = "white" )

    separator <- peak_maker( file = peak_path,
                            name = "",
                            gen = gen,
                            chromosome = chr,
                            shape = "box",
                            line_color = "white",
                            line_width = 1,
                            fill_color = "white",
                            bg.title = "white",
                            connector_col = "white"  )

    
    return( list( track, peak, separator ) )

}

## Given a list of tracks it create the objects to plot with Gviz
plotlist_maker <- function( meta_df, pos, genome, cores, aux_folder, cex = 1, axis_bg = "grey", line_color = "black" ){

    ## If to handle pos that start from 0
    if( pos$from == 0 ){ pos$from <- 1 }

    ## Split the string of position in single elements
    pos_string <- paste0( pos$chr, '_', pos$from, '-', pos$to )
    pos_gr <- GRanges( pos$chr, IRanges( as.numeric( pos$from ), as.numeric( pos$to ) ) )

    chr <- pos$chr
    from <- pos$from
    to <- pos$to
    
    ## Make, if not exists, the folder where store the temp data
    system( command = paste( "mkdir -p", aux_folder ) )

    ## Calculate the number of elements to plot
    meta_lines__n <- dim( meta_df )[ 1 ]
    elements_vect <- seq( from = 1, to = meta_lines__n  )

    ## Make the plot objects
    plot_list <- mclapply( elements_vect, mc.cores = cores, function( element_n ){

        ## Set the format
        metaline <- meta_df[ element_n, ]

        ## Set the format
        track_format <- metaline$Track_type

        ## Check if is needed to calculate the log or the mean
        log_scale <- log_scale_check( metaline$Extra )
        less_mean <- less_mean_check( metaline$Extra )
        
        ### Sample name
        sample_name <- metaline$Track_name

        if( log_scale == TRUE ){ sample_name <- paste( sample_name, "(log)" ) }
        if( less_mean == TRUE ){ sample_name <- paste( sample_name, "(-mean)" ) }

        ### Track path
        track_path <- make_track_to_plot( metaline, pos_string, pos_gr, chr, to, aux_folder )
        
        ### Track colors
        color1 <- metaline$Track_downcolor
        color2 <- metaline$Track_upcolor
        track_colors <- c( color1, color2 )
        
        ### Peak path
        peak_path <- make_bed_to_plot( metaline, aux_folder )
        
        ### Peak color
        peak_color <- metaline$Peaks_col

        ### Window
        window <- metaline$Window

        ### Ylim
        ylim_string <- metaline$Ylim
        ylim <- as.numeric( unlist( strsplit( ylim_string, split = ',' ) ) )

        ### Calculate the plot
        plot <- track_and_peak(
            
            sample_name = sample_name,
            track_path = track_path,
            track_format = "bigWig",
            track_colors = track_colors,
            log_scale = log_scale,
            peak_path = peak_path,
            peak_color = peak_color,
            chr = chr,
            ylim = ylim,
            window = window,
            cex_title = cex,
            axis_bg = axis_bg,
            gen = genome,
            less_mean = less_mean
            
        )
        
    })
    names( plot_list ) <- meta_df$Track_name

    ### Make the cytoband
    cytoband <- make_cytoband( genome, chr, aux_folder, cex )
    
    ### Make genome axis track
    axis_track <- make_genomeaxis( cex )

    ### Merge cytoband, axis track and plots in a single list
    full_plot_list <- list()
    full_plot_list[[ "cytoband" ]]  <- cytoband
    full_plot_list[[ "genome_axis" ]] <- axis_track
    
    for( plot in names( plot_list ) ){
        
        full_plot_list[[ plot ]] <- plot_list[[ plot ]]

    }

    print( "Plot object created!" )

    return( full_plot_list )

}

overlay <- function( tracks_df, pos, genome, aux_folder, cex, axis_bg = "grey" ){

    ## If to handle pos that start from 0
    if( pos$from == 0 ){ pos$from <- 1 }

    ## Split the string of position in single elements
    pos_string <- paste0( pos$chr, '_', pos$from, '-', pos$to )
    pos_gr <- GRanges( pos$chr, IRanges( as.numeric( pos$from ), as.numeric( pos$to ) ) )

    chr <- pos$chr
    from <- pos$from
    to <- pos$to

    ## Make, if not exists, the folder where store the temp data
    system( command = paste( "mkdir -p", aux_folder ) )

    ## Make the plot objects
    tracks <- mclapply( tracks_df$Track_name, mc.cores = length( tracks_df$Track_name ), function( element_name ){

        metaline <- tracks_df[ which( tracks_df$Track_name == element_name ), ]

        ## Set the format
        track_format <- metaline$Track_type

        ## Check if is needed to calculate the log or the mean
        log_scale <- log_scale_check( metaline$Extra )
        less_mean <- less_mean_check( metaline$Extra )
        
        ### Sample name
        sample_name <- get_overlay_groups( metaline$Extra )

        if( log_scale == TRUE ){ sample_name <- paste( sample_name, "(log)" ) }
        if( less_mean == TRUE ){ sample_name <- paste( sample_name, "(-mean)" ) }

        ### Track path
        track_path <- make_track_to_plot( metaline, pos_string, pos_gr, chr, to, aux_folder )
        
        ### Track colors
        color <- metaline$Track_upcolor
        
        ### Window
        window <- metaline$Window

        ### Ylim
        ylim_string <- metaline$Ylim
        ylim <- as.numeric( unlist( strsplit( ylim_string, split = ',' ) ) )

        if( length( ylim ) != 2 ){

                ticks <- ylim
                ylim <- c( min( ylim ), max( ylim ) )

        }

        
        ### Calculate the plot
        range <- import( track_path, format = "BigWig" )

        ### Create the track
        single_track <- DataTrack(

            range = range,
            name = sample_name,
            genome = genome,
            chromosome = chr,
            type = 'a',
            window = window,
            col = color,
            col.axis = "black",
            col.title = "black",
            background.title = axis_bg,
            cex.title = cex,
            cex.axis = cex * 0.7,
            ylim = ylim,
            legend = TRUE,
            cex.legend = cex

        )

        ### Set extra parameters
        displayPars( single_track ) <- list(

            groups = factor( element_name,
                            levels = tracks_df$Track_name ),
            col = tracks_df$Track_upcolor, legend = TRUE

        )

        if( exists( 'ticks' ) ){

            displayPars( obj ) <- list( yTicksAt = ticks )

        }

        
        return( single_track )
        
    })

    ## Overlay the tracks and set extra parameter to the new track
    overlay_track <- OverlayTrack( trackList = tracks )
    displayPars( overlay_track ) <- list(

        cex.title = cex,
        cex.axis = cex * 0.7,
        cex = cex,
        background.title = axis_bg,
        cex.legend = cex

    )
    
    ## Make the separator track
    peak_path <- paste0( aux_folder, "/fake_peaks.bed" )
    system( command = paste0( "echo -e 'chrNULL\t0\t1\tnovalues' > ", peak_path ) )
    separator <- peak_maker( file = peak_path,
                            name = "",
                            gen = genome,
                            chromosome = chr,
                            shape = "box",
                            line_color = "white",
                            line_width = 1,
                            fill_color = "white",
                            bg.title = "white",
                            connector_col = "white"  )

    print( "Tracks overlayed" )
    
    return( list( overlay_track, separator ) )
    
}
    
## Divide in overlay tracks in groups of overlay and overlay each one separately
overlay_all <- function( to_overlay, pos, genome, aux_folder, cex, axis_bg ){

    all_groups <- get_overlay_groups( to_overlay$Extra )
    groups <- unique( all_groups )

    to_overlay <- set_colors( to_overlay )
    
    overlayed_tracks <- mclapply( groups, mc.cores = length( groups ), function( group ){

        tracks_df <- to_overlay[ which( all_groups == group ), ]
        track <- overlay( tracks_df, pos, genome, aux_folder, cex, axis_bg )

        return( track )
        
    })
    names( overlayed_tracks ) <- groups

    return( overlayed_tracks )
    
}
