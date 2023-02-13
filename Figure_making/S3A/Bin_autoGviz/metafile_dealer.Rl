## It checks if the file is an autoGviz valid metafile
check_file <- function( file, standard_header = "#gViz_metafile", standard_fields = "Track_name\tTrack_path\tTrack_upcolor\tTrack_downcolor\tTrack_type\tYlim\tWindow\tPeaks_bool\tPeaks_path\tPeaks_col\tExtra" ){

    header <- as.character( unlist( fread( cmd = paste( "head -1 ", file ), header = F ) ) )

    if( header != standard_header ){

        print( "This is not a gViz_metafile!!!" )
        quit()

    }


    fields <- paste( as.character( unlist( fread( cmd = paste( "grep ^Track_name ", file ), header = F ) ) ), collapse = "\t" )

    if( fields != standard_fields ){

        print( "Wrong fields numbers or names!!!" )
        quit()

    }

    print( "File check: Passed!" )

}

## It trasforms the genomic coordinates in a script readable form (it take in input UCSC format coordinate string)
split_coord_info_UCSC <- function( coord ){

    chr <- strsplit( coord, ':' )[[ 1 ]][ 1 ]
    range <- strsplit( coord, ':' )[[ 1 ]][ 2 ]

    from_string <- strsplit( range, '-' )[[ 1 ]][ 1 ]
    to_string <-  strsplit( range, '-' )[[ 1 ]][ 2 ]

    from <- as.numeric( gsub( ',', '', x = from_string ) )
    to <- as.numeric( gsub( ',', '', x = to_string ) )

    return( list( chr = chr, from = from, to = to ) )

}


## It extract the position information from the metafile
get_pos <- function( file ){

    pos_line <- as.character( unlist( fread( cmd = paste( "grep '^#pos' ", file ), header = F, sep = '\t' ) ) )

    if( length( pos_line ) == 0 ){

        print( "No position detected!!!" )
        quit()

    }

    pos <- pos_line[ 2 ]

    print( "Coordinates loaded!" )
    
    return( split_coord_info_UCSC( coord = pos ) )

}

## It extract the genome information from the metafile
get_genome <- function( file ){

    genome_line <- as.character( unlist( fread( cmd = paste( "grep '^#genome' ", file ), header = F ) ) )

    if( length( genome_line ) == 0 ){

        print( "No genome detected!!!" )
        quit()

    }

    print( "Genome info got!" )

    return( genome_line[ 2 ] )

}

## It extracts from the metafile the cex dimension to use
get_cex <- function( file ){

    cex_line <- as.character( unlist( fread( cmd = paste( "grep '^#cex' ", file ), header = FALSE ) ) )

    if( length( cex_line ) == 0 ){

        print( "No cex detected!!!" )
        return( 1 )

    }

    print( "Cex got!" )
    
    return( as.numeric( cex_line[ 2 ] ) )

}

## It extracts from the metafile the axis baground color to use
get_axis_background <- function( file ){

    axis_background_line <- as.character( unlist( fread( cmd = paste( "grep '^#ax_bg' ", file ), header = F ) ) )

    if( length( axis_background_line ) == 0 ){

        print( "No axis background detected!!!" )
        return( "white" )
        
    }

    print( "Axis background got!" )
    
    return( axis_background_line[ 2 ] )
    
}

## They extract the info about the output plot dimension
get_height <- function( file ){

    height_line <- as.character( unlist( fread( cmd = paste( "grep '^#height' ", file ), header = FALSE ) ) )

    if( length( height_line ) == 0 ){

        print( "No height detected!!!" )
        return( 8 )

    }

    print( "Height got!" )
    return( as.integer( height_line[ 2 ] ) )

}

get_width <- function( file ){

    width_line <- as.character( unlist( fread( cmd = paste( "grep '^#width' ", file ), header = F ) ) )

    if( length( width_line ) == 0 ){

        print( "No width detected!!!" )
        return( 10 )

    }

    print( "Width got!" )

    return( as.integer( width_line[ 2 ] ) )

}

## Get the sizes...
## ... for the ideotracks
get_ideosizes <- function( file ){

    ideosizes_line <- as.character( unlist( fread( cmd = paste( "grep '^#ideosizes' ", file ), header = F, sep = '\t' ) ) )

    if( length( ideosizes_line ) == 0 ){

        return( c( 1, 1 ) )

    }

    ideosizes_num_string <- ideosizes_line[ 2 ]
    ideosizes <- as.numeric( unlist( strsplit( ideosizes_num_string, ',' ) ) )

    if( length( ideosizes ) < 2 ){

        return( c( 1, 1 ) )

    } else{

        print( "Ideosizes got!" )
        return( ideosizes )

    }

}

### ... for all the tracks
get_sizes <- function( field, ideo_sizes = c( 1, 1 ) ){

    field_string <- paste0( field, collapse = '' )
    
    if( grepl( ";sizes:", field_string ) ){ 
        
        field_length <- length( field ) + 2

        ## Select only "sizes" information, discarding everything else like "elm_type"
        separates <- unlist( strsplit( field_string, ';' ) )
        only_sizes_index <- grep( "sizes:", separates )
        sizes_rows <- separates[ only_sizes_index ]

        ## Collapse everything in a string
        field_string <- paste0( ';', sizes_rows, collapse = ';' )

        ## Extract the info of size
        new_field_string <- gsub( ";sizes:", ',', field_string )
        new_field_string <- gsub( ';', '', new_field_string )
        new_field_string <- gsub( "^,", '', new_field_string )
        new_field <- as.numeric( unlist( strsplit( new_field_string, ',' ) ) )
        
        new_field_sizes <- c( ideo_sizes, new_field )

        ## If some info is missing it replace it with 1
        if( length( new_field_sizes ) < field_length ){

            missing <- field_length - length( new_field_sizes )
            missing_sizes <- rep( 1, missing )
            
            new_field_sizes <- c( new_field_sizes, missing_sizes )
        }

        ## Usually all the tracks have the same size
    } else{

        field_length <- length( field )
        field_sizes <- rep( 1, field_length )

        new_field_sizes <- c( ideo_sizes, field_sizes )

    }

    return( new_field_sizes )

}

## Get the colors
get_colors <- function( field, color_path = "/path/to/colors_db.txt" ){

    ## Get official colors
    if( !file.exists( color_path ) ){

        print( "The color database does not exist!" )
        quit()

    }

    color_df <- fread( file = color_path, header = TRUE, data.table = FALSE )  

    ## Get track info from file
    field_string <- paste0( field, collapse = '' )

    change_list <- list()
    
    if( grepl( ";elm_type:", field_string ) ){
        
        line_n <- 0
        for( line in field ){

            line_n <- line_n + 1

            ## If in that line there is an info about what the track represent...
            if( grepl( ";elm_type:", line ) ){

                ## ... select that info
                line_els <- unlist( strsplit( line, ';' ) )

                track_info_index <- grep( "elm_type:", line_els )
                track_info <- line_els[ track_info_index ]

                tracks <- gsub( "elm_type:", '', track_info )
                uptrack <- gsub( ",.*", '', tracks )
                dntrack <- gsub( ".*,", '', tracks )

                ## Convert that info in the color of the track
                uptrack_color <- color_df[ which( color_df$Mark == uptrack ), "Color1" ]
                dntrack_color <- color_df[ which( color_df$Mark == dntrack ), "Color1" ]

                element_n <- paste0( "line_n:", line_n )
                change_list[[ element_n ]] <- c( line_n, uptrack_color, dntrack_color )

            }
                
        }

    }

    return( change_list )

}

## Set the colors considering even if there are official ones to use that have been set up in a database
set_colors <- function( meta_df ){

    ### If the tracks have official colors override the normal one
    replace_colors <- get_colors( meta_df$Extra )
    if( length( replace_colors ) > 0 ){
    
        for( replace in replace_colors ){

            replace <- unlist( replace )
            replace_n <- as.numeric( replace[ 1 ] )
            
            meta_df$Track_upcolor[ replace_n ] <- replace[ 2 ]
            meta_df$Track_downcolor[ replace_n ] <- replace[ 3 ]
        
        }

    }

    return( meta_df )
    
}

## Check if in the metafile there are information to make the plot logartmic or apply the mean
log_scale_check <- function( field ){

    return( grepl( ";log_scale=TRUE;", field ) )

}

less_mean_check <- function( field ){

    return( grepl( ";less_mean=TRUE;", field ) )

}

## Determine which are the groups to overlay
get_overlay_groups <- function( col ){

    all_groups <- gsub( ".*;O:", '', col )
    all_groups <- gsub( ";.*", '', all_groups )
 
    return( all_groups )

}

## Get which are the tracks to overlay (if there are) and return it
extract_to_overlay <- function( meta_df ){

    overlay_rows <- grep( ";O:.*;", meta_df$Extra )

    if( length( overlay_rows ) == 0 ){

        return( "No" )

    }

    overlay_df <- meta_df[ overlay_rows, ]

    return( overlay_df )
    
}

## If there are tracks to overalay in the meta_df it create a new list of Track_name take in consideration that multiple tracks should be merged in a single one (so we have just one name for them)
merge_overlay_names <- function( meta_df, uniq = TRUE ){

    original_names <- c( "cytoband", "genome_axis", meta_df$Track_name )
    group_names <- c( '', '', get_overlay_groups( meta_df$Extra ) )

    new_names <- original_names 
    new_names[ which( group_names != '' ) ] <- group_names[ which( group_names != '' ) ]

    new_names <- gsub( "\\(.*\\)", '', new_names ) ## To manage group extra labels
    
    if( uniq ){

        new_names <- unique( new_names )

    }

    return( new_names )

}

## It calculates the new vector of sizes (in case of overlay)
update_overlay_sizes <- function( meta_df, sizes ){

    original_names <- c( "cytoband", "genome_axis", meta_df$Track_name )
    new_names <- merge_overlay_names( meta_df, FALSE )
    
    names_db <- data.frame( new_names, original_names )
    names_db$group_bool <- sapply( 1:dim( names_db )[ 1 ], function( x ){ ifelse( names_db[ x, 1 ] == names_db[ x, 2 ], 'N', 'Y' ) })

    sizes_names <- c( "cytoband", "genome_axis", sapply( meta_df$Track_name, function( x ){ rep( x, 3 ) }) )

    sizes_df <- data.frame( names = sizes_names, sizes = sizes )

    sizes_df[ , c( "names", "group_bool" ) ] <- names_db[ sapply( sizes_df$names, function( x ){ grep( paste0( '^', x, '$' ), names_db$original_names ) }), c( "new_names", "group_bool" ) ]

    sizes_df$peak_bool <- c( 'N', 'N', sapply( meta_df$Track_name, function( x ){ c( 'N', 'Y', 'N' ) }) )
    
    sizes_df <- sizes_df[ which( sizes_df$group_bool != 'Y' | sizes_df$peak_bool != 'Y' ), ] 

    sizes <- unlist( unique( sizes_df )[ 2 ] )

    return( sizes )
    
}

## It removes the tracks to overlay from a metafile
remove_overlay <- function( meta_df ){

    overlay_rows <- grepl( ";O:.*;", meta_df$Extra )
    meta_df <- meta_df[ !overlay_rows, ]

    return( meta_df )
    
}

## Look for tracks to group and in case store the names in a list of vector
group_dealer <- function( to_overlay ){

    if( to_overlay == "No" ){ return( NULL ) }
    
    ## Look for groups
    find <- sum( grepl( "^.*;O:.*\\(.*\\);.*;$", to_overlay$Extra ) )

    if( find > 0 ){

        ## Try to deal this part with a dictionary
        ## If group/s has/ve been found...
        ## ... take the group name
        group_labels <-  gsub( ".*\\(", '', to_overlay$Extra )
        group_labels <- gsub( "\\).*", '', group_labels )

        ## ... and the label names
        overlay_labels <- gsub( "\\(.*", '', to_overlay$Extra )
        overlay_labels <- gsub( ";O:", '', overlay_labels )

        ## Create a dictionary to associate overlay tracks and group names (i.e. {tooverlayX:togroupD,togroupB,togroupG})
        overlay_db <- data.frame( overlay = overlay_labels, group = group_labels )

        group_lists <- list()
        
        for( overlay_label in unique( overlay_labels ) ){

            temp_list <- overlay_db[ which( overlay_db$overlay == overlay_label ), "group" ]
            group_lists[[ overlay_label ]] <- temp_list

        }

        ## Rename the to_overlay$Extra field (removing the group names)
        to_overlay$Extra <- gsub( "\\(.*\\)", ';', to_overlay$Extra )

        return( list( group_lists = group_lists, to_overlay = to_overlay ) )
        
    } else{

        ## If not group has been found return NA for group and to_overlay (without doing any changes)
        return( list( group_lists = NA, to_overlay = to_overlay ) )

    }
}

check_horiz <- function( metaline ){

    return( grepl( ";horiz;", metaline$Extra ) )
    
}
