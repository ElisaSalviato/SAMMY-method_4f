## LIBRARIES
require( "parallel" )
require( "data.table" )
require( "rtracklayer" )
require( "preprocessCore" )

## FUNCTIONS
bins_calculator <- function( genome_bedfile, binsize, genome, dropping_chrs = c( "chrM", "M", "MT", "chrMT") ){

    ## Import genome
    hg38_df <- fread( genome_bedfile, data.table = FALSE )
    names( hg38_df ) <- c( "seqnames", "start", "end" )

    hg38_gr <- makeGRangesFromDataFrame( df = hg38_df, start.field = "start", end.field = "end", seqnames.field = "seqnames" )

    ## Add information to the genome gr object
    seqlevels( hg38_gr ) <- hg38_df$seqnames
    seqlengths( hg38_gr ) <- hg38_df$end
    genome( hg38_gr ) <- genome

    ## Remove unuseful chromosomes
    hg38_gr <- keepStandardChromosomes( hg38_gr, pruning.mode = "coarse" )
    hg38_gr <- dropSeqlevels( hg38_gr, dropping_chrs, pruning.mode = "coarse" )

    ## Store in a variable the chromosome lengths to use for calculate the bins
    chr_seqlengths <- seqlengths( hg38_gr )

    ## Create the bins to use for the following correlation analysis
    bin_list <- tileGenome( chr_seqlengths, tilewidth = binsize, cut.last.tile.in.chrom = TRUE )

    return( bin_list )
    
}

import_and_rebin__bw <- function( files, bin_list, genome, names, cores ){

    ## cores <- length( files )

    bws <- mclapply( files, mc.cores = cores, function( file ){
        
        ## Import the bigwig as RleLists
        bwR <- import.bw( file, as = "RleList" )

        ## Sort the "bin_list" seqlevels names to make them coincide with the bigwig imported
        bin_list__names <- seqlevels( bin_list )
        bwR__grepnames <- paste0( '^', names( bwR ), '$' )
    
        chr_bwR_order <- sapply( bwR__grepnames, function( chr ){

            grep( chr, bin_list__names )

        })
        chr_bwR_order <- unlist( chr_bwR_order[ !is.na( as.integer( chr_bwR_order ) ) ] ) ## remove extra chromosomes from the list of analysis chromosomes
        
        seqlevels( bin_list ) <- bin_list__names[ chr_bwR_order ]
        bwR <- bwR[ seqlevels( bin_list ) ] ## remove extra chromosomes from Rle object

        ## Rebin the imported bigwig according to the previous calculated bins
        bw <- binnedAverage( bins = bin_list, numvar = bwR, varname = "score" )

        ## Add extra information to the rebinned bigwig
        genome( bw ) <- genome
    
        return( bw )

    })
    names( bws ) <- names

    return( bws )

}

correlator_by_chr <- function( bws1, names1, bws2, names2, chrs, correlation_type, cores ){

    ## Cycle1 on all the bigwig imported
    corr_list <- mclapply( names1, mc.cores = cores, function( name1 ){

        bw1 <- bws1[[ name1 ]]

        ## Cycle2 on all the bigwig imported
        name1_corr_list <- lapply( names2, function( name2 ){

        bw2 <- bws2[[ name2 ]]

            ## Cycle on each chromosome
            corr_by_chr__list <- lapply( chrs, function( chr ){

                bw1_chr_score <- bw1[ seqnames( bw1 ) == chr ]$score
                bw2_chr_score <- bw2[ seqnames( bw2 ) == chr ]$score

                ## Calculate the correlation among the two bigwig for the chr selected
                corr_val <- cor( x = bw1_chr_score, y = bw2_chr_score, method = correlation_type )

                ## Store the results in a vector
                df_row <- c( row = name1, col = name2, chr = chr, corr_val = corr_val ) 

            return( df_row )

            })

            ## Join the chromosome results in a dataframe
            corr_by_chr__df <- do.call( rbind, corr_by_chr__list )

        })
        ## Join the cycle2 results in dataframe
        name1_corr_df <- do.call( rbind, name1_corr_list )

    return( name1_corr_df )
    
    })

    ## Join the list of correlation in a dataframe
    corr_df <- as.data.frame( do.call( rbind, corr_list ) )

    corr_df$row <- as.character( corr_df$row )
    corr_df$col <- as.character( corr_df$col )
    corr_df$chr <- as.character( corr_df$chr )
    corr_df$corr_val <- as.numeric( as.character( corr_df$corr_val ) )

    return( corr_df )

}

perchr_to_genome_corr <- function( chr_df ){

    df <- data.frame()
    
    for( col in unique( chr_df$col ) ){

        for( row in unique( chr_df$row ) ){

            corr_mean <- mean( chr_df[ which( chr_df$col == col & chr_df$row == row ), "corr_val" ] )

            row_df <- data.frame( cbind( row = row, col = col, corr_val = corr_mean ) )

            df <- data.frame( rbind( df, row_df ) )
        }
        
    }

    df$corr_val <- as.numeric( as.character( df$corr_val ) )

    return( df )
    
}

ztransform <- function( data ){

    return( data - mean( data ) ) / sd( data )

}

qnorm_BWlist <- function( bws, cores = 4 ){

    bw <- leftjoin_granges( bws )
    bw_df <- as.data.frame( bw )
        
    bw_score_mat <- as.matrix( bw_df[ , seq( 6, dim( bw_df )[ 2 ] ) ] )
    colnames( bw_score_mat ) <- names( bw_df[ , seq( 6, dim( bw_df )[ 2 ] ) ] )
        
    bw_normq <- normalize.quantiles( bw_score_mat )
    names( bw_normq ) <- names( bw_score_mat )
        
    bw_normq_df <- cbind( bw_df[ , seq( 5 ) ], as.data.frame( bw_normq ) )
    names( bw_normq_df ) <- c( names( bw_df[ , seq( 5 ) ] ), colnames( bw_score_mat ) )

    normq_bws <- mclapply( names( bws ), mc.cores = cores, function( bw_name ){

        score( bws[[ bw_name ]] ) <- bw_normq_df[ , bw_name ]
        return( bws[[ bw_name ]] )
        
    })
    names( normq_bws ) <- names( bws )

    return( normq_bws )
    
}

fromBWlist_to_ggdf <- function( bws ){

    bws_df <- lapply( names( bws ), function( bw_name ){

        bw_df <- as.data.frame( bws[[ bw_name ]] )
        bw_df$name <- bw_name
        bw_df$zscore <- ztransform( bw_df$score )
        
        return( bw_df )
        
    })

    df <- as.data.frame( do.call( "rbind", bws_df ) )

    return( df )
    
}


leftjoin_granges <- function( bws ){

    df <- as.data.frame( bws[[ 1 ]] )[ , c( 1, 2, 3, 4, 5 ) ]

    for( elm in names( bws ) ){

        vals <- as.data.frame( bws[[ elm ]] )[ , 6 ]

        df <- cbind( df, vals )
        names( df ) <- c( names( df )[ -length( names( df ) ) ], elm )

    }

    joint_bws <- makeGRangesFromDataFrame( df, keep.extra.columns = TRUE )

    return( joint_bws )
    
}

calculate_variance <- function( bws, peaks, cores = 4 ){

    bw <- leftjoin_granges( bws )

    sd_vect <- apply( as.data.frame( bw )[ , seq( 6, 5 + length( names( bws ) ) ) ], 1, sd )
    var_vect <- apply( as.data.frame( bw )[ , seq( 6, 5 + length( names( bws ) ) ) ], 1, var )
    mean_vect <- apply( as.data.frame( bw )[ , seq( 6, 5 + length( names( bws ) ) ) ], 1, mean )

    bw$sd <- sd_vect
    bw$variance <- var_vect
    bw$mean <- mean_vect
    
    domvar_list <- mclapply( names( peaks ), function( peak_name ){

        peak <- peaks[[ peak_name ]]
        bins_in_peaks__indices <- findOverlaps( peak, bw )@to

        bins_in_peaks__gr <- bw[ bins_in_peaks__indices, ]
        bins_in_peaks__df <- as.data.frame( bins_in_peaks__gr )

        hist_df <- cbind( bins_in_peaks__df, histone = peak_name )

        return( hist_df )

    })
    domvar_df <- as.data.frame( do.call( rbind, domvar_list ) )
    
    return( domvar_df )
    
}

check_variability <- function( comp_db, peak_db, bin_list, binsize, varplot_title, genome = "hg38", cores = 8 ){

    ### Import track data
    S2LvsS3_bws <- import_and_rebin__bw(
        files = comp_db$File,
        bin_list = bin_list,
        genome = genome,
        names = comp_db$Patient_name,
        cores = cores
    )

    comp_ggdf <- fromBWlist_to_ggdf( S2LvsS3_bws )

    qnorm_S2LvsS3_bws <- qnorm_BWlist( S2LvsS3_bws ) 
    comp_qnorm_ggdf <- fromBWlist_to_ggdf( qnorm_S2LvsS3_bws )

    ### Import domains data
    peak_beds <- mclapply( peak_db$Patient_name, mc.cores = cores, function( name ){

        path <- peak_S2LvsS3_db$File[ which( peak_db$Patient_name == name ) ]

        return( import.bed( path ) )

    })
    names( peak_beds ) <- peak_db$Patient_name

    ### Calculate variance per bin across patients per histone
    var_df <- calculate_variance( qnorm_S2LvsS3_bws, peak_beds, cores )
    var_df$histone <- gsub( "vs.*", '', var_df$histone )
    var_df$domain <- paste0( var_df$seqnames, ':', var_df$start, '-', var_df$end )

    ## Make plots
    dist <- ggplot( comp_ggdf, aes( x = patient, y = score ) ) +
        geom_boxplot() +
        ggtitle( label = "mle score" ) +
        theme_classic()
    distz <- ggplot( comp_ggdf, aes( x = patient, y = zscore ) ) +
        geom_boxplot() +
        ggtitle( label = "ztransformed( mle score )" ) +
        theme_classic()
    distq <- ggplot( comp_qnorm_ggdf, aes( x = patient, y = score ) ) +
        geom_boxplot() +
        ggtitle( label = "quantile normalizated" ) +
        theme_classic()

    varplot <- ggplot( var_df, aes( x = histone, y = variance ) ) +
        geom_boxplot() +
        geom_text_repel( aes( label = ifelse(  variance > 0.019, as.character( domain ), '' ) ) ) + 
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "After quantile normalization" ) +
        theme_classic()

    return( list(

        dist = dist,
        distz = distz,
        distq = distq,
        varplot = varplot,

        dist_df = comp_ggdf,
        distq_ggdf = comp_qnorm_ggdf,
        var_df = var_df,

        domains_bed = peak_beds,
        comp_bws = S2LvsS3_bws,
        comp_qnorm_bws = qnorm_S2LvsS3_bws

    ) )
    
}

tracks_import__and__distplot <- function( comp_db, bin_list, genome, cores = 8 ){

    S2LvsS3_bws <- import_and_rebin__bw(
        files = comp_db$File,
        bin_list = bin_list,
        genome = genome,
        names = comp_db$Patient_name,
        cores = cores
    )

    comp_ggdf <- fromBWlist_to_ggdf( S2LvsS3_bws )

    qnorm_S2LvsS3_bws <- qnorm_BWlist( S2LvsS3_bws ) 
    comp_qnorm_ggdf <- fromBWlist_to_ggdf( qnorm_S2LvsS3_bws )

    dist <- ggplot( comp_ggdf, aes( x = patient, y = score ) ) +
        geom_boxplot() +
        ggtitle( label = "mle score" ) +
        theme_classic()
    distz <- ggplot( comp_ggdf, aes( x = patient, y = zscore ) ) +
        geom_boxplot() +
        ggtitle( label = "ztransformed( mle score )" ) +
        theme_classic()
    distq <- ggplot( comp_qnorm_ggdf, aes( x = patient, y = score ) ) +
        geom_boxplot() +
        ggtitle( label = "quantile normalizated" ) +
        theme_classic()

    
    return(

        list(
            comp_ggdf = comp_ggdf,
            qnorm_bws = qnorm_S2LvsS3_bws,
            comp_qnorm_ggdf = comp_qnorm_ggdf,
            dist_plot = dist,
            distz_plot = distz,
            distq_plot = distq
        )
        
    )
    
}

import_peaks <- function( peak_db, cores = 8 ){

    peak_beds <- mclapply( peak_db$Patient_name, mc.cores = cores, function( name ){

        path <- peak_S2LvsS3_db$File[ which( peak_db$Patient_name == name ) ]

        return( import.bed( path ) )

    })
    names( peak_beds ) <- peak_db$Patient_name

    return( peak_beds )
    
}

variance_per_domain <- function( qnorm_bws, peak_beds, cores = 8 ){

    var_df <- calculate_variance( qnorm_bws, peak_beds, cores )
    var_df$histone <- gsub( "vs.*", '', var_df$histone )
    var_df$domain <- paste0( var_df$seqnames, ':', var_df$start, '-', var_df$end )
    
    return( var_df )

}

plot_domvariance_per_group <- function( var_df, outpath, varplot_title = '' ){

    colvect <- var_df$color
    names( colvect ) <- var_df$histone

    ## Plots using variance
    histx_varplot <- ggplot( var_df, aes( x = histone, y = variance, fill = histone ) ) +
        geom_boxplot() +
        ## geom_text_repel( aes( label = ifelse(  variance > 0.019, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ group, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing variance per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    groupx_varplot <- ggplot( var_df, aes( x = group, y = variance, fill = histone ) ) +
        geom_boxplot() +
        ## geom_text_repel( aes( label = ifelse(  variance > 0.019, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ histone, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing variance per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    ## Plots using sd
    histx_sdplot <- ggplot( var_df, aes( x = histone, y = sd, fill = histone ) ) +
        geom_boxplot() +
        geom_text_repel( aes( label = ifelse( sd > 0.2, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ group, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing standard deviation per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    groupx_sdplot <- ggplot( var_df, aes( x = group, y = sd, fill = histone ) ) +
        geom_boxplot() +
        geom_text_repel( aes( label = ifelse( sd > 0.2, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ histone, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing standard deviation per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    ## Plots using mean
    histx_meanplot <- ggplot( var_df, aes( x = histone, y = mean, fill = histone ) ) +
        geom_boxplot() +
        ## geom_text_repel( aes( label = ifelse(  variance > 0.019, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ group, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing mean per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    groupx_meanplot <- ggplot( var_df, aes( x = group, y = mean, fill = histone ) ) +
        geom_boxplot() +
        ## geom_text_repel( aes( label = ifelse(  variance > 0.019, as.character( domain ), '' ) ) ) +
        facet_grid( . ~ histone, ) +
        ggtitle( label = paste0( varplot_title, " binsize: ", binsize ), subtitle = "Comparing mean per group" ) +
        scale_fill_manual( values = colvect ) +
        theme_classic()

    ## Create a toy samples of the plots
    pdf( outpath )
    
    print( histx_varplot )
    print( groupx_varplot )

    print( histx_sdplot )
    print( groupx_sdplot )

    print( histx_meanplot )
    print( groupx_meanplot )

    dev.off()

    return( list(
        histx_varplot = histx_varplot,
        groupx_varplot = groupx_varplot,
        histx_sdplot = histx_sdplot,
        groupx_sdplot = groupx_sdplot,
        histx_meanplot = histx_meanplot,
        groupx_meanplot = groupx_meanplot
    ) )
    
}

### Transform the score of the bws in zscore values
zscor_BWlist <- function( bws, cores = 4 ){

    zbws <- mclapply( bws, mc.cores = cores, function( bw ){

        zscore <- ztransform( bw$score )
        bw$score <- zscore

        return( bw )

    })

    names( zbws ) <- names( bws )

    return( zbws )

}
