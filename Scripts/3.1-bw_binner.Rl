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

import_and_rebin__bw <- function( files, bin_list, genome, names ){

    cores <- length( files )
    bws <- mclapply( files, mc.cores = cores, function( file ){
        
        ## Import the bigwig as RleLists
        bwR <- import.bw( file, as = "RleList" )
        bwR <- bwR[ seqnames( bin_list )@values ]
        
        ## Sort the "bin_list" seqlevels names to make them coincide with the bigwig imported
        bin_list__names <- seqlevels( bin_list )
        bwR__grepnames <- paste0( '^', names( bwR ), '$' )
    
        chr_bwR_order <- sapply( bwR__grepnames, function( chr ){

            grep( chr, bin_list__names )

        })
    
        seqlevels( bin_list ) <- bin_list__names[ chr_bwR_order ]

        ## Rebin the imported bigwig according to the previous calculated bins
        bw <- binnedAverage( bins = bin_list, numvar = bwR, varname = "score" )

        ## Add extra information to the rebinned bigwig
        genome( bw ) <- genome
    
        return( bw )

    })
    names( bws ) <- names

    return( bws )

}
