## LIBRARIES
require( "data.table" )
require( "ggplot2" )

source( "figS1F_maker.Rl" )

## FUNCTIONS
make_corr_ggdf <- function( bws1, bws2, cores ){

    dflist <- mclapply( names( bws1 ), function( bw1_name ){

        bw1 <- bws1[[ bw1_name ]]
        df1 <- as.data.frame( bw1 )

        tmp1_df <- data.frame()
        for( bw2_name in names( bws2 ) ){

            bw2 <- bws2[[ bw2_name ]]
            df2 <- as.data.frame( bw2 )

            tmp2_df <- cbind( df1, df2 )
            names( tmp2_df ) <- c( paste0( names( df1 ), '1' ), paste0( names( df1 ), '2' ) )
            tmp2_df$name1 <- bw1_name
            tmp2_df$name2 <- bw2_name

            tmp2_df$corr <- round( cor( x = df1$score, y = df2$score, method = "spearman" ), digits = 3 )
            
            tmp1_df <- rbind( tmp1_df, tmp2_df )
            
        }

        return( tmp1_df )
    
    })

    df <- do.call( rbind, dflist )

    return( df )
    
}

add_info_on_ggdf <- function( ggdf ){

    ggdf$fraction1 <- gsub( ".*_", '', ggdf$name1 ) 
    ggdf$protocol1 <- gsub( "_.*", '', gsub( ".*_4f", '4f', ggdf$name1 ) )
    ggdf$sample1 <- gsub( "_.*", '', ggdf$name1 )

    ggdf$fraction2 <- gsub( ".*_", '', ggdf$name2 ) 
    ggdf$protocol2 <- gsub( "_.*", '', gsub( ".*_4f", '4f', ggdf$name2 ) )
    ggdf$sample2 <- gsub( "_.*", '', ggdf$name2 )

    return( ggdf )
    
}

## INPUT
bws_files <- c(

    "path/to/4f-C002_S2L__genome-coverage.bw",
    "path/to/4f-C002_S2S__genome-coverage.bw",
    "path/to/4f-C002_S3__genome-coverage.bw",
    "path/to/4f-C002_S4__genome-coverage.bw", ## C002-r1_4f

    "path/to/4f-C004-r1_S2L__genome-coverage.bw",
    "path/to/4f-C004-r1_S2S__genome-coverage.bw",
    "path/to/4f-C004-r1_S3__genome-coverage.bw",
    "path/to/4f-C004-r1_S4__genome-coverage.bw", ## C004-r1_4f
    
    "path/to/4f-C004-r2_S2L__genome-coverage.bw",
    "path/to/4f-C004-r2_S2S__genome-coverage.bw",
    "path/to/4f-C004-r2_S3__genome-coverage.bw",
    "path/to/4f-C004-r2_S4__genome-coverage.bw", ## C004-r2_4f

    "path/to/10Kh-C002-r1_S2__genome-coverage.bw",
    "path/to/10Kh-C002-r1_S3__genome-coverage.bw",
    "path/to/10Kh-C002-r1_S4__genome-coverage.bw", ## C002-r1_10Kh

    "path/to/10Kh-C004-r1_S2L__genome-coverage.bw",
    "path/to/10Kh-C004-r1_S2S__genome-coverage.bw",
    "path/to/10Kh-C004-r1_S3__genome-coverage.bw",
    "path/to/10Kh-C004-r1_S4__genome-coverage.bw", ## C004-r1_10Kh

    "path/to/10Kh-C004-r2_S2L__genome-coverage.bw",
    "path/to/10Kh-C004-r2_S2S__genome-coverage.bw",
    "path/to/10Kh-C004-r2_S3__genome-coverage.bw",
    "path/to/10Kh-C004-r2_S4__genome-coverage.bw", ## C004-r2_10Kh

)

bws_names <- c(
    
    "C002-r1_4f_S2L",
    "C002-r1_4f_S2S",
    "C002-r1_4f_S3",
    "C002-r1_4f_S4",
    
    "C004-r1_4f_S2L",
    "C004-r1_4f_S2S",
    "C004-r1_4f_S3",
    "C004-r1_4f_S4",
    
    "C004-r2_4f_S2L",
    "C004-r2_4f_S2S",
    "C004-r2_4f_S3",
    "C004-r2_4f_S4",

    "C002-r1_10Kh_S2",
    "C002-r1_10Kh_S3",
    "C002-r1_10Kh_S4",
    
    "C004-r1_10Kh_S2L",
    "C004-r1_10Kh_S2S",
    "C004-r1_10Kh_S3",
    "C004-r1_10Kh_S4",
    
    "C004-r2_10Kh_S2",
    "C004-r2_10Kh_S3",
    "C004-r2_10Kh_S4"

)

binsize <- 50 * 10^3
genome <- "hg38"
bedgenome <- "/path/to/chromsize.bed"

cores <- 6

## PROCESSING
### Import SAMMY values
bin_list <- bins_calculator( bedgenome, binsize, genome )

bws <- import_and_rebin__bw(
    files = bws_files,
    bin_list = bin_list,
    genome = genome,
    names = bws_names,
    cores = cores
)

f4_bws <- bws[ grepl( "_4f_", names( bws ) ) ]
k10_bws <- bws[ grepl( "_10Kh_", names( bws ) ) ]

### Make the ggdf for the correlation plot
full_f4_ggdf <- make_corr_ggdf( bws1 = f4_bws, bws2 = f4_bws, cores = cores )
full_k10_ggdf <- make_corr_ggdf( bws1 = k10_bws, bws2 = k10_bws, cores = cores )
full_f4vk10_ggdf <- make_corr_ggdf( bws1 = f4_bws, bws2 = k10_bws, cores = cores )

full_f4_ggdf <- add_info_on_ggdf( full_f4_ggdf )
full_k10_ggdf <- add_info_on_ggdf( full_k10_ggdf )
full_f4vk10_ggdf <- add_info_on_ggdf( full_f4vk10_ggdf )

### Reduce ggdf size and keep only data interesting for the plot
comparison_keep <- c( "X: C002-r1\nvs\nY: C004-r1", "X: C002-r1\nvs\nY: C004-r2", "X: C004-r1\nvs\nY: C004-r2" )

#### 4f
f4_ggdf <- full_f4_ggdf[
    which( full_f4_ggdf$fraction1 == full_f4_ggdf$fraction2 &
           full_f4_ggdf$sample1 != full_f4_ggdf$sample2 ),
    ]
f4_ggdf$comparison <- paste0( "X: ", f4_ggdf$sample1, "\nvs\nY: ", f4_ggdf$sample2 )
f4_ggdf <- f4_ggdf[ which( f4_ggdf$comparison %in% comparison_keep ), ]

#### 4f10K
k10_ggdf <- full_k10_ggdf[
    which( full_k10_ggdf$fraction1 == full_k10_ggdf$fraction2 &
           full_k10_ggdf$sample1 != full_k10_ggdf$sample2 ),
    ]
k10_ggdf$comparison <- paste0( "X: ", k10_ggdf$sample1, "\nvs\nY: ", k10_ggdf$sample2 )
k10_ggdf <- k10_ggdf[ which( k10_ggdf$comparison %in% comparison_keep ), ]

#### 4f vs 4f10K
f4vk10_ggdf <- full_f4vk10_ggdf[
    which( full_f4vk10_ggdf$fraction1 == full_f4vk10_ggdf$fraction2 &
           full_f4vk10_ggdf$sample1 == full_f4vk10_ggdf$sample2 &
           full_f4vk10_ggdf$protocol1 == "4f" & full_f4vk10_ggdf$protocol2 == "4f10K" ),
    ]
f4vk10_ggdf$comparison <- paste0( "X: ", f4vk10_ggdf$sample1, '-', f4vk10_ggdf$protocol1, "\nvs\nY: ", f4vk10_ggdf$sample2, '-', f4vk10_ggdf$protocol2 )

f4_ggdf$fraction1 <- factor( f4_ggdf$fraction1, levels = c( "S2", "S2S", "S2L", "S3", "S4" ) )
k10_ggdf$fraction1 <- factor( k10_ggdf$fraction1, levels = c( "S2", "S2S", "S2L", "S3", "S4" ) )
f4vk10_ggdf$fraction1 <- factor( f4vk10_ggdf$fraction1, levels = c( "S2", "S2S", "S2L", "S3", "S4" ) )


## OUTPUT
### 4f
scatter <- ggplot( f4_ggdf, aes( x = log10( score1 + 1 ), y = log10( score2 + 1 ) ) ) +
    geom_point( alpha = 0.01, size = 0.1 ) +
    geom_smooth( method = "lm" ) +
    geom_text( aes( label = paste0( "Spearman: ", corr ), x = 0.5, y = 1 ) ) +
    facet_grid( fraction1 ~ comparison ) +
    ylim( c( 0, 1 ) ) +
    xlim( c( 0, 1 ) ) +
    xlab( "log10( RPKM + 1 )" ) +
    ylab( "log10( RPKM + 1 )" ) +
    ggtitle( "4fvs4f" ) +
    theme_bw()
ggsave( "./figS1F.pdf", width = 10, height = 10 )
