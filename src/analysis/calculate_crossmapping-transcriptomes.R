#!/usr/bin/env Rscript

library(tidyverse)
library(argparser)
library(Rsamtools) %>% suppressPackageStartupMessages()

### 
# Given two BAM directories compare all the read fates for the shared BAM names
# and output a summary
###

# Create a parser
p <- arg_parser("Evaluates crossmapping between two mapped transcriptomes against the same samples")

# Add command line arguments
p <- add_argument(p, "--ref-path",
                  help="the path to the dir with the reference BAMs",
                  type="character",
                  default = '')
p <- add_argument(p, "--query-path", help="the path to the dir with the query BAMs", default='')
p <- add_argument(p, "--out", help="the name of the outputfile, specifying the directory", default= '')

argv <- parse_args(p)

path_ref <- argv$ref_path
path_query <- argv$query_path

bams_transcriptome_ref <- list.files(path = path_ref,
                                  recursive = TRUE,
                                  pattern = '.bam',
                                  full.names = T)

bams_transcriptome_query <- list.files(path = path_query,
                                  recursive = TRUE,
                                  pattern = '.bam',
                                  full.names = T)


samples_shared <- intersect( basename(bams_transcriptome_ref), 
                             basename(bams_transcriptome_query) )


print( str_c("There are ", length(samples_shared), " shared samples here."))

extract_reads_and_query_names <- function(path){
    # 
    bam = scanBam(path)
    
    tibble( qname = bam[[1]]$qname,
            refname = bam[[1]]$rname) %>% 
        return()
    
}


summary_crossmapping <- function(sample){
    
    bamreads1 <- extract_reads_and_query_names(str_c(path_ref, sample, sep = '/'))
    bamreads2 <- extract_reads_and_query_names(str_c(path_query, sample, sep = '/'))
    
    
    comparison <-  bamreads1 %>% 
        left_join(bamreads2, by = 'qname') %>%  
        mutate( crossmapping = ifelse(!is.na(refname.y), TRUE, FALSE)) %>% 
        mutate( refname.x = as.character(refname.x), 
                refname.y = as.character(refname.y))
    
    statistic_crossmapping <-  comparison  %>% 
        group_by( refname.x, crossmapping) %>% 
        summarize( count_crossmapped_reads = n(), 
                   genes_cross_set = list(unique(refname.y))) %>% 
        ungroup() %>% 
        arrange(refname.x) 
    
    return(statistic_crossmapping)
}


print("Calculating the summary") 

allsample_comparison_df <- tibble( sample = samples_shared) %>% 
    mutate( statistic = map(sample, ~summary_crossmapping(.x)))  %>% 
    unnest(statistic)

print("Outputting final statistics")

classification_df <- allsample_comparison_df %>% 
    group_by(refname.x, crossmapping) %>% 
    summarize( totals = sum(count_crossmapped_reads)) %>% 
    pivot_wider( names_from = crossmapping,
                 values_from = totals,
                 id_cols = refname.x) %>% 
    mutate( ratio = ifelse(is.na(`TRUE`), 0, `TRUE` / `FALSE`), 
            category_cross =  case_when(
                ratio == 0 ~ 'no crossmapping', 
                ratio <= 0.1 ~ 'small amount of crossmapping (<10%)',
                ratio >= 0.5 ~ 'large crossmapping (>50%)', 
                is.na(ratio) ~ 'only reads crossmapped!',
                TRUE ~ 'intermediate situation (10-50%)'
            )
    ) 



print("The genes present the crossmapping as follows:")

print(
classification_df %>% 
    pull(category_cross) %>% 
    table()
)

write_rds(x = allsample_comparison_df, 
        file = str_c( argv$out, '_df.rds'))

write_tsv(x = classification_df, 
        file = str_c( argv$out, '.tsv'))


    