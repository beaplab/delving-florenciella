library(tidyverse)
library(Rsamtools)
library(UpSetR)

allbams <- list.files(path = 'data/quantification/alignment',
           recursive = TRUE,
           pattern = '.filtered.bam',
           full.names = T)

extract_reads_and_query_names <- function(path){
    
    bam = scanBam(path)
    tibble( qname = bam[[1]]$qname,
            refname = bam[[1]]$rname) %>% 
        return()
    
}

query_reads_df <- tibble(path = allbams) %>% 
    mutate( transcriptome = map_chr(path, ~str_split_1(.x, pattern = '/')[4])) %>% 
    mutate( sample = map_chr(path, ~str_split_1(.x, pattern = '/')[5]))

unique_samples <- query_reads_df$sample %>% unique() 
unique_transcriptomes <- query_reads_df$transcriptome%>% unique() 


samples_w_7_transcriptomes <- query_reads_df %>% 
    group_by(sample) %>% 
    # Let's work only with the samples that present all the 7 transcriptomes
    filter( n() == 7) %>% 
    # ok and for starters only some of them
    ungroup() %>% 
    filter( sample %in% sample( unique(sample), size = 4)) %>%
    mutate( query_readnames = map(path, extract_reads_and_query_names))

comparison_bathy_pelagomonas_phaeocystis <- query_reads_df %>% 
    filter( transcriptome %in% c("EP00893_Bathycoccus_prasinos",
                                 "EP00908_Phaeocystis_cordata",
                                 "Pelagomonas_calceolata")) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    # Let's work only with the samples that present all the 7 transcriptomes
    filter( n() == 3) %>% 
    ungroup() %>% 
    filter( sample %in% sample( unique(sample), size = 4)) 
    # mutate( query_readnames = map(path, extract_reads_and_query_names))


comparison_florenciellas <- query_reads_df %>% 
    filter( transcriptome %in% c("EP00618_Florenciella_parvula",
                                 "EP00619_Florenciella_sp_RCC1007",
                                 "EP00620_Florenciella_sp_RCC1587",
                                 "EP00621_Florenciella_sp_RCC1693") ) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    # Let's work only with the samples that present all the 7 transcriptomes
    filter( n() == 4) %>% 
    ungroup() %>% 
    filter( sample %in% sample( unique(sample), size = 10)) 

    
# Overall distribution ----------------------------------------------------


create_list_by_transcriptome <- function(df){
    
    df %>% 
        select(-path) %>% 
        unnest(query_readnames) %>% 
        split( ~transcriptome) %>% 
        map(~.x %>% pull(qname)) %>% 
        return()
    
}



## Very abundant transcriptomes  -------------------------------------------

upset(nsets = 7,
      
      fromList(
          comparison_bathy_pelagomonas_phaeocystis %>% 
              mutate( query_readnames = map(path, extract_reads_and_query_names))  %>%
              create_list_by_transcriptome()
      ) 
)



## Florenciella transcriptomes ---------------------------------------------

upset(nsets = 7,
      
      fromList(
          comparison_florenciellas %>% 
              mutate( query_readnames = map(path, extract_reads_and_query_names))  %>%
              create_list_by_transcriptome()
      ) 
)


## Florenciella removing the most crossmapping genes -----------------------

crossmapping_parvula <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_Pelagomonas.tsv')

largecross_genes_floren <- crossmapping_parvula %>% 
    filter(category_cross %in% "large crossmapping (>50%)") %>% 
    pull(refname.x)

reads_df_floren <- comparison_florenciellas %>% 
    mutate( query_readnames = map(path, extract_reads_and_query_names)) %>% 
    select(-path) %>% 
    unnest(query_readnames) 

reads_to_not_consider <- reads_df_floren  %>% 
    filter(refname %in% largecross_genes_floren) %>% 
    pull(qname)


upset(nsets = 7,
      
      fromList(
          comparison_florenciellas %>% 
              mutate( query_readnames = map(path, extract_reads_and_query_names))  %>%
              select(-path) %>% 
              unnest(query_readnames) %>% 
              filter(!qname %in% reads_to_not_consider) %>% 
              split( ~transcriptome) %>% 
              map(~.x %>% pull(qname)) 
      ) 
)
    


# Testing read missclassification by gene ---------------------------------

allinfo <- samples_w_7_transcriptomes %>% 
    # filter( sample %in% sample( unique(sample), size = 10)) %>% 
    select(-path) %>% 
    unnest(query_readnames) 
    

count_by_reads <- allinfo %>% 
    group_by( sample, qname) %>% 
    summarize( ntimes = n()) 

crossmapping_df <- allinfo %>% 
    filter( transcriptome %in% 'EP00618_Florenciella_parvula') %>% 
    left_join(count_by_reads, by = c('sample', 'qname')) %>% 
    group_by(refname) %>% 
    summarize( multimapping = mean(ntimes)) 

crossmapping_df %>% 
    ggplot( aes(multimapping)) + 
    geom_histogram() + 
    ggtitle('Histogram of the crossmapping between the 7 Florenciella transcriptomes in Florenciella parvula',
            subtitle = 'Only 10 samples are being considered') + 
    xlab('Number of transcriptomes')

eggnog.df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

crossmapping_df %>% 
    filter(multimapping > 7) %>% 
    left_join(eggnog.df, by = c('refname' = 'gene_name')) %>% 
    View()
