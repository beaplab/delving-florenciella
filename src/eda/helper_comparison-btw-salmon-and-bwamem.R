library(tidyverse)

profiles <- list.files('data/quantification/alignment/EP00618_Florenciella_parvula/',
                       pattern = 'profile',
                       full.names = T)

profile_df <- tibble( path = profiles, 
        sample  = basename(profiles) %>% str_remove(pattern = '.profile.txt.gz')) %>% 
    mutate( 
        dataset = map(profiles, ~read_tsv(.x, skip = 10, col_names = c('gene', 'count_bwa')))
    ) %>% 
    select(sample,dataset) %>% 
    unnest(dataset)

profile_by_salmon <- read_tsv('data/quantification/Numreads_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv')

profile_salmon_long <- pivot_longer(profile_by_salmon, names_to = 'gene', values_to = 'count_salmon', cols = -sample)


allcounts <- left_join(profile_df, profile_salmon_long, by = c('sample', 'gene'))
comparison <- allcounts %>% 
    filter( (count_bwa != 0) & count_salmon != 0)

# General plot ------------------------------------------------------------

ggplot( comparison, aes(count_bwa, count_salmon)) + 
    geom_hex() + 
    geom_smooth() + 
    geom_abline(  slope = 1, intercept = 0, color = 'red') + 
    scale_y_log10() + 
    scale_x_log10() + 
    coord_equal() 

comparison %>% 
    summarize( total_bwa = sum(count_bwa), 
               total_salmon = sum(count_salmon), 
               ratio = total_bwa / total_salmon)

# General stats --------------------------------------------------------

corrs <- comparison %>% 
    group_by(gene) %>% 
    summarize( correlation = cor(count_bwa, count_salmon))

ratios <- comparison %>% 
    mutate( ratio = count_bwa / count_salmon) %>% 
    group_by(gene) %>% 
    summarize( m_ratio = median(ratio),
               total_ratio = sum(count_bwa) / sum(count_salmon))



corrs %>% 
    ggplot( aes(correlation)) +
    geom_histogram()

ratios %>% 
    filter(m_ratio <= 1) %>% 
    ggplot( aes(m_ratio)) +
    geom_histogram()


# Genes presenting a lot of reads  ----------------------------------------

selection_lotsreads  <- comparison %>% 
    group_by(gene) %>% 
    summarize( total = sum(count_bwa)) %>% 
    filter(total >= 1000)

ratios %>% 
    filter(gene %in% selection_lotsreads$gene) %>% 
    filter(m_ratio <= 1) %>% 
    mutate( putative_housekeeping = ifelse(!gene %in% gene_selection$gene, TRUE, FALSE)) %>% 
    ggplot( aes(m_ratio)) +
    geom_histogram(aes(fill = putative_housekeeping)) + 
    ggtitle( 'Genes presenting >1000 total counts')
    
(selection_lotsreads$total %>% sum() ) / sum(allcounts$count_bwa)


comparison %>% 
    group_by(gene) %>% 
    # at least in 20 samples 
    filter( n() > 20) %>% 
    summarize( total_bwa = sum(count_bwa), 
               total_salmon = sum(count_salmon), 
               ratio = total_bwa / total_salmon) %>% 
    arrange(-ratio) %>% 
    filter(ratio < 1) %>% 
    ggplot( aes(ratio)) +
    geom_histogram()


# Other -------------------------------------------------------------------


# other 
ocurrence_genes <- left_join(profile_df, profile_salmon_long, by = c('sample', 'gene'))

gene_selection <- ocurrence_genes %>% 
    group_by(gene) %>% 
    summarize( ocurrence = sum(count_bwa != 0)) %>% 
    filter( ocurrence <= (114 * 0.9)) 

(selection_lotsreads %>%
        filter(gene %in% gene_selection$gene) %>%
        pull(total) %>% sum()) / sum(allcounts$count_bwa)

comparison %>% 
    filter(gene %in% gene_selection$gene) %>% 
    mutate( ratio = count_bwa / count_salmon) %>% 
    group_by( gene ) %>% 
    summarize( total_bwa = sum(count_bwa),
               total_salmon = sum(count_salmon),
               m_ratio = mean(ratio)) %>% 
    left_join(corrs, by = 'gene') %>% 
    filter(m_ratio <= 1) %>% 
    ggplot( aes(m_ratio, correlation)) + 
    geom_point()



