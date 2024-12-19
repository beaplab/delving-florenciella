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

profile_by_salmon <- read_tsv(c(#'data/quantification/mapping/Numreads_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv',
                                'data/quantification/mapping/Numreads_EP00618_Florenciella_parvula_quant-over_2021_tara_polar.tsv')
                                )

profile_salmon_long <- pivot_longer(profile_by_salmon, names_to = 'gene', values_to = 'count_salmon', cols = -sample)

allcounts <- left_join(profile_salmon_long, profile_df, by = c('sample', 'gene'))

comparison <- allcounts %>% 
    filter( (count_bwa != 0) & count_salmon != 0)

corr_genes_to_eggnog <- read_tsv('~/projects/gene_explorer/data/gene_data/linkage_eggnogdb/correspondence_annotations.tsv') %>% 
    #TODO do it more elegantly
    distinct(Preferred_name, .keep_all = T)

annotation <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula.emapper.annotations', 
                       comment = "##") %>% 
    rename( gene_name = `#query`) %>% 
    left_join(corr_genes_to_eggnog, by = 'Preferred_name')
 

# General plot ------------------------------------------------------------

stats_bygene <- comparison %>% 
    group_by(gene) %>% 
    summarize(ocurrence_salmon = sum(count_salmon != 0),
              ocurrence_bwa = sum(count_bwa != 0),
              median_count_salmon = median(count_salmon),
              median_count_bwa = median(count_bwa)
              ) %>% 
    arrange(-median_count_bwa) %>% 
    mutate( rank_bwa = 1:nrow(.))


stats_bygene %>% 
    # filter(rank_bwa <= 50) %>%
    filter(rank_bwa >= 50 & rank_bwa <= 100) %>%
    # filter(rank_bwa >= 100 & rank_bwa <= 200) %>% 
    pivot_longer( names_to = 'method',
                  values_to = 'median_count',
                  cols = contains("median_count")) %>% 
    ggplot( aes(rank_bwa, median_count)) + 
    geom_bar(stat = 'identity', aes(fill = method), position = 'dodge') + 
    ylab('Mean count') 

genes.dist.df <- comparison %>% 
    group_by(gene) %>% 
    summarize(total = sum(count_bwa)) %>% 
    arrange(-total) %>% 
    filter(total > 100) %>% 
    mutate( differentiation = case_when( total >= 900000 ~ 'most crossmapped genes',
                                         total < 900000 & total >= 30000 ~ 'middle ground',
                                         total < 300000 ~ 'the rest'))

genes.dist.df$differentiation %>% table()
    


library(ggpubr)

topones <- genes.dist.df %>% filter(differentiation == 'most crossmapped genes')
middle <- genes.dist.df %>% filter(differentiation == 'middle ground')
rest <- genes.dist.df %>% filter(differentiation == 'the rest')

ribosomes <- annotation %>% 
    filter(label == 'Ribosomal protein') %>% 
    pull(gene_name)


plot_relationship <- function(df){
    
    df %>% 
    ggplot( aes(count_bwa, count_salmon)) + 
        geom_point() + 
        geom_smooth() + 
        geom_abline(  slope = 1, intercept = 0, color = 'red') +
        stat_cor(method = 'pearson') + 
        facet_wrap(~str_c(gene,
                          '; rank:',
                          rank_bwa, '; ',
                          Preferred_name), scales = 'free')
        # scale_y_log10() + 
        # scale_x_log10() + 
        
}

comparison_w_info <-  comparison %>% 
    left_join(stats_bygene, by = 'gene') %>% 
    left_join(annotation, by = c('gene' = 'gene_name')) 
    
comparison_w_info %>% 
    filter(gene %in% topones$gene ) %>%
    plot_relationship() 

comparison_w_info %>% 
    filter(gene %in% sample(middle$gene, 20) ) %>% 
    plot_relationship() 

comparison_w_info %>% 
    filter(gene %in% sample(rest$gene, 20) ) %>% 
    plot_relationship() 

comparison_w_info %>% 
    filter(gene %in% sample(ribosomes, 20) ) %>% 
    plot_relationship() + 
    stat_cor(method = 'pearson') + 
    facet_wrap(~gene, scales = 'free')


comparison %>% 
    filter(gene %in% topones$gene ) %>%
    # filter(gene %in% sample(ribosomes, 10) ) %>% 
    left_join(stats_bygene, by = 'gene') %>% 
    left_join(annotation, by = c('gene' = 'gene_name')) %>% 
    ggplot( aes(count_bwa, count_salmon)) + 
        geom_point() + 
        geom_smooth() + 
    geom_abline(  slope = 1, intercept = 0, color = 'red') +
    stat_cor(method = 'pearson') + 
    facet_wrap(~str_c(gene,'; rank:', rank_bwa, '; ', Preferred_name), scales = 'free')



# General stats --------------------------------------------------------

# BWA presents half of the reads than Salmon
comparison %>% 
    filter(!gene %in% topones$gene) %>%
    summarize( total_bwa = sum(count_bwa), 
               total_salmon = sum(count_salmon), 
               ratio = total_bwa / total_salmon)
    
corrs <- comparison %>% 
    group_by(gene) %>% 
    summarize( correlation = cor(count_bwa, count_salmon))

ratios <- comparison %>% 
    mutate( ratio = count_bwa / count_salmon) %>% 
    group_by(gene) %>% 
    summarize( m_ratio = median(ratio),
               total_ratio = sum(count_bwa) / sum(count_salmon))

library(patchwork)

corrs %>% 
    filter(correlation < 0.9) %>% 
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

ocurrence_genes <- left_join(profile_df, profile_salmon_long, by = c('sample', 'gene'))

gene_selection <- ocurrence_genes %>% 
    group_by(gene) %>% 
    summarize( ocurrence = sum(count_bwa != 0)) %>% 
    filter( ocurrence >= (114 * 0.9)) 

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




