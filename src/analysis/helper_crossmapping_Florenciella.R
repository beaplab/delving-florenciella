library(tidyverse)
library(broom)


crossmap_w_RCC1007 <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_RCC1007.tsv')
crossmap_w_RCC1587 <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_RCC1587.tsv')
crossmap_w_pelagomonas <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_Pelagomonas.tsv')
crossmap_w_RCC1693 <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_RCC1693.tsv')



pull(crossmap_w_RCC1007, category_cross) %>% table() %>% as_tibble()
pull(crossmap_w_RCC1587, category_cross) %>% table() %>% as_tibble()
pull(crossmap_w_pelagomonas, category_cross) %>% table() %>% as_tibble()
pull(crossmap_w_RCC1693, category_cross) %>% table() %>% as_tibble()

eggnog <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')
    

crossmap_w_RCC1007 %>% 
    filter(category_cross == 'large crossmapping (>50%)') %>% 
    left_join(eggnog, by = c('refname.x' = 'gene_name')) %>% 
    pull(Description) %>% table() %>% 
    as_tibble() %>% 
    arrange(-n)

crossmap_w_RCC1587 %>% 
    filter(category_cross == 'large crossmapping (>50%)') %>% 
    left_join(eggnog, by = c('refname.x' = 'gene_name')) %>% 
    pull(Description) %>% table() %>% 
    as_tibble() %>% 
    arrange(-n)

allcomparisons <- list( RCC1007 = crossmap_w_RCC1007,
      RCC1587 = crossmap_w_RCC1587, 
      RCC1693 = crossmap_w_RCC1693, 
      Pelagomonas = crossmap_w_pelagomonas) %>% 
    bind_rows(.id = 'comparison') %>% 
    left_join(eggnog, by = c('refname.x' = 'gene_name'))

allcomparisons %>% 
    rename( reads_nocrossmap = `FALSE`, 
            reads_crossmap = `TRUE`, 
            gene_name = refname.x) %>% 
    write_tsv('data/statistics/crossmapping_florenciella_against_all_transcriptomes.tsv')


allcomparisons %>% 
    filter(category_cross == 'large crossmapping (>50%)')  %>% 
    group_by(refname.x) %>% 
    filter( n() >= 3) %>% 
    mutate( total = n()) %>% 
    ungroup() %>% 
    arrange( comparison, total) %>% 
    mutate( refname.x = fct_inorder(refname.x)) %>% 
    ggplot( aes(y = refname.x, x = 1)) + 
    geom_text(aes(label = Description))


comparison.gg <- allcomparisons %>% 
    mutate( category_cross = factor(category_cross, levels = c("no crossmapping",
                                                               "small amount of crossmapping (<10%)",
                                                               "intermediate situation (10-50%)", 
                                                               "large crossmapping (>50%)", 
                                                               "only reads crossmapped!")),
            comparison = factor(comparison, levels = rev(c("RCC1587",
                                                       "RCC1693", 
                                                       "RCC1007", 
                                                       "Pelagomonas")))
            
            ) %>% 
    group_by(comparison,category_cross) %>% 
    summarize(count = n()) %>% 
    ggplot( aes(count, comparison) ) + 
    geom_bar(stat = 'identity', aes(fill = category_cross)) + 
    scale_fill_brewer() + 
    theme_bw()

ggsave(filename = 'results/figures/crossmapping_florenciella-multiple.pdf',
       plot = comparison.gg,
       width = 7, height = 3)

# how many of the genes are still present in the TPM filtered 

tpm_filt <- read_tsv('data/quantification/alignment/tpm_all_filtered.tsv')

# we will consider it against the RCC1007 since it is the outermost transcriptome
# from this transcriptome https://github.com/beaplab/transcriptome_metaT_quantification/issues/1

crossmap_w_RCC1007 %>% 
    filter(category_cross == 'large crossmapping (>50%)')   %>% 
    filter(refname.x %in% unique(tpm_filt$gene_name)) 


genes_w_large_cross <- allcomparisons %>% 
    filter(category_cross == 'large crossmapping (>50%)')  %>% 
    filter(refname.x %in% unique(tpm_filt$gene_name)) %>% 
    group_by(refname.x) %>% 
    filter( n() >= 3) %>% 
    pull(refname.x) %>% unique()

length(genes_w_large_cross) / length(unique(tpm_filt$gene_name) )


# Relationship between crossmapping and the abundance of the gene 

tpm_filt %>% 
    group_by(gene_name) %>% 
    summarize( total = sum(TPM)) %>% 
    filter(total < 100000) %>%  
    # left_join(eggnog, by = 'gene_name') %>% 
    # View()
    left_join(crossmap_w_RCC1007, by = c('gene_name' ='refname.x' )) %>% 
    ggplot( aes( category_cross, total)) + 
    geom_point()


# Plot with abundance instead of count genes ------------------------------

datasetab <- allcomparisons %>% 
    mutate( category_cross = factor(category_cross, levels = c("no crossmapping",
                                                               "small amount of crossmapping (<10%)",
                                                               "intermediate situation (10-50%)", 
                                                               "large crossmapping (>50%)", 
                                                               "only reads crossmapped!")),
            comparison = factor(comparison, levels = rev(c("RCC1587",
                                                           "RCC1693", 
                                                           "RCC1007", 
                                                           "Pelagomonas")))
            
    ) %>% 
    rename( reads_nocrossmap = `FALSE`, 
            reads_crossmap = `TRUE`,
            gene_name = refname.x) 


datasetab %>% 
    group_by(comparison, category_cross) %>% 
    summarize( total_nocrossmap = sum(reads_nocrossmap, na.rm = T),
               total_crossmap = sum(reads_crossmap, na.rm = T)) %>% 
    pivot_longer( names_to = 'type', values_to = 'reads', cols = c(total_nocrossmap, total_crossmap)) %>% 
    group_by(comparison) %>% 
    mutate( relab_category = reads / sum(reads)) %>% 
               # relab_crossmap = total_crossmap / sum(total_nocrossmap, total_crossmap)) %>% 
    ggplot( aes(relab_category, comparison) ) + 
    geom_bar(stat = 'identity', aes(color = type, fill = category_cross)) + 
    xlab('Total relative abundance (%)') + 
    scale_fill_brewer() + 
    theme_bw() 


