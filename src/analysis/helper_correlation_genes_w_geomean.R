library(tidyverse)
library(corrr)

# Data importing
florenciella <- read_tsv('data/quantification/mapping/TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv') 

annotation <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula.emapper.annotations', 
                       comment = "##") %>% 
    rename( gene_name = `#query`) %>% 
    left_join(corr_genes_to_eggnog, by = 'Preferred_name')

corr_genes_to_eggnog <- read_tsv('~/projects/gene_explorer/data/gene_data/linkage_eggnogdb/correspondence_annotations.tsv') %>% 
    #TODO do it more elegantly
    distinct(Preferred_name, .keep_all = T)

quant.transcriptome <- florenciella %>% 
    pivot_longer( names_to = 'gene', values_to = 'count', cols = -sample) %>% 
    group_by(gene) %>% 
    filter(!sum(count == 0) == n()) 

geomean.transcriptome <- quant.transcriptome %>% 
    group_by(sample) %>% 
    summarize( geomean = exp(mean(log( count + 0.1))))

denom.genes <- annotation %>% 
    filter(Preferred_name %in% corr_genes_to_eggnog$Preferred_name) %>% 
    pull(gene_name)

correlations_possible_denoms <- quant.transcriptome %>% 
    filter(gene %in% denom.genes) %>% 
    pivot_wider( names_from = 'gene', values_from = 'count', id_cols = sample) %>% 
    bind_cols(geomean.transcriptome %>% select(-sample) ) %>% 
    correlate(method = 'pearson')


# sampled wiwth genes presenting at least a total TPM of 2. 
random.genes <- quant.transcriptome %>% 
    group_by(gene) %>% 
    summarize( total = sum(count)) %>% 
    filter(total > 2, 
           !gene %in% denom.genes ) %>% 
    pull(gene) %>% unique()

correlation_random_genes <- quant.transcriptome %>% 
    filter(gene %in% sample(size = 800, random.genes)) %>%
    filter(gene %in% random.genes) %>%
    pivot_wider( names_from = 'gene',
                 values_from = 'count',
                 id_cols = sample) %>% 
    .[,-1] %>% 
    map(~correlate(.x,
                   geomean.transcriptome %>% select(-sample), 
                   method = 'pearson')) %>% 
    bind_rows( .id = 'term')

allcors <- correlations_possible_denoms %>% 
    select(term, geomean) %>% 
    bind_rows( correlation_random_genes) %>% 
    left_join(annotation %>% select(gene_name, category, label), 
              by = c('term' = 'gene_name')) %>% 
    rename( corr_w_geomean = geomean) 

median.counts <- quant.transcriptome %>%
    group_by(gene) %>%
    summarize(median.count = median(count),
              occurrence = sum(count > 0))


# From here on visualization ----------------------------------------------


allcors %>% 
    left_join(median.counts, by = c('term' = 'gene')) %>%
    # filter(!is.na(label)) %>% 
    ggplot( aes(corr_w_geomean, label)) + 
    # geom_jitter(aes(color = label), show.legend = F) + 
    geom_jitter(aes(size = occurrence, color = occurrence)) + 
    geom_violin(draw_quantiles = 0.5, alpha = 0.9) +
    xlab('Correlation') + 
    theme_bw()

# selection over 0.6 correlation
allcors %>% 
    filter(corr_w_geomean > 0.6) %>% 
    left_join(median.counts, by = c('term' = 'gene')) %>% 
    ggplot( aes(label, corr_w_geomean)) + 
    geom_jitter(aes(size = occurrence, color = occurrence)) + 
    geom_violin(draw_quantiles = 0.5, alpha = 0.9) +
    xlab('Correlation')
    

relationships_w_correlation <- allcors %>% 
    left_join(annotation, by = c('term' = 'gene_name')) %>% 
    left_join(median.counts, by = c('term' = 'gene'))

relationships_w_correlation %>% 
    # filter(!is.na(label.x)) %>%
    ggplot(aes(occurrence, corr_w_geomean)) + 
    geom_point(aes(color = label.x, size = log10(median.count + 1))) + 
    ggrepel::geom_text_repel(data = relationships_w_correlation %>%
                                 filter(!is.na(label.x)) %>%
                                 filter(occurrence > 150),
                             aes(label = Preferred_name, color = label.x)) +
    theme_bw() + 
    scale_size_continuous(name = "Median count") + 
    ylab('Correlation') + 
    ggtitle('Correlation between putative housekeeping genes') 


relationships_w_correlation %>% 
    filter(is.na(label.x)) %>% 
    ggplot(aes(occurrence, corr_w_geomean)) + 
    geom_point(aes(size = log10(median.count + 1))) + 
    ggrepel::geom_text_repel(data = relationships_w_correlation %>%
                                 filter(is.na(label.x)) %>%
                                 filter(occurrence > 150),
                             aes(label = Preferred_name)) +
    theme_bw() + 
    ggtitle('Only the other genes') + 
    ylab('Correlation')


selection_genes <- read_tsv('~/projects/gene_explorer/data/stats/selection_houseg_denoms_bycorr-07.tsv')

selection_by_degree_and_corr <- annotation %>% filter(Preferred_name %in%
                                           (corr_genes_to_eggnog %>%
                                                filter(houseg %in% selection_genes$houseg) %>% 
                                                pull(Preferred_name))) %>% 
    pull(gene_name)



dataonly <- relationships_w_correlation %>% 
    filter(!is.na(label.x))

dataonly %>% 
    ggplot(aes(occurrence, corr_w_geomean)) + 
    geom_point(aes(color = label.x, size = log10(median.count + 1))) + 
    geom_point(data = dataonly %>% filter(term %in% selection_by_degree_and_corr),
               color = 'black') + 
    ggrepel::geom_text_repel(data = relationships_w_correlation %>%
                                 filter(!is.na(label.x)) %>%
                                 filter(occurrence > 150),
                             aes(label = Preferred_name, color = label.x)) +
    theme_bw() + 
    scale_size_continuous(name = "Median count") + 
    ylab('Correlation') + 
    ggtitle('Correlation between putative housekeeping genes') 




# Overall -----------------------------------------------------------------

# Distribution all correlations 
corr.possible.denoms.annotated <- correlations_possible_denoms %>% 
    pivot_longer(names_to = 'term2', values_to = 'corr', cols = -term) %>% 
    filter(! term %in% 'geomean', !term2 %in% 'geomean', !is.na(corr)) %>% 
    left_join(annotation %>% select(gene_name, category, label), 
              by = c('term' = 'gene_name')) %>% 
    left_join(annotation %>% select(gene_name, category, label), 
              by = c('term2' = 'gene_name'))  

# Only chaperones 
corr.possible.denoms.annotated %>% 
    filter( label.x == 'Chaperone' & label.y == 'Chaperone') %>% 
    ggplot(aes(corr)) + 
    geom_histogram() + 
    theme_bw()

corr.possible.denoms.annotated %>% 
    filter( label.x %in% c('proteasome', 'Ribosomal protein', 'Chaperone') &
                label.y %in% c('proteasome', 'Ribosomal protein', 'Chaperone')) %>% 
    filter(label.x != label.y) %>% 
    ggplot(aes(corr)) + 
    geom_histogram() + 
    theme_bw()

