library(tidyverse)
library(corrr)


# evaluation of sample nomenclature: 
tpm_table <- readxl::read_xlsx('data/quantification/statistics/orciraptor_suppdata/Gene_expression/Gene_expression.xlsx') 

# a vs g (attacking vs gliding),
# log2FC, for m.29664_type-3prime_partial_len159_TR10001_c0_g2_i1, should be -0.323

log2(tpm_table[8,3:5] %>% rowMeans() / tpm_table[8,6:8] %>% rowMeans() )
log2(tpm_table[8,3:5] %>% rowMeans() / tpm_table[8,9:11] %>% rowMeans() )
log2(tpm_table[8,6:8] %>% rowMeans() / tpm_table[8,9:11] %>% rowMeans() )


annotation.eggnog <- readxl::read_xlsx('data/quantification/statistics/orciraptor_suppdata/Annotation/Annotation_eggNOG.xlsx')
tpm_table <- readxl::read_xlsx('data/quantification/statistics/orciraptor_suppdata/Gene_expression/Gene_expression.xlsx') %>% 
    select(transcript, V1S4_TPM:V1S3_TPM) %>% 
    pivot_longer(names_to = 'sample', values_to = 'TPM', cols = -transcript)

corr_genes_to_eggnog <- read_tsv('~/projects/gene_explorer/data/gene_data/linkage_eggnogdb/correspondence_annotations.tsv') %>% 
    distinct(Preferred_name, .keep_all = T)

geomean.transcriptome <- tpm_table %>% 
    group_by(sample) %>% 
    summarize( geomean = exp(mean(log( TPM + 0.1))))


denom.genes <- annotation.eggnog %>% 
    filter(Preferred_name %in% corr_genes_to_eggnog$Preferred_name) %>% 
    pull(X.query_name)

correlations_possible_denoms <- tpm_table %>% 
    filter(transcript %in% denom.genes) %>% 
    pivot_wider( names_from = 'transcript', values_from = 'TPM', id_cols = sample) %>% 
    bind_cols(geomean.transcriptome %>% select(-sample) ) %>% 
    correlate(method = 'spearman')

random.genes <- tpm_table %>% 
    filter(!transcript %in% denom.genes ) %>% 
    pull(transcript) %>% 
    unique() %>% 
    sample(261)


correlations_random_genes <- tpm_table %>% 
    filter(transcript %in% random.genes) %>% 
    pivot_wider( names_from = 'transcript', values_from = 'TPM', id_cols = sample) %>% 
    bind_cols(geomean.transcriptome %>% select(-sample) ) %>% 
    correlate(method = 'spearman')


data.frame( denoms = correlations_possible_denoms$geomean,
            random = correlations_random_genes$geomean) %>% 
    pivot_longer(names_to = 'type', values_to = 'corr', cols = c(denoms, random)) %>% 
    ggplot( aes(corr)) + 
    geom_histogram(aes(fill = type))

ggplot(correlations_possible_denoms, aes(geomean)) + geom_histogram()
ggplot(correlations_random_genes, aes(geomean)) + geom_histogram()

