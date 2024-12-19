library(tidyverse)
library(corrr)


exp.list <- readRDS('~/datasets/genomes_assemblies_and_genes/2023_obiol-experiment_assemblies-and-data/quantification/species52/quantification_tmm_species25.rds')
annotation <- read_tsv('~/datasets/genomes_assemblies_and_genes/2023_obiol-experiment_assemblies-and-data/annotation/species52/species_52.emapper.annotations.gz')
corr_genes_to_eggnog <- read_tsv('~/projects/gene_explorer/data/gene_data/linkage_eggnogdb/correspondence_annotations.tsv') %>% 
    distinct(Preferred_name, .keep_all = T)

annotation.sel <- annotation %>% 
    filter(str_detect(`#query`, 'EP00616|EP00893|EP00908'))
    
oneexp <- exp.list$Jul17 
    

annotation.sel <- annotation %>% 
    filter(str_detect(`#query`, 'EP00908'))

tpm_table <- oneexp %>% 
    filter(EukProt_ID == 'EP00908') %>% 
    select(Name, Sample, Abundance)


calculate_correlations <- function(tpm_table, annotation){
    
    geomean <- tpm_table %>% 
        group_by(Sample) %>% 
        summarize( geomean = exp(mean(log( Abundance + 0.1))))
    
    denom.genes <- annotation %>% 
        filter(Preferred_name %in% corr_genes_to_eggnog$Preferred_name) %>% 
        pull(`#query`)
    
    correlations_possible_denoms <- tpm_table %>% 
        filter(Name %in% denom.genes) %>% 
        pivot_wider( names_from = 'Name', values_from = 'Abundance', id_cols = Sample) %>% 
        bind_cols(geomean %>% select(-Sample) ) %>% 
        correlate(method = 'spearman')
    
    random.genes <- tpm_table %>% 
        filter(!Name %in% denom.genes ) %>% 
        pull(Name) %>% 
        unique() %>% 
        sample(nrow(correlations_possible_denoms) - 1)
    
    correlations_random_genes <- tpm_table %>% 
        filter(Name %in% random.genes) %>% 
        pivot_wider( names_from = 'Name', values_from = 'Abundance', id_cols = Sample) %>% 
        bind_cols(geomean %>% select(-Sample) ) %>% 
        correlate(method = 'spearman')
    
    function.denoms <- correlations_possible_denoms %>% 
        select(term, geomean) %>% 
        filter(term != 'geomean') %>% 
        left_join(annotation, by = c('term' = "#query")) %>% 
        left_join(corr_genes_to_eggnog, by = "Preferred_name") 
    
    overall <- data.frame( denoms = correlations_possible_denoms$geomean,
                                random = correlations_random_genes$geomean) %>% 
        pivot_longer(names_to = 'type', values_to = 'corr', cols = c(denoms, random)) 
    
    return( list( functions.denom = function.denoms,
                  overall = overall))
}

allinfo_merged <-  exp.list %>%
    bind_rows(.id = 'experiment') %>%
    group_by(experiment, EukProt_ID) %>%
    nest() %>%
    mutate(correlations = map(data,
                              ~calculate_correlations(.x, annotation = annotation)))
    select(-data) %>% 
    mutate( function.denoms = map(correlations, ~ .x$functions.denom),
            overall = map(correlations, ~ .x$overall)) %>% 
    select(-correlations)
            
    
allinfo_merged %>% 
    unnest(function.denoms) %>% 
    ggplot( aes( geomean)) + 
    geom_histogram(aes(fill = EukProt_ID)) + 
    facet_wrap(~label) + 
    ggtitle('Correlations with geometric mean for each species') + 
    xlab('Correlation') + 
    theme_minimal()

allinfo_merged %>% 
    unnest(function.denoms) %>% 
    filter(label == 'Ribosomal protein') %>% 
    ggplot(aes(geomean)) + 
    geom_histogram(aes(fill = EukProt_ID)) + 
    facet_wrap(~houseg)


allinfo_merged %>% 
    unnest(function.denoms) %>% 
    filter(label == 'Other') %>% 
    ggplot(aes(geomean)) + 
    geom_histogram(aes(fill = EukProt_ID)) + 
    facet_wrap(~houseg)

phaeocystis <- calculate_correlations(  tpm_table = oneexp %>% 
                                            filter(EukProt_ID == 'EP00908') %>% 
                                            select(Name, Sample, Abundance),
                                        annotation =  annotation %>% 
                                            filter(str_detect(`#query`, 'EP00908')))


rhizochromulina <- calculate_correlations(  tpm_table = oneexp %>% 
                                                filter(EukProt_ID == 'EP00616') %>% 
                                                select(Name, Sample, Abundance),
                                            annotation =  annotation %>% 
                                                filter(str_detect(`#query`, 'EP00616')))

MAST <- calculate_correlations(  tpm_table = oneexp %>% 
                                     filter(EukProt_ID == 'MAST-8B-sp1') %>% 
                                     select(Name, Sample, Abundance),
                                 annotation =  annotation %>% 
                                     filter(str_detect(`#query`, 'MAST-8B-sp1')))

acanthoecidae <- calculate_correlations(  tpm_table = oneexp %>% 
                                              filter(EukProt_ID == 'EP00036') %>% 
                                              select(Name, Sample, Abundance),
                                          annotation =  annotation %>% 
                                              filter(str_detect(`#query`, 'EP00036')))

isochrysidales <-  calculate_correlations(  tpm_table = oneexp %>% 
                                                filter(EukProt_ID == 'EP00312') %>% 
                                                select(Name, Sample, Abundance),
                                            annotation =  annotation %>% 
                                                filter(str_detect(`#query`, 'EP00312')))

emiliana <-  calculate_correlations(  tpm_table = oneexp %>% 
                                          filter(EukProt_ID == 'EP00314') %>% 
                                          select(Name, Sample, Abundance),
                                      annotation =  annotation %>% 
                                          filter(str_detect(`#query`, 'EP00314')))

triparma <-  calculate_correlations(  tpm_table = oneexp %>% 
                                          filter(EukProt_ID == 'EP00518') %>% 
                                          select(Name, Sample, Abundance),
                                      annotation =  annotation %>% 
                                          filter(str_detect(`#query`, 'EP00518')))


MAST04 <-  calculate_correlations(  tpm_table = oneexp %>% 
                                          filter(EukProt_ID == 'EP00942') %>% 
                                          select(Name, Sample, Abundance),
                                      annotation =  annotation %>% 
                                          filter(str_detect(`#query`, 'EP00942')))


agrupated_dataset <- list( phaeo = phaeocystis, 
      emiliana =  emiliana, 
      triparma = triparma, 
      MAST04 = MAST04,
      isochrysidales = isochrysidales, 
      MAST8 = MAST, 
      rhizochromulina = rhizochromulina) %>% 
    map( ~ .x$functions.denom$data) %>% 
    bind_rows(.id = 'species')


agrupated_dataset %>% 
             # select(geomean, label)) %>% 
    ggplot( aes( geomean)) + 
    geom_histogram(aes(fill = species)) + 
    facet_wrap(~label) + 
    ggtitle('Correlations with geometric mean for each species') + 
    xlab('Correlation') + 
    theme_minimal()
    
agrupated_dataset %>% 
    filter(label == 'Ribosomal protein') %>% 
    ggplot(aes(geomean)) + 
    geom_histogram(aes(fill = species)) + 
    facet_wrap(~houseg)


agrupated_dataset %>% 
    filter(label == 'Other') %>% 
    ggplot(aes(geomean)) + 
    geom_histogram(aes(fill = species)) + 
    facet_wrap(~houseg)

