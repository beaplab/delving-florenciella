library(tidyverse)

tpm_filt <- read_tsv('data/quantification/alignment/tpm_all_filtered.tsv')%>% 
    filter(group %in% c('2021_tara_polar', '2012_carradec_tara')) %>% 
    mutate( sample = str_replace(sample, pattern = '0o8', replacement = '0.8')) %>% 
    mutate( sample = str_replace(sample, pattern = 'SUR', 'SRF'),
            sample = str_replace(sample, pattern = 'inf', '>'))

envdir <- '/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_polar_protists/Tara_Oceans_Pangaea_context.rds'
panagea_context_list <- readRDS(envdir)
correspondence_panagea_tarapolar <- read_tsv('/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_polar_protists/correspondence_pangaea_tara_polar.tsv')

iron_df <- panagea_context_list$watercolumn$values %>% 
    as_tibble() 

alliron_df <- iron_df %>% 
    left_join(correspondence_panagea_tarapolar, by = c('sample_id_pangaea' = 'panagea_id')) %>% 
    select(sample_id_pangaea, sample_material,
           run_accession, env_feature, iron_darwin_5m, iron_darwin_5m_std) %>% 
    mutate(sample = ifelse( is.na(run_accession),
                            map_chr(sample_material, 
                                    ~str_split_1(.x, pattern = "_") %>% 
                                    .[c(2,4,3)] %>% 
                                    str_flatten(collapse =  '_')),
           run_accession))  %>% 
    distinct(sample, .keep_all = T)


polarenv <- read_tsv(str_c(envdir,
                          '/tara_polar_protists/envdata_tara_polar_protist.tsv'))

eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

genes_sel <- eggnog_df %>% 
    filter( str_detect( Description, 'ferredoxin|Flavodoxin|peptidase')) %>% 
    select(gene_name, Description) %>% 
    mutate( Description = case_when( str_detect( Description, 'ferredoxin') ~'Ferredoxin',
                                     str_detect( Description, 'Flavodoxin') ~'Flavodoxin',
                                     str_detect( Description, 'peptidase') ~'Peptidase',
                                     TRUE ~ Description
                                     
            ))

tpm_filt %>% 
    filter(gene_name %in% genes_sel$gene_name | gene_name %in% sample(gene_name, 20)) %>% 
    left_join(alliron_df, by = 'sample') %>% 
    left_join(genes_sel, by = 'gene_name') %>% 
    # filter(is.na(iron_darwin_5m)) %>% 
    ggplot( aes(iron_darwin_5m, TPM)) + 
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~Description, scales = 'free_y') + 
    theme_minimal() 


genes_sel2 <- eggnog_df %>% 
    filter( str_detect( Description, 'fructose')) 

tpm_filt %>% 
    filter(gene_name %in% genes_sel2$gene_name ) %>% 
    left_join(genes_sel2, by = 'gene_name') %>% 
    group_by(Description, sample) %>% 
    summarize( total_TPM = sum(TPM)) %>% 
    left_join(alliron_df, by = 'sample') %>% 
    # filter(is.na(iron_darwin_5m)) %>% 
    ggplot( aes(iron_darwin_5m, total_TPM)) + 
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~Description, scales = 'free_y') + 
    theme_minimal() 

