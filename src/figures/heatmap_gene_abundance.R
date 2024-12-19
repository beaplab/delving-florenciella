library(tidyverse)
library(patchwork)


tpm_all <- read_tsv('data/quantification/mapping/tpm_all_filtered.tsv')
eggnog.df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv') %>% 
    distinct(gene_name, .keep_all = T) %>% 
    # TODO what happens with this 
    filter(!is.na(gene_name))
envdata <- read_tsv('/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_Alldata_including_polar_protists/envdata_Tara_All_downloaded_metaT_nisaba.tsv')

source('src/figures/helper_heatmap-tara-samples.R')


matrix <- tpm_all %>%
    filter(group %in% c('2012_carradec_tara', '2021_tara_polar')) %>% 
    mutate( TPM = log10(TPM + 0.001)) %>% 
    select(-group) %>% 
    pivot_wider( names_from = sample, values_from = TPM, id_cols = gene_name) %>% 
    column_to_rownames(var = 'gene_name') %>% 
    as.matrix() 


# General  ----------------------------------------------------------------

peptidases <- read_tsv('/scratch/datasets_symbolic_links/genomes_assemblies_and_genes/2023_obiol-experiment_assemblies-and-data/annotation/peptidases_kos.tsv') %>% 
    pull(KEGG) %>% 
    str_c(., collapse = '|')
    
    
annotation_row <- eggnog.df %>% 
    mutate( function_interest = case_when(
        str_detect( Description, 'ferredoxin') ~ 'ferredoxin',
        str_detect( Description, 'Flavodoxin') ~ 'flavodoxin',
        str_detect( KEGG_ko, peptidases) ~ 'peptidase',
        str_detect( Preferred_name, 'PSM[B,D]') ~ 'proteasome',
        str_detect( Description, 'photosystem') ~ 'photosystem',
        str_detect( Description, 'ribosome') ~ 'ribosome',
        TRUE ~ NA)) %>% 
    select(gene_name, function_interest) %>% 
    column_to_rownames(var = 'gene_name')


genesofinterest <- annotation_row %>% 
    filter(!is.na(function_interest)) %>% 
    rownames()

heatmap_tara(tpm_table = matrix,
             genesel = genesofinterest,
             gene_annot = annotation_row)

heatmap_tara(tpm_table = matrix,
             genesel = rownames(matrix),
             gene_annot = annotation_row)




# What is the first cluster with the biggest abundance --------------------

heatmap_tara(tpm_table = matrix,
             genesel = allgenes$tree_row$labels[allgenes$tree_row$order][1:130],
             gene_annot = annotation_row)



vannelid_alex <- read_tsv('data/quantification/mapping/TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara.tsv')

longi <- vannelid_alex %>% 
    pivot_longer( names_to = 'genename', values_to = 'tpm', cols = -sample)


top1000 <- longi %>% 
    group_by(genename) %>% 
    summarize( meantpm = mean(tpm)) %>% 
    top_n(n = 1000, wt = meantpm) %>% 
    pull(genename)

vannmatrix <- vannelid_alex %>%
    column_to_rownames(var = 'sample') %>%
    as.matrix() %>%
    t()
    

heatmap_tara(tpm_table = log10(vannmatrix + 0.01),
             genesel = top1000)

    # general_expression <- tpm_all %>%
    #     filter(group %in% c('2012_carradec_tara', '2021_tara_polar')) %>% 
    #     # filter(TPM < 1000) %>% 
    #     mutate( TPM = log10(TPM + 0.001),
    #             gene_name = factor(gene_name,
    #                                levels = hclustord_rows),
    #             sample = factor( sample,
    #                              levels = hclustord_cols)
    #             ) %>% 
    #     arrange(gene_name, sample) %>% 
    #     ggplot( aes(x =  sample, y = gene_name)) + 
    #     geom_tile(aes(fill = TPM)) + 
    #     scale_fill_viridis_c(option = 'inferno') +
    #     theme_void() + 
    #     ggtitle( 'Expression patterns of Florenciella over Tara',
    #              subtitle = 'Genes are in Y axis, Samples are in X axis')
    # 
    # 
    # annotatedplot <- tibble( gene_name = hclustord_rows) %>% 
    #     left_join( eggnog.df, by = 'gene_name') %>% 
    #     mutate( gene_name = factor(gene_name, levels = hclustord_rows)) %>% 
    #     ggplot( aes(x = 1)) + 
    #     geom_tile(aes( x = 1, y = gene_name, fill = is.na(seed_ortholog))) + 
    #     scale_fill_manual(  values = c('white', 'black'), 
    #                         name = 'Annotated\nKegg') +
    #     theme_void()
    # 
    # 
    # allset <- envdata %>% 
    #     filter( name_downloaded %in% hclustord_cols) %>% 
    #     mutate(name_downloaded = factor(name_downloaded, levels = hclustord_cols)) %>% 
    #     arrange(name_downloaded)
    #     
    # 
    # fractionplot <- allset %>% 
    #     ggplot() + 
    #     geom_tile(aes(x = name_downloaded, y = 1, fill = size_fraction)) + 
    #     theme_void() 
    # 
    # regionplot <- allset %>% 
    #     ggplot() + 
    #     geom_tile(aes(x = name_downloaded, y = 1, fill = ocean_region)) + 
    #     theme_void() 
    # 
    # depthplot <- allset %>% 
    #     ggplot() + 
    #     geom_tile(aes(x = name_downloaded, y = 1, fill = depth)) + 
    #     theme_void() 
    # 
    # 
    # 
    # general_expression /
    #     fractionplot /
    #     regionplot /
    #     depthplot +
    #     plot_layout(heights = c(0.9, 0.033, 0.033, 0.033),
    #                 guides = 'collect') 
    # 
