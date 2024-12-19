library(tidyverse)


tf = tempfile(fileext = ".xlsx")
curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-022-03939-z/MediaObjects/42003_2022_3939_MOESM4_ESM.xlsx", tf)
pelagomonas_metagenome <- readxl::read_excel(tf, sheet = 10, skip = 1) %>% 
    select(`Sample Name`, `Relative abundance (%)`)


# to relative abundances
taxa.df <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm_taxo.tsv.gz') %>% 
    filter(str_detect(closest_hits_lca, 'Pelagomonas' )) %>% 
    pull(amplicon)

pelagomonas_v9 <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm-Mumu_table.tsv.gz')
pelagomonas_v9[,6:ncol(pelagomonas_v9)] <- pelagomonas_v9[,6:ncol(pelagomonas_v9)] / colSums(pelagomonas_v9[,6:ncol(pelagomonas_v9)] )

pelagomonas_v9 <- pelagomonas_v9 %>% 
    filter(amplicon %in% taxa.df) %>% 
    select(-maxrelabund, -total, -spread, -sequence) %>% 
    pivot_longer( names_to = 'sample_id_pangaea', values_to = 'relab', cols = -amplicon) %>% 
    group_by(sample_id_pangaea) %>% 
    summarize( relab = sum(relab))

envdata.df <- read_tsv('/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_Alldata_including_polar_protists/envdata_Tara_All_downloaded_metaT_nisaba.tsv') %>% 
    select(sample_id_pangaea, name_downloaded, station_label, depth, size_fraction) %>% 
    filter(!is.na(name_downloaded)) %>% 
    mutate( station_label = str_remove(station_label, pattern = 'TARA_')) %>% 
    mutate(alternative_name = str_c(station_label, '_', depth, '_', size_fraction)) %>% 
    left_join(pelagomonas_v9, by = 'sample_id_pangaea')

pelagomonas_metaT <- read_tsv(c('data/quantification/mapping/TPM_Pelagomonas_calceolata_quant-over_2012_carradec_tara.tsv', 
                                'data/quantification/mapping/TPM_Pelagomonas_calceolata_quant-over_2021_tara_polar.tsv'),
                              num_threads = 1) %>% 
    pivot_longer( names_to = 'gene', values_to = 'count', cols = -sample) %>% 
    group_by(gene) %>% 
    filter(! sum(count == 0) == n())


n.transcribed.at.least.once <- pelagomonas_metaT$gene %>% unique() %>% length()

# Breadth of transcriptome in each sample
breadth.transcriptome <- pelagomonas_metaT %>% 
    filter(count > 0) %>% 
    group_by(sample) %>% 
    # amount of sequences divided per total
    summarize(amount.diff.sequences = n() / ( n.transcribed.at.least.once))


geomean.transcriptome <-  pelagomonas_metaT %>% 
    group_by(sample) %>% 
    summarize( geomean = exp(mean(log( count + 0.1))))

allinfo_gathered <-  geomean.transcriptome %>% 
    left_join(envdata.df, by = c('sample' = 'name_downloaded')) %>% 
    left_join(pelagomonas_metagenome, by = c('alternative_name' = 'Sample Name'))

naniar::vis_miss(allinfo_gathered)
allinfo_gathered %>% 
    select(sample, alternative_name, `Relative abundance (%)`) %>% 
    filter(is.na(`Relative abundance (%)`))

allinfo_gathered %>% 
    # filter(geomean > 1)
    ggplot(aes(geomean, `Relative abundance (%)`)) + 
    geom_point() +
    geom_smooth(method = 'lm') 

allinfo_gathered %>% 
    select(geomean, relab, `Relative abundance (%)`) %>% 
    GGally::ggpairs(columnLabels = c('Geometric mean (metaT)', '% (v9 18S rRNA gene)', '% (metaG)'))
    

