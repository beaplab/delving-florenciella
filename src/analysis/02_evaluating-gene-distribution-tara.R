library(tidyverse)

numreads_filt <- read_tsv('data/quantification/alignment/numreads_all_filtered.tsv')
denoms_geomean <- read_tsv('data/statistics/housekeeping-gene-denomdenom_geometricmean_reads.tsv')

eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

ratio <- numreads_filt %>% 
    left_join(denoms_geomean, by = 'sample') %>% 
    mutate(ratio = reads / read_geomean) %>% 
    select(-read_geomean) %>% 
    filter(group %in% c('2012_carradec_tara', '2021_tara_polar'))


ratio  %>% 
    group_by(gene_name) %>% 
    summarize( mean = mean(ratio),
               variance = var(ratio)) %>% 
    View()
