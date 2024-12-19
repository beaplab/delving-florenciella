library(tidyverse)

tpm_filt <- read_tsv('data/quantification/alignment/tpm_all_filtered.tsv')
tpm_nonfilt <- read_tsv('data/quantification/alignment/TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv')

TPM.df <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'TPM_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group', num_threads = 1) %>% 
    mutate(group = basename(group) %>%
               str_remove('TPM_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

gene_names <- colnames(tpm_nonfilt)[-1]
n.transcribed.at.least.once <- length(gene_names)


breadth <- TPM.df %>% 
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "TPM"
    )  %>% 
    filter(TPM > 0) %>% 
    group_by(sample) %>% 
    # amount of sequences divided per total
    summarize(amount.diff.sequences = n() / ( n.transcribed.at.least.once)) %>% 
    arrange(-amount.diff.sequences)

sams.w.breadth <- breadth %>% 
    filter( amount.diff.sequences >= 0.1)  %>% 
    pull(sample)


breadth %>% 
    filter( sample %in% sams.w.breadth) %>% 
    ggplot(  aes(amount.diff.sequences)) +
    geom_histogram()


tpm_filt %>% 
    filter(gene_name %in% sample(x = intersect(gene_name,gene_names), size = 5)) %>% 
    mutate( condition = ifelse(sample %in% sams.w.breadth, 'yes', 'no')) %>% 
    ggplot( aes( condition, TPM)) +
    geom_violin( aes( fill = condition)) + 
    geom_point( alpha = 0.7) + 
    coord_cartesian( ylim = c(0, 2)) + 
    facet_wrap(~gene_name, scales = 'free') 
    scale_y_log10()
    
    

# OK let's try with something else ----------------------------------------

eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

ferredoxins <- eggnog_df %>% 
    filter( str_detect( Description, 'ferredoxin')) %>% 
    pull(gene_name)


tpm_filt %>% 
    filter(gene_name %in% ferredoxins) %>% 
    mutate( condition = ifelse(sample %in% sams.w.breadth, 'yes', 'no')) %>% 
    ggplot( aes( condition, TPM)) +
    geom_violin( aes( fill = condition)) + 
    geom_point( alpha = 0.7) + 
    coord_cartesian( ylim = c(0, 2)) + 
    facet_wrap(~gene_name, scales = 'free') 
    scale_y_log10()

