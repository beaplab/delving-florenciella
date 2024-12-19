library(tidyverse)


crossmapping <- read_tsv('data/statistics/crossmapping_Florenciella-parvula_vs_RCC1007.tsv') %>% 
    filter(category_cross %in% c('large crossmapping (>50%)',
                                 'only reads crossmapped!',
                                 'intermediate situation (10-50%)')) %>% 
    pull(refname.x)

## we will calculate it both for the BWA and Salmon
# Data import -------------------------------------------------------------
numreads.df <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'Numreads_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group', num_threads = 1) %>% 
    mutate(group = basename(group) %>%
               str_remove('Numreads_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

TPM.df <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'TPM_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group', num_threads = 1) %>% 
    mutate(group = basename(group) %>%
               str_remove('TPM_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

numreads_long <- numreads.df  %>%
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "reads"
    ) 

TPM_long <- TPM.df  %>%
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "TPM"
    ) 
    

# Filter genes based on abundance  ----------------------------------------
genes_present <- numreads_long %>% 
    # by experiment to be more restrictive in the presence
    group_by(group, gene_name) %>% 
    filter((sum(reads >= 10) >= 10) & (sum(reads > 0) <= n()* 0.9 )) %>% 
    # removing crossmapped genes 
    filter(!gene_name %in% crossmapping) %>% 
    pull(gene_name) %>% 
    unique()

numreads_filt <- numreads_long %>% 
    filter(gene_name %in% genes_present)

TPM_filt <- TPM_long %>% 
    filter(gene_name %in% genes_present) 

write_tsv(x = numreads_filt, 
          file = 'data/quantification/alignment/numreads_all_filtered.tsv')

write_tsv(x = TPM_filt, 
          file = 'data/quantification/alignment/tpm_all_filtered.tsv')

## Salmon calculation
# Data import -------------------------------------------------------------
numreads.df <- read_tsv(list.files('data/quantification/mapping/',
                              pattern = 'Numreads_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group', num_threads = 1) %>% 
    mutate(group = basename(group) %>%
               str_remove('Numreads_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

TPM.df <- read_tsv(list.files('data/quantification/mapping/',
                              pattern = 'TPM_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group', num_threads = 1) %>% 
    mutate(group = basename(group) %>%
               str_remove('TPM_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

numreads_long <- numreads.df  %>%
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "reads"
    ) 

TPM_long <- TPM.df  %>%
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "TPM"
    ) 
    

# Filter genes based on abundance  ----------------------------------------
genes_present <- numreads_long %>% 
    # by experiment to be more restrictive in the presence
    group_by(group, gene_name) %>% 
    filter((sum(reads >= 10) >= 10) & (sum(reads > 0) <= n()* 0.9 )) %>% 
    # removing crossmapped genes 
    filter(!gene_name %in% crossmapping) %>% 
    pull(gene_name) %>% 
    unique()

numreads_filt <- numreads_long %>% 
    filter(gene_name %in% genes_present)

TPM_filt <- TPM_long %>% 
    filter(gene_name %in% genes_present) 

write_tsv(x = numreads_filt, 
          file = 'data/quantification/mapping/numreads_all_filtered.tsv')

write_tsv(x = TPM_filt, 
          file = 'data/quantification/mapping/tpm_all_filtered.tsv')
