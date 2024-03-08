library(tidyverse)


# Data import -------------------------------------------------------------
numreads.df <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'Numreads_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group') %>% 
    mutate(group = basename(group) %>%
               str_remove('Numreads_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))

TPM.df <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'TPM_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group') %>% 
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
        values_to = "reads"
    ) 
    

# Filter genes based on abundance  ----------------------------------------
genes_present <- numreads_long %>% 
    # by experiment to be more restrictive in the presence
    group_by(group, gene_name) %>% 
    filter((sum(reads >= 10) >= 10) & (sum(reads > 0) <= n()* 0.9 )) %>% 
    pull(gene_name) %>% 
    unique()

numreads_filt <- numreads_long %>% 
    filter(gene_name %in% genes_present)

TPM_filt <- TPM_long %>% 
    filter(gene_name %in% genes_present) 


# Disentangle gene information --------------------------------------------

summary_genes <- TPM_filt %>% 
    group_by(group, gene_name) %>% 
    summarize( mean.reads = mean(reads), 
               ocurrence = sum(reads > 0))

summary_genes %>% 
    ggplot( aes( ocurrence, mean.reads)) + 
    geom_point() + 
    facet_wrap(~group, scales = 'free') + 
    scale_y_log10() 




## Evaluate functional annotation for some groups --------------------------

# we need to connect the eggnog annotation with the nucleotide name
fasta_aa <- Biostrings::readAAStringSet('data/genomic_data/transcriptomes/EP00618_Florenciella_parvula.fasta')

parse_correspondence <- function(string){
    relationship <- string %>% str_split_1(pattern = ' ') %>% .[c(1,10)]
    tibble( prot_name = relationship[1],
            nuc_name = str_remove(relationship[2], pattern = '/DNA_ID=')) %>% 
        return()
    
}

correspondence_df <- names(fasta_aa) %>%
    map_df(~ parse_correspondence(.x))  

eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula.fasta.emapper.annotations',
                      comment = "##") %>% 
    rename( prot_name = `#query`)  %>% 
    left_join( correspondence_df, by = 'prot_name') %>% 
    select(nuc_name, everything(), -prot_name)

summary_genes %>%
    group_by(group) %>%
    top_n(n = 20, wt = mean.reads)  %>%
    group_by(gene_name) %>%
    summarize(ocurrence_in_datasets = n())  %>%
    left_join(eggnog_df, by = c('gene_name' = 'nuc_name')) %>% 
    pull(PFAMs) %>% table()
    View()
    
summary_genes 
eggnog_df$PFAMs %>% table() %>% sort() %>% tail()


# Other -------------------------------------------------------------------


correlations <- table_filt  %>% 
    select(sample, gene_name, reads) %>% 
    pivot_wider( id_cols = sample, values_from = reads, names_from = gene_name) %>% 
    .[,-1] %>% 
    cor(method = 'spearman')


pc <-  cmdscale( as.dist(correlations), k = 4 ) %>% 
    tibble(MDS1 = .[,1],
           MDS2 = .[,2],
           gene = rownames(.)) %>% 
    select(gene, MDS1, MDS2)


ggplot(pc, aes(MDS1, MDS2)) + 
    geom_point() + 
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(-0.2,0.2))

%>% 
    as_tibble() %>% 
    as.dist()
dim(correlations)



# table_filt %>% 
#     group_by(group) %>% 
#     summarize( overallab = sum(reads)) %>% 
#     arrange(-overallab)
# 
# table_filt %>% 
#     group_by(gene_name) %>% 
#     summarize( ocurrence = sum(reads > 0)) %>% 
#     arrange(-ocurrence)
# 
# table_long$sample %>% unique()    %>% length()
        
    

