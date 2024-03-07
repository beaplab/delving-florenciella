library(tidyverse)

tables <- read_tsv(list.files('data/quantification/alignment/',
                              pattern = 'Numreads_EP00618_Florenciella_parvula_quant-over',
                              full.names = T), id = 'group')

table_long <- tables  %>%
    pivot_longer(
        cols = -c(sample, group),
        names_to = "gene_name",
        values_to = "reads"
    ) 

genes_present <- table_long %>% 
    # by experiment to be more restrictive in the presence
    group_by(group, gene_name) %>% 
    filter((sum(reads >= 10) >= 10) & (sum(reads > 0) <= n()* 0.9 )) %>% 
    pull(gene_name) %>% 
    unique()

table_filt <- table_long %>% 
    filter(gene_name %in% genes_present) %>% 
    mutate(group = basename(group) %>%
               str_remove('Numreads_EP00618_Florenciella_parvula_quant-over_') %>% 
               str_remove('.tsv'))


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
        
    

