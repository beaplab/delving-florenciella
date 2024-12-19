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

eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

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
 

# Establish denom ---------------------------------------------------------

# stats for filtering afterwards
statistics <- numreads_long %>% 
    left_join(TPM_long, by = c('group', 'sample', 'gene_name')) %>% 
    group_by(gene_name) %>% 
    mutate( occ = sum(reads > 0), # occurrence
            cv = (sd(TPM, na.rm = T) / mean(TPM, na.rm = T)) * 100, #coeff var
            fcfrommean = abs( log2( TPM / mean(TPM) )) ) #difference in fold change

n.transcribed.at.least.once <- numreads_long$gene_name %>% unique() %>% length()

# Breadth of transcriptome in each sample
breadth.transcriptome <- TPM_long %>% 
    filter(TPM > 0) %>% 
    group_by(sample) %>% 
    # amount of sequences divided per total
    summarize(amount.diff.sequences = n() / ( n.transcribed.at.least.once))


samples.w.phaeo <- breadth.transcriptome  %>% 
    # at least 10% of the transcriptome present
    filter(amount.diff.sequences >= 0.1) %>% 
    pull(sample)


# the genes that are expressed in all the samples with Phaeocystis
genes.present.samples <-  numreads_long %>% 
    mutate( pres.phaeo = ifelse(sample %in% samples.w.phaeo, 'yes', 'no')) %>% 
    group_by( pres.phaeo, gene_name) %>% 
    mutate( occ = sum(reads > 0)) %>% 
    filter( pres.phaeo == 'yes' & occ == length(samples.w.phaeo)) %>% 
    pull(gene_name) %>% unique()

# calculate which genes present:
# Coefficient of variation below 200% 
# A fold change below mean of 2
denoms <- statistics %>% 
    filter(gene_name %in% genes.present.samples) %>% 
    filter(sample %in% samples.w.phaeo) %>%
    group_by(gene_name) %>% 
    filter(  cv <= 300, 
             mean(fcfrommean) <= 2)  %>% 
    pull(gene_name) %>% 
    unique()

# calculate the correlation between all the genes to obtain a subset that presents a 
# coherent similar somehow pattern
tpm.corr.mat <- TPM_long %>% 
    filter(gene_name %in% denoms) %>% 
    select(gene_name, sample, TPM) %>% 
    pivot_wider( names_from = 'gene_name', values_from = 'TPM', id_cols = sample) %>% 
    select(-sample) %>% 
    cor() 

# Viz
heats <- tpm.corr.mat %>% 
    pheatmap::pheatmap()

# Viz with kmeans of 3 on the basis of the last viz, bc I saw one big component 
# and its the one I want to catch
heats.w.kmeans <- tpm.corr.mat %>% 
    pheatmap::pheatmap(kmeans_k = 3 )

# the kmeans with a lot of iterations to have stable components 
kmeans <- kmeans(tpm.corr.mat, centers = 3, iter.max = 9999)

# check the n seqs for each cluster
table <- kmeans$cluster %>% table()
# mine was the bigger one
cluster.choosen <- table[table == max(table)] %>% names()

# a subset of the denoms I can trust
denoms.denoised <- names(kmeans$cluster[kmeans$cluster == cluster.choosen])

# viz for double check bc we are paranoid people here
pheatmap::pheatmap(tpm.corr.mat[denoms.denoised, denoms.denoised])

# pattern over the samples 
TPM_long %>% 
    filter(gene_name %in% denoms.denoised,
           sample %in% sample(sample, size = 100)
           
           ) %>% 
    ggplot(aes(sample, TPM)) + 
    geom_line(aes(group = gene_name)) 



# Saving results ----------------------------------------------------------

housedir <- 'data/statistics/housekeeping-gene-denom'
dir.create(housedir)

tibble( gene_name = denoms.denoised) %>% 
    write_tsv( str_c(housedir, 'names_denominators.tsv'))

TPM_long %>% 
    filter(gene_name %in% denoms.denoised) %>% 
    group_by(sample) %>% 
    # adding a pseudocount to avoid the zeroes
    mutate(TPM = TPM + 0.01) %>% 
    summarize(tpm_geomean = exp(mean(log(TPM))))  %>% 
    write_tsv( str_c( housedir, 'denom_geometricmean_TPM.tsv') )

numreads_long %>% 
    filter(gene_name %in% denoms.denoised) %>% 
    group_by(sample) %>% 
    # adding a pseudocount to avoid the zeroes
    mutate(reads = reads + 0.01) %>% 
    summarize(read_geomean = exp(mean(log(reads))))  %>% 
    write_tsv( str_c( housedir, 'denom_geometricmean_reads.tsv') )


eggnog_df %>% 
    filter( gene_name %in% denoms.denoised) %>% 
    write_tsv( str_c( housedir, 'functional-annotation_denominators.tsv'))

