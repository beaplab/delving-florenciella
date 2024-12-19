library(tidyverse)


# Data importing ----------------------------------------------------------
florenciella <- read_tsv('data/quantification/mapping/TPM_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv') 

corr_genes_to_eggnog <- read_tsv('~/projects/gene_explorer/data/gene_data/linkage_eggnogdb/correspondence_annotations.tsv') %>% 
    #TODO do it more elegantly
    distinct(Preferred_name, .keep_all = T)

annotation <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula.emapper.annotations', 
                       comment = "##") %>% 
    rename( gene_name = `#query`) %>% 
    left_join(corr_genes_to_eggnog, by = 'Preferred_name')
 

crossmapping <- read_tsv('data/statistics/crossmapping_florenciella_against_all_transcriptomes.tsv')


# Rank abundance ----------------------------------------------------------

themostcross <- crossmapping %>% 
    filter(comparison == 'Pelagomonas') 

quant.florenciella <- florenciella %>% 
    pivot_longer( names_to = 'gene_name', values_to = 'tpm', cols = -sample)

rank.ab.gene <- quant.florenciella %>% 
    group_by(gene_name) %>% 
    summarize( median.tpm = median(tpm), 
               variance.tpm = var(tpm)) %>% 
    arrange(-median.tpm) %>% 
    mutate( rank = 1:nrow(.))

rank.ab.gene %>% 
    filter(rank < 25) %>% 
    ggplot(aes(rank, median.tpm)) + 
    geom_bar(stat = 'identity') + 
    ylab('Median TPM') + 
    xlab('Rank') + 
    theme_bw()
        
putative_housekeeping <- annotation %>% 
    filter(!is.na(label)) %>% 
    select(gene_name, label) 

putative_housekeeping %>% 
    left_join(themostcross, by = 'gene_name') %>% 
    left_join(rank.ab.gene, by = 'gene_name') %>% 
    filter(category_cross == 'large crossmapping (>50%)') %>% 
    View()
    
    

nice.rank.df <- rank.ab.gene %>% 
    filter(rank < 400) %>% 
    mutate( which_rank = case_when(
        rank <= 25 ~ 'rank 1-25',
        rank <= 100 ~ 'rank 25 - 100', 
        rank <= 200 ~ 'rank 100 - 200', 
        rank <= 300 ~ 'rank 200 - 300',
        TRUE ~ 'rank 300 - 400'
        ) %>% 
            fct_inorder()) %>% 
    left_join(themostcross, by = 'gene_name') %>% 
    mutate( category_cross = factor(category_cross, levels = c("no crossmapping",
                                                               "small amount of crossmapping (<10%)",
                                                               "intermediate situation (10-50%)", 
                                                               "large crossmapping (>50%)", 
                                                               "only reads crossmapped!"))) %>% 
    left_join(putative_housekeeping, by = 'gene_name')

plot_rank <- nice.rank.df %>% 
    ggplot(aes( y = median.tpm, x = rank)) +
    geom_bar(stat = 'identity', aes(fill = category_cross, color = label)) + 
    facet_wrap(~which_rank, ncol = 1, scales = 'free') + 
    ylab('Median TPM') + 
    xlab('Rank') +
    scale_fill_brewer() + 
    scale_color_discrete(na.value = 'transparent', name = 'Functions') + 
    theme_bw()


ggsave(filename = 'results/figures/rank-abundance-w-crossmapping-Florenciella.png', 
       plot = plot_rank, width = 14, height = 8)

obtain_coverage_df <- function(path){
    
    lines.df <- read_lines(path)
    
    fasta <- list()
    
    for (line in lines.df) {
        if(str_detect(line, pattern = '>')){
            
            fasta[[line]] <- ""
            lastname <- line
            
        }
        else{
            
            fasta[[lastname]] <- str_c(fasta[[lastname]], line, sep = ' ')
            
        }
        
    }
    
    
    tibble( gene_name = names(fasta) %>% str_remove(., '>'), 
            coverage = unlist(fasta)) %>% 
        mutate( coverage = str_trim(coverage) %>% str_split(pattern = ' '),
                position = map(coverage, ~1:length(.x))) %>% 
        return()
    
}

cov.df1 <- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/082_0o8-5_SUR_cov.txt.gz')
cov.df4 <- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/083_0o8-5_SUR_cov.txt.gz')
cov.df8 <- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/085_0o8-5_SUR_cov.txt.gz')

cov.dfsomenumber <- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/ERR6349778_cov.txt.gz')
cov.another <- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/046_0o8-5_SUR_cov.txt.gz')
cov.038_085_dcm<- obtain_coverage_df('data/quantification/alignment/EP00618_Florenciella_parvula/038_0o8-5_DCM_cov.txt.gz')


ribosome_genes <- annotation %>% 
                filter(label == 'Ribosomal protein') %>% 
                pull(gene_name)

thetopribosomes <- rank.ab.gene %>% 
    filter(rank < 100) %>% 
    filter( gene_name %in% ribosome_genes) %>% 
    pull(gene_name)

other_genes <- rank.ab.gene %>% 
    filter(rank < 100) %>% 
    filter( !gene_name %in% ribosome_genes) %>% 
    pull(gene_name)


plot_coverage <- function(df){
    
    df %>%  
        ggplot( aes(x = position, y = as.double(coverage))) + 
        geom_bar(stat = 'identity') +
        facet_wrap(~gene_name, scales = 'free') + 
        ylab('Coverage') + 
        xlab('Gene position')
    
}

listallcovs <- list(sample1 = cov.df1, 
                         sample4 = cov.df4,
                         sample8 = cov.df8,
                         sample_idk = cov.dfsomenumber,
                         sample_idk2 = cov.another,
                         sample_038_085_dcm = cov.038_085_dcm
                    
                    )

ribosomes_cov_df <- listallcovs %>% 
    map_df(~.x %>% 
               filter(gene_name %in% thetopribosomes) %>% 
               unnest(cols =c(coverage, position)), .id = 'samples') 

others_cov_df <- listallcovs %>% 
    map_df(~.x %>% 
               filter(gene_name %in% other_genes) %>% 
               unnest(cols =c(coverage, position)), .id = 'samples') 

selection_random_ribosomes <- sample(unique(ribosomes_cov_df$gene_name), 6)

ribosomes_cov_df %>% 
    filter(gene_name %in% selection_random_ribosomes) %>% 
    filter(samples == 'sample1') %>%
    plot_coverage()  + 
    ggtitle('Coverage of ribosome genes in one sample')

ribosomes_cov_df %>%
    filter(gene_name  == 'CAMNT_0007505551') %>% 
    plot_coverage() + 
    facet_wrap(~samples, scales = 'free_y') + 
    ggtitle('Ribosomal gene #43 over the 4 samples')


selection_others <- sample(other_genes, 6)

others_cov_df %>% 
    filter(gene_name %in% selection_others) %>% 
    filter(samples == 'sample1') %>% 
    plot_coverage() 

others_cov_df %>% 
    # filter(samples == 'sample1') %>% 
    plot_coverage() 

others_cov_df %>%
    filter(gene_name  == 'CAMNT_0007509633') %>% 
    plot_coverage() +
    facet_wrap(~samples, scales = 'free_y') + 
    ggtitle('Gene ranked #24, annotated as an Ubiquitin')

others_cov_df %>%
    filter(gene_name  == 'CAMNT_0007507209') %>% 
    plot_coverage() +
    facet_wrap(~samples, scales = 'free_y') + 
    ggtitle('#36, a HSP90, annotated as HSP83-1')


others_cov_df %>% 
    filter(samples == 'sample_038_085_dcm') %>% 
    pull(coverage) %>% as.double() %>%  sum(.)/200


genes_favored_bwa <-  str_c('CAMNT_000', c( 7589389, 7517805, 7577019,
                                            7584511, 7511655 ))

genes_favored_salmon <-  str_c('CAMNT_000', c(7589635, 7567723,
                                           7587191, 7589635, 7589177))

listallcovs %>% 
    map_df(~.x %>% 
               filter(gene_name %in% genes_favored_bwa) %>% 
               unnest(cols =c(coverage, position)), .id = 'samples') %>% 
    plot_coverage() +
    facet_grid(gene_name~samples, scales = 'free_y') 

listallcovs %>% 
    map_df(~.x %>% 
               filter(gene_name %in% genes_favored_salmon) %>% 
               unnest(cols =c(coverage, position)), .id = 'samples') %>% 
    plot_coverage() +
    facet_wrap(gene_name~samples, scales = 'free') 

