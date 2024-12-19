library(tidyverse)
library(Biostrings)

# we need to connect the eggnog annotation with the nucleotide name
fasta_aa <- readAAStringSet('data/genomic_data/transcriptomes/EP00618_Florenciella_parvula.fasta')
eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula.fasta.emapper.annotations',
                      comment = "##") %>% 
    dplyr::rename( prot_name = `#query`)  

parse_correspondence <- function(string){
    relationship <- string %>% str_split_1(pattern = ' ') %>% .[c(1,10)]
    list( prot_name = relationship[1],
            gene_name = str_remove(relationship[2], pattern = '/DNA_ID=')) %>% 
        return()
    
}

correspondence_df <- names(fasta_aa) %>%
    map_df(~ parse_correspondence(.x))  

eggnog_joined_df <- eggnog_df %>% 
    left_join( correspondence_df, by = 'prot_name') %>% 
    select(gene_name, everything(), -prot_name)

dir <- 'data/genomic_data/functional_annotation/'
path <- str_c(dir, 'EP00618_Florenciella_parvula_nucnames.tsv')
write_tsv(x = eggnog_joined_df, 
          file = path)

