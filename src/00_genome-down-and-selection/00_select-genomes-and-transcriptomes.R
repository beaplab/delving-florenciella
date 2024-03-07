library(tidyverse)
library(fs)


# from the metadata from TOPAZ, SMAGs, and transcriptomes (eukprot) we
# select the ones presenting taxonomic information related to Florenciella


# TOPAZ genomes -----------------------------------------------------------
taxonomy.topaz.df <- read_csv('~/scratch4/TOPAZ_MAGs/TableS02_EukaryoticMAG.csv')

# only selecting a subset of the genomes that can potentially be Florenciellas
selection.topaz <- taxonomy.topaz.df %>% 
  filter(str_detect(mmseqs_taxonomy, 'Florencie') |
           str_detect(eukulele_taxonomy, 'Florencie')) %>% 
  pull(1) %>% 
  unique()


seltopaz_path <- str_c('~/scratch4/TOPAZ_MAGs/Eukaryotic_TOPAZ_MAGs/', 
                      selection.topaz)

dir_copy(seltopaz_path,
         new_path = str_c('~/projects/delving_into_florenciella/data/genomic_data/genomes/',
                          selection.topaz), 
         overwrite = T)



# SMAGs genomes -----------------------------------------------------------
report.mags.df <- readxl::read_excel('~/scratch4/sMAGs/supp_tables/Table_S03_statistics_nr_SMAGs_METdb.xlsx',
                                     skip = 2)

selection.smags <- report.mags.df %>% 
  filter(Best_taxonomy_GENRE == 'Florenciella', 
         Database == 'TARA_SMAGs' ) %>% 
  pull(`Genome_Id final names`)
  

# we will also save part of the information
report.mags.df %>% 
  filter(`Genome_Id final names` %in% selection.smags) %>% 
  select(1, total_length, Estimated_length, Nombre_de_genes,
         ANVIO_completion, ANVIO_redundancy, BUSCO_completion, BUSCO_redundancy)

selsmags_path <- str_c('~/scratch4/sMAGs/Contigs/', 
                      selection.smags,
                      '.fa')

file_copy(selsmags_path,
         new_path = str_c('~/projects/delving_into_florenciella/data/genomic_data/genomes/',
                          selection.smags,
                          '.fa'), 
         overwrite = T)

# we will also copy the proteins
selsmags_aa_path <- str_c('~/scratch4/sMAGs/Genes/Peptides/', 
                      selection.smags,
                      '.gmove.pep.fa')

file_copy(selsmags_aa_path,
          new_path = str_c('~/projects/delving_into_florenciella/data/genomic_data/genomes/',
                           selection.smags,
                           '.gmove.pep.faa'), 
          overwrite = T)

# Eukprot transcriptomes --------------------------------------------------

sel_transcriptomes <- list.files('~/scratch4/eukprot/proteins',
                                 pattern = 'Florenciella',
                                 full.names = T)


file_copy(sel_transcriptomes,
         new_path = str_c('~/projects/delving_into_florenciella/data/genomic_data/transcriptomes/',
                          list.files('~/scratch4/eukprot/proteins',
                                 pattern = 'Florenciella')), 
         overwrite = T)

# that's it! We have the most up to date possible genomes in the data
# lets play with them now 
