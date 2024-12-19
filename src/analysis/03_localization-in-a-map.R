library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")


# Data import -------------------------------------------------------------

tpm_filt <- read_tsv('data/quantification/alignment/tpm_all_filtered.tsv')

taraenv <- read_tsv(str_c('/scratch/datasets_symbolic_links/environmental_data_from_datasets/',
                          'tara_Alldata_including_polar_protists/envdata_Tara_All_downloaded_metaT_nisaba.tsv'))

aleixenv <- read_tsv( str_c('/scratch/datasets_symbolic_links/environmental_data_from_datasets',
                            '/obiol_bacterivory_experiments_metadata'))

world <- ne_countries(scale = "medium", returnclass = "sf")


# Plotting patterns -------------------------------------------------------
nseqs <- tpm_filt$gene_name %>% unique() %>% length()

tpm_filt_summ <- tpm_filt %>%  
    group_by(group, sample) %>% 
    summarize( breadth = sum(TPM > 0) / nseqs,
               mean_tpm = mean(TPM), 
               geomean_tpm = exp(mean(log(TPM + 0.00001)))
)  %>% 
    ungroup()


tara <- tpm_filt_summ %>% 
    filter(group %in% c('2012_carradec_tara', '2021_tara_polar')) %>% 
    left_join( taraenv, by = c('sample' = 'name_downloaded' ))
    
alldatamap <- tara %>% 
    mutate( depth = case_when( 
        depth == 'MIX' ~ 'SRF', 
        depth == 'FSW' ~ NA, 
        depth == 'ZZZ' ~ 'SRF',
        TRUE ~ depth )) %>% 
    filter(!is.na(depth)) %>% 
    mutate( depth = factor(depth, levels = c('SRF','DCM', 'MES')),
            size_fraction = factor(size_fraction, c("0.8->", "0.8-5", '3-20',
                                                    '5-20', '20-180', '180-2000'))) %>% 
    filter( size_fraction != '180-2000')


sampling_loc <- taraenv %>% 
    filter(depth %in% c('SRF','DCM', 'MES')) %>% 
    select(event_longitude, event_latitude, depth) %>% 
    distinct() %>% 
    mutate( depth = factor(depth, levels = c('SRF','DCM', 'MES')))

word_plot <- ggplot(world) + 
    geom_sf() + 
    geom_point( data = sampling_loc,
                aes( x= event_longitude, y = event_latitude), 
                shape = 21,
                size = 0.5,
                color = 'black')

world_tpm <- world_plot + 
    geom_point( data = alldatamap,
                aes( x= event_longitude,
                     y = event_latitude,
                     size = mean_tpm , 
                     color = size_fraction)) + 
    facet_grid(size_fraction~depth) + 
    ylab('') + 
    xlab('') + 
    scale_size_continuous( name = 'mean TPM') + 
    theme_minimal()

dir.create('results/figures/map_metat')

ggsave(filename = 'results/figures/map_metat/world_tpm.pdf',
       plot = world_tpm, width = 8, height = 7)


world_geotpm <- world_plot + 
    geom_point( data = alldatamap,
                aes( x= event_longitude,
                     y = event_latitude,
                     size = geomean_tpm ,
                     color = size_fraction)) +
    facet_grid(size_fraction~depth) + 
    ylab('') + 
    xlab('') + 
    scale_size_continuous( name = 'geomean TPM') + 
    theme_minimal()

ggsave(filename = 'results/figures/map_metat/world_geotpm.pdf',
       plot = world_geotpm, width = 8, height = 7)


world_breadth <- world.plot + 
    geom_point( data = alldatamap,
                aes( x= event_longitude,
                     y = event_latitude,
                     size = breadth, 
                     color = size_fraction)) + 
    facet_grid(size_fraction~depth) + 
    ylab('') + 
    xlab('') + 
    scale_size_continuous( name = 'Breadth') + 
    theme_minimal()

ggsave(filename = 'results/figures/map_metat/world_breadth.pdf',
       plot = world_breadth, width = 8, height = 7)


# What about Aleix experiment? --------------------------------------------
eggnog_df <- read_tsv('data/genomic_data/functional_annotation/EP00618_Florenciella_parvula_nucnames.tsv')

tpm_filt %>% 
    filter(group == '2023_obiol_blanes-experiment') %>% 
    left_join(aleixenv, by = c('sample' = 'RunAccession')) %>% 
    left_join(eggnog_df, by = 'gene_name') %>% 
    mutate(Description = ifelse( TPM > 1000, Description, NA)) %>% 
    ggplot( aes( Time_h, TPM)) + 
    geom_point( aes( color = Experiment)) +
    geom_text( aes( label = Description))  

# niente    


# Ferredoxin --------------------------------------------------------------
ferredoxins <- eggnog_df %>% 
    filter( str_detect( Description, 'ferredoxin')) %>% 
    pull(gene_name)

ferr <- tpm_filt %>% 
    filter(gene_name %in% ferredoxins) %>% 
    left_join( alldatamap %>% 
    select(sample, event_longitude, event_latitude, depth, size_fraction) %>% 
    unique(), by = 'sample') %>% 
    filter( !is.na(size_fraction),  !is.na(depth)) 

world_ferredoxin <- world_plot + 
    geom_point( data = ferr,
                aes( x= event_longitude,
                     y = event_latitude,
                     size = TPM, 
                     color = size_fraction)) + 
    facet_grid(size_fraction~depth) + 
    ylab('') + 
    xlab('') + 
    theme_minimal()


ggsave(filename = 'results/figures/map_metat/world_ferredoxin.pdf',
       plot = world_ferredoxin, width = 8, height = 7)


# Relationship between breadth and geomean and tpm ------------------------

rel_plot <- alldatamap %>% 
    select(sample, breadth, mean_tpm, geomean_tpm) %>% 
    pivot_longer( names_to = 'type', values_to = 'mean', cols = -c(sample, breadth)) %>% 
    ggplot( aes( breadth, mean)) + 
    geom_point() +
    facet_wrap( ~type, scales = 'free') +
    ggtitle('Geometric mean captures more efficiently our expectation regarding breadth vs average expression')

ggsave(filename = 'results/figures/statttttts/relationship_breadth_vs_means.png',
       plot = rel_plot, width = 6, height = 6)
