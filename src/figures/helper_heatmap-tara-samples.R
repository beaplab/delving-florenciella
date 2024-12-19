library(tidyverse)
library(pheatmap)
library(viridis)


# heatmap_tara(matrix,
#              genesofinterest = genesofinterest,
#              annotation_row = annotation_row )


heatmap_tara <- function(tpm_table, genesel, gene_annot = NULL, silent = F){
    
    envdata <- read_tsv('/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_Alldata_including_polar_protists/envdata_Tara_All_downloaded_metaT_nisaba.tsv')
    
    # Annotation samples ------------------------------------------------------
    size_fraction_order <- c("0.8->", "0.8-5",  "3-20",  "5-20",  "20-180", "180-2000")
    ocean_region_order <- c('North Atlantic Ocean',  'South Atlantic Ocean',
                            'Mediterranean Sea', "Red Sea",
                            'Indian Ocean', "Southern Ocean", 
                            "North Pacific Ocean", "South Pacific Ocean",
                            "Arctic Ocean")
    depth_order <- c('SRF', 'DCM', 'MES', 'FSW')
    
    
    annotation_col <- envdata %>% 
        select(name_downloaded,depth, size_fraction, ocean_region) %>% 
        filter(!is.na(name_downloaded)) %>%
        distinct(name_downloaded, .keep_all = T) %>% 
        mutate( size_fraction = factor(size_fraction, levels = size_fraction_order),
                ocean_region = str_remove(ocean_region,
                                          pattern = '\\[.*\\] ') %>% 
                    str_remove(pattern = ' \\(.*\\)') %>% 
                    factor(., levels = ocean_region_order), 
                depth = ifelse(depth %in% c('MIX', 'ZZZ'), NA, depth) %>% 
                    factor(., levels = depth_order)
        ) %>% column_to_rownames(var = 'name_downloaded')
    
    
    depth_colors <- c("#00BFFF", "#1874CD", "#092E52", 'black')
    names(depth_colors) <- depth_order
    
    ocean_region_color <- c('#94C5CCFF', '#278B9AFF', '#5E81ACFF',
                            '#E75B64FF', '#D8AF39FF',
                            '#E8C4A2FF', "orange", '#E67D65', 'skyblue')
    names(ocean_region_color) <- ocean_region_order
    
    size_fraction_color <- c('#B4DAE5FF', 'dodgerblue', '#F7EABDFF',
                             '#F0D77BFF', '#403369FF', '#1D2645FF')
    names(size_fraction_color) <- size_fraction_order
    
    ann_colors <- list( depth = depth_colors,
                        ocean_region = ocean_region_color,
                        size_fraction = size_fraction_color)
    
    pheatmap(tpm_table[rownames(tpm_table) %in% genesel,],
             border_color = NA,
             drop_levels = F,
             color = inferno(n = 100),
             annotation_colors = ann_colors,
             annotation_col = annotation_col,
             annotation_row = gene_annot,
             show_colnames = F, silent = silent)
    
    
}



