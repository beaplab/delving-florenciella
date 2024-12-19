#!/usr/bin/env Rscript

library(tidyverse)
library(argparser)
library(Biostrings)

# Create a parser
p <- arg_parser("Finds by taxonomy or by sequence a taxonomic group and plots \
                the V9/V4 presence vs transcriptome presence")

# Add command line arguments
p <- add_argument(p, "--transcriptome_quantification_matrix",
                  help="the matrix to work on",
                  type="character",
                  default = '/home/agalvez/projects/V_biogeo_alex-isolates/data/quantification/mapping/BEAP0079trimmed_clean/Numreads_BEAP0079trimmed_clean_quant-over_2012_carradec_tara.tsv')
p <- add_argument(p, "--rrna_region", help="Either the V4 or V9 of the 18S rRNA gene", default='V9')
p <- add_argument(p, "--taxonomic_sel", help="Taxonomy to compare against", default=NA)
p <- add_argument(p, "--sequence_sel", help="Exact 18S sequence", default= NA)
p <- add_argument(p, "--outdir", help="the directory output", default= 'Tara_comparison')

argv <- parse_args(p)


# argv$sequence_sel = 'AGTGTCAGTATCACGATTTTACCCTCGTGTGACTGCAGAAGGCTCATTATATCAGTTTTAAGCAATGGTGTGAAATTCGGCCGGGCCTTCGGGGCCGGTCGACGGATATCCGATCGAATTAATTGGCTAATACGTGTATTAGAAGAGCGTTGCGTGAGCGACGCAGGCGGAGGCAACTCAAGGCTGATGCTTTGTGTGAAAAAGGCTTCCGTATCGGATCGAGAATGGGATACCATCGACGTTCACCATGCGAGACTGCCCTATTAGCTTTGGAAGGGAGTATAGAGCACTCTCTTGGCTTCGATGGGTACGGGGAATCAGGGTTCGGTTCCGGAGAGGGAGCTTGAGAAACGGCTCCCACTTCTAAGGAAGGCAGCAGGCGCGTAACTTACCCAATGCAAACTCTGTGAGGTAGTGACAAGAAATAACAATATGGTATTCTATACGTTTACTGTAATTGGAATGAGTACAATCTAAAACCATTAACGAGCACCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCGAGAGCATATGAAAAAGTTGTTGCGATTAAAACACTCGTAGTTGGACTATTGGCGAGGGCCTTTGACGGTCCAGGTCGATAAGGCCGGAACTGTCGGGACTCGCCTTGTTGGCGAGATTGCACCCGTGGCATTAAGTTGTTGTGGGACCGCGCTCGTCTACTTCACTTTGAATAAATTGGAGTGTTCAAAGCAGACACTATCATAGTCATTGCACATTTTAGCATAGGATGATAGAACACGACCCCCGCCGACTTCTGTTAGCCTTGGCGGCGAGGGGTAATGGTTAAAGGGGATGGTGGAAGGTCCTTGTATGGAGCTCGTTAGAGGTGAAATTCTTAGATCGGCTCTAGACAAACTGCTGCGAAGGCGTTGTACCAACGACGTTTTCATCAATCAAGAATGAAAGTTCGGGGCGCGAAGATGATCAGATACCGTTGTAGTCCGAACCGTAAACTATGCCGTCTAAGAATGTGGGGACGTTACATATCTTATGAAGACTCCCTCCGAGCTTTGTGCGAAAGCAAAGATTTTGGGCTTCGGGGGTAGTACGGCCGCAAGACTGAAACTTAAAGGAATTGACGGAGGGGCACCACCAGGCGTGGAGTTTGCGGCTCAATTTGACTCAACACGGGGAAACTCACCAGGTCCGGACAGAGGAAGGATTGACAGATTGATAGCTCTTTCTTGATTCTTTGGGAGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGATAACGAACGAGATCTTTACCTGCTTACTAGACACGCGTATCTTTACCATACAAGCCTTGGTCGCGGCGTTTGGTCCCGCGCCCCTTCGGGGGAGCGTGGACCCCTGCGGCCTCGGTGGACTGATACGTGGATGAGCTTCTTAGAGGGACTATCTGCGTTTAGCAGATGGAAGTTGAAGGCAATAACAGGTCTGTGATGCCCTTTGATGTCCTGGGCTGCACGCGTACTACAATGACTACCGCATCAAGCAACATCACCTTGGCCGAAGGGTCTGGGTAATCTTAGAAAGGTAGTCGTACTACGGATTGCTGGCTGTAACGCCGGCATGAAGCAGGAATGCCTAGTAAACGCAAGTCAACAACTTGCATTGATTACGTCCCTGCCCCTTGTACACACCGCCCGTCGCTTCAACCGATTGGCTGTAATGGTGAGGTCTCCTGAGACTGAGGCAATATGGGGTCCGTAAGGTTCCCATGCCATCTATCAAAGTTGATCGAATCCTTGCAGTTAGAGGATGAAGAAGTCGTAACAAGGTACCT'
argv$transcriptome_quantification_matrix = 'data/quantification/alignment/Numreads_EP00618_Florenciella_parvula_quant-over_2012_carradec_tara.tsv'
argv$rrna_region = 'V9'
argv$taxonomic_sel = 'Florenciella'

# Importing data ----------------------------------------------------------

if ( argv$rrna_region == 'V9' ) {
    
    mumu.df <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm-Mumu_table.tsv.gz')
    taxa.df <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v9_swarm/TARA-Oceans_18S-V9_Swarm_taxo.tsv.gz')
    
}

if ( argv$rrna_region == 'V4' ) {
    
    mumu.df <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v4_swarm/TARA-Oceans_18S-V4_Swarm-Mumu_table.tsv.gz')
    taxa.df <- read_tsv('/scratch/datasets_symbolic_links/common_databases/tara_v4_swarm/TARA-Oceans_18S-V4_Swarm-Mumu_taxo.tsv.gz')
    
}

tara.samples.df <- read_rds('/scratch/datasets_symbolic_links/environmental_data_from_datasets/tara_general_environmental_data/context_general.rds') %>% 
    as_tibble()


# Match and find either the sequences or the taxonomy --------------------------

# Check if there is a definition for taxonomy
if ( !is.na(argv$taxonomic_sel) ) {
    
    amplicon.hash <- taxa.df %>%
        filter(str_detect(closest_hits_lca, argv$taxonomic_sel )) %>% 
        pull(amplicon)
    
}

evaluate_matches_against_mumu <- function(subject_seq = argv$sequence_sel,
                                          mumu.df = mumu.df){
    
    mumu.filt <- mumu.df %>% 
        filter(nchar(sequence) >= 80)
    
    patt <- mumu.filt$sequence %>% 
        map(~matchPattern(pattern = .x,
                          subject = subject_seq,
                          max.mismatch = 5,
                          with.indels = T))
    
    names(patt) <- mumu.filt$amplicon 
    
    
    evaluate_width <- function(val){
        
        value = val@ranges@width
        
        if(is.integer(value) && length(value) == 0L){
            return(0)
        }else{
            return(sum(value))
        }
    }
    
    
    matches.widths.list <- map(patt, ~evaluate_width(.x))
    
    matches <- matches.widths.list[matches.widths.list != 0] %>%
        as_tibble() %>% 
        pivot_longer( names_to = 'amplicon',
                      values_to = 'matchlength',
                      cols = everything()) %>% 
        arrange(-matchlength)
    
    return(matches)
    
}


# Check all the sequences presenting matches
if ( !is.na(argv$sequence_sel) ) {
    
    matches <- evaluate_matches_against_mumu(subject_seq = argv$sequence_sel,
                                             mumu.df = mumu.df)
                                             
    amplicon.hash <-  matches$amplicon
    
}


# Filter out the results --------------------------------------------------

dataset <-  mumu.df %>% 
    filter(amplicon %in% amplicon.hash) %>% 
    select(amplicon, maxrelabund, total, spread,
           one_of(tara.samples.df$sample_id_pangaea))


samples.w.counts <- dataset %>% 
    select(-maxrelabund, -total, -spread) %>% 
    pivot_longer( names_to = 'sample', values_to = 'count', cols = -amplicon) %>% 
    # we only take into account the results presenting enough reads to say something
    filter( count > 10) %>%
    group_by(sample) %>% 
    summarize(total.counts.tax = sum(count))

## Obtaining Tara counts ---------------------------------------------------

rrna_counts <- samples.w.counts %>% 
    left_join(tara.samples.df %>%  select(sample_id_pangaea, station,
                                          depth, size_fraction),
              by = c('sample' = 'sample_id_pangaea')) %>% 
    arrange(station) %>% 
    mutate( samples = str_c(station,
                            str_replace(size_fraction, pattern = '\\.', 'o'),
                            str_replace(depth, 'SRF', 'SUR'),
                            sep = '_')) 

# Compare with transcriptome ----------------------------------------------

quant.transcriptome <- read_tsv( argv$transcriptome_quantification_matrix, num_threads = 1) %>% 
    pivot_longer( names_to = 'gene', values_to = 'count', cols = -sample)

n.transcribed.at.least.once <- quant.transcriptome$gene %>% unique() %>% length()

# Breadth of transcriptome in each sample
breadth.transcriptome <- quant.transcriptome %>% 
    filter(count > 0) %>% 
    group_by(sample) %>% 
    # amount of sequences divided per total
    summarize(amount.diff.sequences = n() / ( n.transcribed.at.least.once))


# PLOT --------------------------------------------------------------------

thenas <- breadth.transcriptome %>% 
    rename( sample = 'samples') %>% 
    left_join( rrna_counts,
              by = 'samples') %>% 
    filter(is.na(sample))


rrnac <- rrna_counts %>% 
    mutate( samples = case_when(
        samples == "047_0o8-20_DCM" ~ "047_0o8-inf_DCM",
        samples == "076_0o8->_DCM" ~ "076_0o8-inf_DCM",  
        samples == "078_0o8->_DCM" ~ "078_0o8-inf_DCM",  
        samples == "100_0o8->_DCM" ~ "100_0o8-inf_DCM",  
        samples == "110_0o8->_DCM" ~ "110_0o8-inf_DCM",  
        samples == "111_0o8->_DCM" ~ "111_0o8-inf_DCM",  
        samples == "128_0o8->_DCM" ~ "128_0o8-inf_DCM",  
        samples == "131_0o8->_DCM" ~ "131_0o8-inf_DCM",  
        samples == "132_0o8->_DCM" ~ "132_0o8-inf_DCM",  
        samples == "137_0o8->_DCM" ~ "137_0o8-inf_DCM",  
        samples == "151_0o8->_DCM" ~ "151_0o8-inf_DCM",  
        TRUE ~ samples ))

thenas %>% 
    mutate(station = str_sub(samples, 1,3),
           depth = str_sub(samples, -3)) %>% 
    select(samples, station, depth) %>% 
    left_join(rrna_counts, by = c("station", "depth")) %>% 
    select(samples.x, samples.y, station, depth, size_fraction) %>% 
    filter(! is.na(samples.y)) 

breadth.rrna.plot <- breadth.transcriptome %>% 
    rename( sample = 'samples') %>% 
    left_join( rrnac,
              by = 'samples') %>% 
    ggplot(aes(total.counts.tax, amount.diff.sequences)) + 
    geom_point(aes(shape = depth, color = size_fraction)) + 
    geom_hline(yintercept = 0.10, linetype = 2 ) + 
    theme_bw() + 
    scale_x_log10() +
    scale_y_continuous(labels = scales::percent) +
    xlab(str_c('Total ', argv$rrna_region, ' counts (log10 scale)')) + 
    ylab('Breadth of transcriptome in sample') + 
    ggtitle(str_c('Relationship between transcriptome breadth and ', argv$rrna_region, ' 18S rRNA gene counts'),
            subtitle = "The dashed line indicates the breadth we can be sure it's 'present'") + 
    scale_shape_discrete( name = 'Depth') + 
    scale_color_discrete( name = 'Fraction') 
    

dir.create('results/figures/metaT_18S_comparison')
ggsave(filename = 'results/figures/metaT_18S_comparison/v9.pdf',
       plot = breadth.rrna.plot )
