## Delving into Florenciella

A small project to evaluate how many species of Florenciella are abundant in TARA. 

It is going to be used to generate a general pipeline to process our transcriptomes at BEAP lab. 


### Problem definition 

*Florenciella parvula* presents the rank 31 in the *Tara* Oceans dataset for the V9 region of the 18S rRNA gene. For the V4 region however the rank is way lower, the 288. A possible explanation of this difference is that there are multiple species or ecotypes collapsed inside the V9 abundance distribution, splitted in the V4 case. 

### Objective 

Obtain the V9, V4, transcriptomes, housekeeping gene trees and try to disentangle how many species are abundant and what is its prevalence over the diferent regions. 

### Approaches 

1) Compare how are the V9 and V4 datasets distributed. 
2) Obtain the sMAGs from [Delmont et al. 2020](https://www.genoscope.cns.fr/tara/) and evaluate for the Florenciella sMAGs how the reads are distributed. 
3) Evaluate for the phylogenetic trees of the housekeeping genes of interest if the placement presents discriminations between the 4 known species. 

### Relevant literature

- [Vannier et al. 2016](https://www.nature.com/articles/srep37900) analyzed two Bathycoccus genomes presenting the identical V9, and showing that these two ecotypes are differentiated on environmental axes. 
- [Guérin et al. 2022](https://www.nature.com/articles/s42003-022-03939-z) evaluated the distribution of Pelagomonas, its possible niches and its adapatations on a genomic level.
