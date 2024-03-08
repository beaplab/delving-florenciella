#!/bin/bash

### 
# Generate the signatures for transcriptomes, genomes 
# and other material, and compare it between them 
# to obtain the relationship
###

conda activate smash

indir='data/genomic_data/'
outdir='data/genomic_data/sourmash-signatures'
mkdir $outdir

# first the genomes
mkdir $outdir/genomes

sourmash sketch dna data/genomic_data/genomes/TARA*.fa --output-dir $outdir/genomes
sourmash sketch dna data/genomic_data/genomes/*/*.fna.gz --output-dir $outdir/genomes

# second the transcriptomes
mkdir $outdir/transcriptomes
# the proteins from the genomes
sourmash sketch protein data/genomic_data/genomes/*/*.faa.gz --output-dir $outdir/transcriptomes
sourmash sketch protein data/genomic_data/genomes/*.gmove.pep.faa --output-dir $outdir/transcriptomes
# the eukprot
sourmash sketch protein data/genomic_data/transcriptomes/*.fasta --output-dir $outdir/transcriptomes

# calculate and compare the genomes
mkdir data/statistics

sourmash compare --csv data/statistics/compar_genomes.csv --ani $outdir/genomes/*.sig -o data/statistics/compar_genomes.dist
sourmash compare --csv data/statistics/compar_genomes.csv --ani $outdir/transcriptomes/*.sig -o data/statistics/compar_trans.dist

sourmash plot data/statistics/compar_genomes.dist --output-dir results/figures
sourmash plot data/statistics/compar_trans.dist --output-dir results/figures


# based on these results we chose the genomes from the SMAGs dataset 
# we will copy them to analyze them later on 

gzip data/genomic_data/genomes/*_CDS.fna
ln -s ~/projects/delving_into_florenciella/data/genomic_data/genomes/TARA_ARC_108_MAG_00262_CDS.fna.gz data/genomic_data/transcriptomes/nucleotide_version/TARA_ARC_108_MAG_00262_CDS.fna.gz 
ln -s ~/projects/delving_into_florenciella/data/genomic_data/genomes/TARA_MED_95_MAG_00475_CDS.fna.gz  data/genomic_data/transcriptomes/nucleotide_version/TARA_MED_95_MAG_00475_CDS.fna.gz 
ln -s ~/projects/delving_into_florenciella/data/genomic_data/genomes/TARA_SOC_28_MAG_00069_CDS.fna.gz  data/genomic_data/transcriptomes/nucleotide_version/TARA_SOC_28_MAG_00069_CDS.fna.gz 


