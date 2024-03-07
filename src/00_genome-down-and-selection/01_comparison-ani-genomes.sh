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




