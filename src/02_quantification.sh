nextflow run beaplab/transcriptome_metaT_quantification  \
    --fastq_sheet /scratch/datasets_symbolic_links/dataset_sheets/metatranscriptomes_datasets.csv \
    --transcriptome "data/genomic_data/transcriptomes/nucleotide_version/*.fna.gz" \
    --outdir data/quantification \
    -r main \
    -resume


# some other genomes 
nextflow run beaplab/transcriptome_metaT_quantification  \
    --fastq_sheet /scratch/datasets_symbolic_links/dataset_sheets/metatranscriptomes_datasets.csv \
    --transcriptome "/home/aauladell/datasets/genomes_assemblies_and_genes/2022_eukprot_v3/some_downloaded_by-hand_datsets/*.fna.gz" \
    --outdir data/quantification \
    -r main -resume
