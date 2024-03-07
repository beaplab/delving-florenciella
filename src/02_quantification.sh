nextflow run beaplab/transcriptome_metaT_quantification  \
    --fastq_sheet /scratch/datasets_symbolic_links/dataset_sheets/metatranscriptomes_datasets.csv \
    --transcriptome "data/genomic_data/transcriptomes/nucleotide_version/EP00618_Florenciella_parvula.fna.gz" \
    --outdir data/quantification \
    -r main \
    -resume
