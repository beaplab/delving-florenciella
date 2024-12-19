

INPATH='data/genomic_data/transcriptomes/nucleotide_version'
DB='/scratch/data1/aauladell/databases/eggnog'

OUTPATH='data/genomic_data/functional_annotation'

mkdir $OUTPATH

conda activate eggnog

for genome in $(ls $INPATH/*.fna.gz);
    do 
    name=$(basename $genome .fna.gz);
    emapper.py --itype CDS --genepred search --cpu 12 -i $genome --data_dir $DB  -o $OUTPATH/$name;
    done


INPATH='/home/aauladell/datasets/genomes_assemblies_and_genes/2022_eukprot_v3/some_downloaded_by-hand_datsets/'


for genome in $(ls $INPATH/*.fna.gz);
    do 
    name=$(basename $genome .fna.gz);
    emapper.py --itype CDS --genepred search --cpu 4 -i $genome --data_dir $DB  -o $OUTPATH/$name;
    done


