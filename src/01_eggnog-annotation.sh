

INPATH='data/genomic_data/transcriptomes'
DB='/scratch/data1/aauladell/databases/eggnog'

OUTPATH='data/genomic_data/functional_annotation'

mkdir $OUTPATH

conda activate eggnog

for genome in $(ls $INPATH/*.fasta);
    do 
    name=$(basename $genome .fasta);
    emapper.py -i $genome --data_dir $DB  -o $OUTPATH/$name;
    done

