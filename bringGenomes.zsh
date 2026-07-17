if [ $# -eq 0 ]; then
    echo "this program needs an organism to work with, like"
    echo "Escherichia coli, Streptomyces, or Actinopterygii"
    exit 1
fi
BICHO="$@"
DATA=$(echo $BICHO | perl -pe 'chomp(); s{\s+}{_}g')
LIST=DataSets/$DATA.list
echo "BICHO is $BICHO"
echo "DATA is $DATA"
echo "downloading $BICHO genomes"
datasets download genome accession --inputfile DataSets/$DATA.list --include all --dehydrated --filename DataSets/$DATA.zip
unzip DataSets/$DATA.zip -d NCBI
echo "Extracting stuff"
datasets rehydrate --gzip --directory NCBI
datasets rehydrate --gzip --directory NCBI

for EXT in cds faa fna gbff gff gtf
do
    mkdir -p ${EXT}-${DATA}
done

for GCF in $(\ls NCBI/ncbi_dataset/data)
do
    echo "working with $GCF"
    mv NCBI/ncbi_dataset/data/$GCF/$GCF*_genomic.fna.gz fna-${DATA}/$GCF:r.fna.gz
    mv NCBI/ncbi_dataset/data/$GCF/cds_from_genomic.fna.gz cds-${DATA}/$GCF:r.cds.gz
    mv NCBI/ncbi_dataset/data/$GCF/protein.faa.gz faa-${DATA}/$GCF:r.faa.gz
    for EXT in gbff gff gtf
    do
        mv NCBI/ncbi_dataset/data/$GCF/genomic.$EXT.gz ${EXT}-${DATA}/$GCF:r.$EXT.gz
    done
done
rm -r NCBI
