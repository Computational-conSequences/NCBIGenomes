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
datasets rehydrate --directory NCBI
datasets rehydrate --directory NCBI

for EXT in cds faa fna gbff gff gtf
do
    mkdir -p ${EXT}-${DATA}
done

for GCF in $(\ls NCBI/ncbi_dataset/data)
do
    echo "working with $GCF"
    mv NCBI/ncbi_dataset/data/$GCF/$GCF*_genomic.fna fna-${DATA}/$GCF:r.fna
    mv NCBI/ncbi_dataset/data/$GCF/cds_from_genomic.fna cds-${DATA}/$GCF:r.cds
    mv NCBI/ncbi_dataset/data/$GCF/protein.faa faa-${DATA}/$GCF:r.faa
    for EXT in gbff gff gtf
    do
        mv NCBI/ncbi_dataset/data/$GCF/genomic.$EXT ${EXT}-${DATA}/$GCF:r.$EXT
    done
done
rm -r NCBI

echo "compressing files"
for EXT in cds faa fna gbff gff gtf
do
    echo "compressing files in ${EXT}-${DATA}"
    gzip -f --best ${EXT}-${DATA}/*.$EXT
done
