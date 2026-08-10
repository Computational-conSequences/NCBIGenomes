if [ $# -eq 0 ]; then
    echo "this program needs an organism to work with, like"
    echo "Escherichia coli, Streptomyces, or Actinopterygii"
    exit 1
fi
BICHO="$@"
DATA=$(echo $BICHO | perl -pe 'chomp(); s{\s+}{_}g')
LIST="NCBIMD/$DATA.list"
DLDIR="NCBIDL"
echo "BICHO is $BICHO"
echo "DATA is $DATA"
echo "downloading $BICHO genomes to $DLDIR"
datasets download genome accession --inputfile $LIST --include all --dehydrated --filename $DATA.zip
unzip $DATA.zip -d $DLDIR
echo "Extracting stuff"
datasets rehydrate --gzip --directory $DLDIR
datasets rehydrate --gzip --directory $DLDIR

for EXT in cds faa fna gbff gff
do
    mkdir -p ${EXT}-${DATA}
done

for GCF in $(\ls $DLDIR/ncbi_dataset/data)
do
    echo "working with $GCF"
    mv $DLDIR/ncbi_dataset/data/$GCF/$GCF*_genomic.fna.gz fna-${DATA}/$GCF:r.fna.gz
    mv $DLDIR/ncbi_dataset/data/$GCF/cds_from_genomic.fna.gz cds-${DATA}/$GCF:r.cds.gz
    mv $DLDIR/ncbi_dataset/data/$GCF/protein.faa.gz faa-${DATA}/$GCF:r.faa.gz
    for EXT in gbff gff gtf
    do
        mv $DLDIR/ncbi_dataset/data/$GCF/genomic.$EXT.gz ${EXT}-${DATA}/$GCF:r.$EXT.gz
    done
done
rm -r $DLDIR
