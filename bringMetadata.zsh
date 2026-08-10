if [ $# -eq 0 ]; then
    echo "this program needs an organism to work with, like"
    echo "Escherichia coli, Streptomyces, or Actinopterygii"
    exit 1
fi
BICHO="$@"
DATA=$(echo $BICHO | perl -pe 'chomp(); s{\s+}{_}g')
ROOTDIR="NCBIMD"
ROOT="$ROOTDIR/$DATA"
echo "BICHO is $BICHO"
echo "DATA is $DATA"
echo "downloading $BICHO metadata"
mkdir -p $ROOTDIR
datasets summary genome taxon "$BICHO" --as-json-lines --annotated --reference > $ROOT.json
echo "bzipping it"
bzip2 -f --best $ROOT.json
echo "producing tsv"
bzcat $ROOT.json.bz2 | dataformat tsv genome | perl -ne 'my @test = split; if ( !$seen{$test[0]} ){print; $seen{$test[0]}++}' > $ROOT.tsv
echo "bzipping it"
bzip2 -f --best $ROOT.tsv
echo "producing list of genome gcfs"
#bzcat $ROOT.tsv.bz2 | awk -F"\t" '$110 == "Complete Genome" && $173 == "SOURCE_DATABASE_REFSEQ" {print $1}' | perl -pe 's{\.\d+}{}' | sort -u > $DATA.list
bzcat $ROOT.tsv.bz2 | awk -F"\t" 'NR > 1 {print $1}' | perl -pe 's{\.\d+}{}' | sort -u > $ROOT.list
echo "done"
