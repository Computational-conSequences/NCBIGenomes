#!/usr/bin/perl

use strict;
my $print = $ARGV[0] eq "print" ? 1 : 0;

## remote
my $ncbiDir       = "rsync://rsync.ncbi.nlm.nih.gov";
my $genomeDir     = $ncbiDir   . "/genomes";
my $refseqDir     = $genomeDir . "/refseq";
my $taxonomyDir   = $ncbiDir   . "/pub/taxonomy";
## here
my $localDir      = "ncbi";
my $localTaxonomy = "ncbi/taxonomy";
my $localLists    = "ncbi/genomeInfo";

#### retrieve reports
if( $print > 0 ) {
    print "rsync -avz  $genomeDir/GENOME_REPORTS/ $localLists\n";
    print "rsync -avzL $refseqDir/assembly_summary_refseq.txt $localLists/\n";
    print "rsync -avzL $refseqDir/README.txt $localLists/README-refseq.txt\n";
    #print "rsync -avzL $taxonomyDir/'taxdump.tar.gz*' $localTaxonomy/\n";
}
else {
    unless( -d "$localDir" ){
        system "mkdir -p $localDir"      unless( -d "$localDir");
        system "mkdir -p $localLists"    unless( -d "$localLists");
        #system "mkdir -p $localTaxonomy" unless( -d "$localTaxonomy");
    }
    system qq(rsync -avz  --exclude="prok_*" --exclude="CLADES/" --exclude="MARKERS/" --delete-excluded $genomeDir/GENOME_REPORTS/ $localLists);
    system "rsync -avzL $refseqDir/assembly_summary_refseq.txt $localLists/";
    system "rsync -avzL $refseqDir/README.txt $localLists/README-refseq.txt";
    #system "rsync -avzL $taxonomyDir/'taxdump.tar.gz*' $localTaxonomy/";
}
