#!/usr/bin/env perl

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
my $rsync = qq(rsync -avzL);
    my $longoptions
        = qq(--exclude="prok_*" )
        . qq(--exclude="CLADES/" )
        . qq(--exclude="MARKERS/" )
        . qq(--delete-excluded )
        . qq($genomeDir/GENOME_REPORTS/);

#### retrieve reports
if( $print > 0 ) {
    print "$rsync $longoptions $localLists\n";
    print "$rsync $refseqDir/assembly_summary_refseq.txt $localLists/\n";
    print "$rsync $refseqDir/README.txt $localLists/README-refseq.txt\n";
    #print "$rsync $taxonomyDir/'taxdump.tar.gz*' $localTaxonomy/\n";
}
else {
    unless( -d "$localDir" ){
        system "mkdir -p $localDir"      unless( -d "$localDir");
        system "mkdir -p $localLists"    unless( -d "$localLists");
        #system "mkdir -p $localTaxonomy" unless( -d "$localTaxonomy");
    }
    print "running:\n$rsync $longoptions $localLists\n";
    my $genomeReports
        = qx($rsync $longoptions $localLists);
    print "running:\n$rsync $refseqDir/assembly_summary_refseq.txt $localLists/\n";
    my $assemblysumm
        = qx($rsync $refseqDir/assembly_summary_refseq.txt $localLists/);
    print "running:\n$rsync $refseqDir/README.txt $localLists/README-refseq.txt\n";
    my $readme
        = qx($rsync $refseqDir/README.txt $localLists/README-refseq.txt);
    #print "$rsync $taxonomyDir/'taxdump.tar.gz*' $localTaxonomy/\n";
    #my $tax = qx($rsync $taxonomyDir/'taxdump.tar.gz*' $localTaxonomy/);
}
