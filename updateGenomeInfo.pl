#!/usr/bin/perl

use strict;
use Bio::DB::EUtilities;
use XML::Simple;
use File::Temp qw( tempfile tempdir );
#use Data::Dumper;

my $print = $ARGV[0] eq "print" ? 1 : 0;

#### for taxonomy stuff
my $setLn   = 10000;
my $takeNap = 5;
my @cleanTaxonomy = qw(
                          group
                          complex
                          cluster
                          subcluster
                          family
                          symbionts
                          subdivisions
                  );
my $cleanTaxonomy = join("|",@cleanTaxonomy);

## remote
my $ncbiDir       = "https://ftp.ncbi.nlm.nih.gov";
my $genomeDir     = $ncbiDir   . "/genomes";
my $refseqDir     = $genomeDir . "/refseq";
my $taxonomyDir   = $ncbiDir   . "/pub/taxonomy";
## here
my $localDir      = "ncbi";
my $localLists    = "ncbi/genomeInfo";
my $ftp = qq(wget -N);

mkdir("$localLists") unless( -d "$localLists" );

#### retrieve reports
if( $print > 0 ) {
    print "$ftp $refseqDir/assembly_summary_refseq.txt -P $localLists/\n";
    print "$ftp $refseqDir/README.txt $localLists/README-refseq.txt\n";
}
else {
    unless( -d "$localDir" ){
        system "mkdir -p $localDir"      unless( -d "$localDir");
        system "mkdir -p $localLists"    unless( -d "$localLists");
    }
    print
        "running:\n"
        . "$ftp $refseqDir/assembly_summary_refseq.txt -P $localLists/\n";
    my $assemblysumm
        = qx($ftp $refseqDir/assembly_summary_refseq.txt -P $localLists/);
    print
        "running:\n"
        . "$ftp $refseqDir/README.txt $localLists/README-refseq.txt\n";
    my $readme
        = qx($ftp $refseqDir/README.txt -O $localLists/README-refseq.txt);
}

my $refTI = findTaxIDs("$localLists/assembly_summary_refseq.txt");
bringTaxonomy($refTI);
print "done with $0\n\n";

sub bringTaxonomy {
    my $txids = $_[0];
    my $count   = 0;
    my @tobring = ();
    my $outtaxid = "$localLists/taxonomy.info";
    print "   will save tax IDs to $outtaxid.bz2\n";
    my $refdone = ( -f "$outtaxid.bz2" ) ? checkTaxID("$outtaxid.bz2") : ();
    my $countback = @{ $txids };
    open( my $OTI,">","$outtaxid.tmp" );
  TXID:
    for my $txid ( @{ $txids } ) {
        #print $txid,"\n";
        if( exists $refdone->{"$txid"} ) {
            print {$OTI} "Main TaxID\t$txid\n",$refdone->{"$txid"};
        }
        else {
            push(@tobring,$txid);
            $count++;
        }
        if( $count == $setLn || $txid == $txids->[-1] ) {
            if( $count > 0 ) {
                print "bringing $count tax IDs of $countback\n";
                my $hashedTI = getTaxonomy(@tobring);
                my $done = @tobring;
                $countback -= $done;
                $count   = 0;
                @tobring = ();
                for my $taxid ( sort { $a <=> $b } keys %{ $hashedTI } ) {
                    print {$OTI} "Main TaxID\t$taxid\n",$hashedTI->{"$taxid"};
                }
                sleep($takeNap);
            }
        }
    }
    close($OTI);
    rename("$outtaxid.tmp","$outtaxid");
    system("bzip2 -f --best $outtaxid");
}

sub checkTaxID {
    my $taxfile = $_[0];
    my %txinfo  = ();
    my $currentid = '';
    print "   reading $taxfile\n";
    open( my $TID,"-|","bzip2 -qdc $taxfile" );
    while(<$TID>) {
        if( m{^Main\s+TaxID\s+(\d+)} ) {
            $currentid = $1;
        }
        else {
            $txinfo{"$currentid"} .= $_;
        }
    }
    close($TID);
    my $count = keys %txinfo;
    print "found $count already downloaded tax IDs\n";
    return(\%txinfo);
}

sub findTaxIDs {
    my $assemblyfile = $_[0];
    open( my $ASS,"<","$assemblyfile" ) or die $!;
    my %taxid = ();
    my $taxfield = 0;
    while(<$ASS>) {
        if( m{^#} ) {
            if( m{assembly_accession} ) {
                s{^#+\s*}{};
                chomp;
                my @headings = split(/\t/,$_);
              FIELD:
                for my $try ( @headings ) {
                    if( $try eq "taxid" ) {
                        last FIELD;
                    }
                    else {
                        $taxfield++;
                    }
                }
                print "The taxid field is $taxfield\n";
            }
        }
        else {
            chomp;
            my @items = split(/\t/,$_);
            #print $items[$taxfield],"<-TAXID\n";
            $taxid{"$items[$taxfield]"}++;
        }
    }
    close($ASS);
    my @taxids = sort { $a <=> $b } keys %taxid;
    my $ctaxids = @taxids;
    print "  found $ctaxids taxids\n";
    return(\@taxids);
}

sub getTaxonomy {
    my @taxIDs = @_;
    my $factory = Bio::DB::EUtilities->new(
        -eutil => 'efetch',
        -email => 'gmorenohagelsieb@wlu.ca',
        -db    => 'taxonomy',
        -id    => \@taxIDs,
    );
    my $fullres  = $factory->get_Response->content;
    my $fulldata = XMLin($fullres);
    #print Dumper($fulldata);
    ## Taxon is an array, thus the ->[0]->:
    ## my $scientific = $fulldata->{Taxon}->[0]->{ScientificName};
    my %fullInfo = ();
    my $taxArray = $fulldata->{Taxon};
    for my $entry ( @{ $taxArray } ) {
        my $scientific = $entry->{ScientificName};
        my $lineage_a  = $entry->{Lineage};
        $lineage_a =~ s{^cellular\s+organisms;\s+}{};
        my $xlineage   = $entry->{LineageEx};
        my $test_xid   = $entry->{TaxId};
        my %lineage    = %{$xlineage};
        my @fullInfo   = ();
        for my $key ( sort keys %lineage) {
            my @onemore = @{$lineage{"$key"}};
            for my $second ( @onemore ) {
                my %onemore = %{$second};
                my @ranks = ();
                for my $third ( sort keys %onemore ) {
                    my $info = $onemore{"$third"};
                    push(@ranks,$info);
                }
                my $allranks = join("\t",@ranks);
                unless( $allranks =~ m{cellular\s+organisms}
                        && $allranks =~ m{no\s+rank} ) {
                    push(@fullInfo,$allranks);
                }
            }
        }
        my @shortened = ();
        for my $taxinf ( split(/;\s+/,$lineage_a) ) {
            if( $taxinf =~ m{cellular\s+organisms}i ||
                $taxinf =~ m{\b($cleanTaxonomy)\b}i ) {
                #print "eliminating $taxinf\n";
            }
            else {
                push( @shortened,$taxinf );
            }
        }
        my $lineage = join("; ",@shortened);
        $lineage =~ s{(\.|,|;|:)$}{};
        my $jointInfo
            = join("\t","ScientificName",$scientific) . "\n"
            . join("\t","FullTaxInfo",$lineage_a) . "\n"
            . join("\t","BasicTaxInfo",$lineage)  . "\n"
            . join("\t","Rank","ScientificName","TaxId") . "\n"
            . join("\n",@fullInfo)."\n";
        $fullInfo{"$test_xid"} = $jointInfo;
        ##return("$lineage_a","$lineage","$scientific","$test_xid",\@fullInfo);
    }
    return(\%fullInfo);
}
