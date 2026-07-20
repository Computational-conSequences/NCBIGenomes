#!/usr/bin/env perl

use strict;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use File::Temp qw( tempfile tempdir );

my $localDir    = "ncbi";
my $localLists  = $localDir . "/genomeInfo";
my $taxidfile = "$localLists/taxonomy.info.bz2";

my @prokas = qw(
                   Archaea
                   Bacteria
           );

my $matchproka = join("|",@prokas);

# Eukaryota
my @eukarya = qw(
                    fungi
                    invertebrate
                    plant
                    protozoa
                    vertebrate_mammalian
                    vertebrate_other
            );
my $matcheukarya = join("|",@eukarya);

my @allstatus = qw(
                      Complete
                      Chromosome
                      Scaffold
                      Contig
              );

my $defDry    = 'T';
my $defNew    = 'T';
my $allstatus = join("|",@allstatus);
my @status    = ();
my $group     = '';
my $dry       = $defDry;
my $new       = $defNew;

my $ownname = $0;
$ownname =~ s{\s+/}{};

my $helpMsg
    = qq(about:\n)
    . qq(  This program downloads large genome collections from NCBI's RefSeq\n)
    . qq(  database\n\n)
    . qq(usage:\n)
    . qq(    $ownname -g [group] [options]\n\n)
    . qq(options:\n)
    . qq(   -g group to download [$matchproka|Eukaryota],\n)
    . qq(      required\n)
    . qq(   -s status to download [$allstatus], can be\n)
    . qq(      more than one, default: @allstatus\n)
    . qq(   -d run a dry run indicating no longer available genomes [T|F],\n)
    . qq(      default: $defDry\n)
    . qq(   -n only bring new genomes, don't update present ones [T|F],\n)
    . qq(      default: $defNew\n)
    . qq(\n)
    ;

my $options = GetOptions(
    "g=s"    => \$group,
    "s=s{,}" => \@status,
    "d=s"    => \$dry,
    "n=s"    => \$new,
) or die "$helpMsg";

if( !$group ) {
    die
        "    You should indicate a group to download:\n"
        . "        [$matchproka|Eukaryota]\n"
        . $helpMsg;
}
if( $group !~ m{^($matchproka|Eukaryota)$}i ) {
    die
        "    $group is not an acceptable option:\n"
        . "        [$matchproka|Eukaryota]\n"
        . $helpMsg;
}
$group = lc($group);

my $matchGroup = ucfirst($group);
my $localGnms  = $localDir . "/$matchGroup";
my $assemblyfile
    = "$localLists/assembly_summary_refseq.txt";
my $logDir     = "ncbi/logs";

if( !@status || grep { m{^all$} } @status ) {
    @status = @allstatus;
}
else {
    my $countpref = @status;
    my @newstatus = ();
    if( $countpref > 0 ) {
        for my $try ( @status ) {
            if( my @matches = grep { m{\b$try\b}i } @allstatus ) {
                print $matches[0], " matched\n";
                push(@newstatus,$matches[0]);
            }
        }
    }
    my $countadded = @newstatus;
    if( $countadded == 0 ) {
        if( $countpref > 0 ) {
            die "  @status are not genome status:\n"
                . "      [$allstatus]\n";
        }
    }
    @status = @newstatus;
}
my $statusMatch = join("|",@status);
print $statusMatch,"<--status to match\n";

$dry = $dry =~ m{^(T|F)$}i ? uc($1): $defDry;
$new = $new =~ m{^(T|F)$}i ? uc($1): $defNew;

if( $dry eq "T" ) {
    print "will print dry commands for downloading $group\n";
}
else {
    print "will download $group\n";
}
if( $new eq 'T' ) {
    print "will download new files, won't update already present\n";
}
else {
    print "will run a complete update\n";
}

########################################################################
############# reading the list of genomes to ensure we know which ones
############# are in the group we want
########################################################################
#print "learning taxIDs:\n  $taxidfile\n";
#my $refTaxID = readTaxID("$taxidfile","$matchGroup");

print "finding corresponding RefSeq genomes:\n";
#my ($heading,$refInfo,$refStatus,$refCount)
#    = readRefSeq($assemblyfile,$refTaxID);
my ($heading,$refInfo,$refStatus)
    = readRefSeq($assemblyfile);
my $total2get = 0;
print "the RefSeq genome database contains:\n";
for my $status ( @status ) {
    if( exists $refStatus->{"$status"} ) {
        my $count = @{ $refStatus->{"$status"} };
        print "   there's ",join(" ",
                         $count,$status,"genomes"),"\n";
        $total2get += $count;
    }
}
if( $total2get < 1 ) {
    die "no genomes to download\n";
}
else {
    print "will download $total2get genomes from RefSeq\n";
}
########################################################################
######### make directories for results
########################################################################
my $tempFolder = tempdir("ncbi/tmp.XXXXXXXXXXXX");
unless( -d "$localDir" ){
    system "mkdir -p $localGnms" unless( -d "$localGnms");
}
unless( -d "$logDir" ){
    mkdir("$logDir") unless( -d "$logDir");
}

########################################################################
######## download genomes:
########################################################################
print "downloading from RefSeq:\n";
for my $status ( @status ) {
    if( exists $refStatus->{"$status"} ) {
        my $count = @{ $refStatus->{"$status"} };
        print "  ",join(" ","downloading",$count,$status,
                        "genomes from RefSeq"),"\n";
        my $resultsdir = "$localGnms/$status";
        system qq(mkdir -p $resultsdir) unless( -d "$resultsdir" );
        bringGenomes($status,\@{ $refStatus->{"$status"} },$resultsdir);
    }
}

#### now let's check if all directories correspond to bichos in use
#### erase otherwise
print "   checking for genomes to erase\n";
my $toerase = 0;
my $tokeep  = 0;
my $erasefl = "$logDir/eraser-$group.log";
if( -f "$erasefl" ) {
    unlink("$erasefl");
}
open( my $BORRADOR,">","$erasefl.tmp" );
for my $status ( @status ) {
    open( my $STATUSF,"|-","bzip2 -9 > $localGnms/$status.info.bz2" );
    print {$STATUSF} $heading;
    my $statusDir = join("/",$localGnms,$status);
    opendir( my $CHECKD,"$statusDir" );
    my @subdirs = grep{ m{^[A-Z]} } readdir($CHECKD);
    closedir($CHECKD);
    for my $subdir ( @subdirs ) {
        if( exists $refInfo->{"$subdir"} ) {
            print {$STATUSF} $refInfo->{"$subdir"},"\n";
            $tokeep++;
        }
        else {
            print {$BORRADOR} "rm -rf $statusDir/$subdir\n";
            $toerase++;
        }
    }
    close($STATUSF);
}
close($BORRADOR);
if( $toerase > 0 ) {
    print "$toerase directories to erase\n";
    print "$tokeep directories to keep\n";
    rename("$erasefl.tmp","$erasefl");
}
else{
    print "nothing to erase\n";
    print "$tokeep directories to keep\n";
    unlink("$erasefl.tmp");
}
print "\n\tdone with $0\n\n";

sub readTaxID {
    my ( $taxfile,$group ) = @_;
    my %txinfo  = ();
    my $currentid = '';
    print "   reading $taxfile\n";
    open( my $TID,"-|","bzip2 -qdc $taxfile" );
    while(<$TID>) {
        if( m{^Main\s+TaxID\s+(\d+)} ) {
            $currentid = $1;
        }
        elsif( m{^(domain|superkingdom)\s+} ) {
            my($label,$checkgrp,$nn) = split;
            if( $checkgrp =~ m{^$group$} ) {
                $txinfo{"$currentid"}++;
            }
        }
    }
    close($TID);
    my $count = keys %txinfo;
    print "found $count tax IDs for $group\n";
    return(\%txinfo);
}

sub readRefSeq {
    my($assemblyFile,$reftaxa) = @_;
    ### open assembly report to learn path to sequences/genome files
    my %fullInfo    = ();
    my %assemblies  = ();
    my $headInfo    = "";
    my $taxfield    = 0;
    my $statusfield = 0;
    my $groupfield  = 0;
    open( my $ASSEM,"<","$assemblyFile" );
  ASSEMBLY:
    while(<$ASSEM>) {
        if( m{^#} ) {
            $headInfo .= $_;
            if( m{assembly_accession} ) {
                s{^#+\s*}{};
                chomp;
                my @headings = split(/\t/,$_);
                for my $index ( 0 .. $#headings ) {
                    if( $headings[$index] eq "taxid" ) {
                        $taxfield = $index;
                        print "   TaxID field is $taxfield\n";
                    }
                    elsif( $headings[$index] eq "assembly_level" ) {
                        $statusfield = $index;
                        print "   Status field is $statusfield\n";
                    }
                    elsif( $headings[$index] eq "group" ) {
                        $groupfield = $index;
                        print "   Group field is $groupfield\n";
                    }
                }
            }
        }
        else {
            chomp;
            my @items = split(/\t/,$_);
            ### we only want those matching the taxIDs of genomes in group
            #if( ! exists $reftaxa->{"$items[$taxfield]"} ) {
            #    next ASSEMBLY;
            #}
            ###
            if( $matchGroup eq "Eukaryota"  ) {
                next ASSEMBLY if( $items[$groupfield] !~ m{$matcheukarya} );
            }
            else {
                my $fixgroup = $groupfield + 1;
                my $testgroup = $items[$groupfield];
                if( $testgroup eq "haploid"
                    && $items[$fixgroup] =~ m{$matchGroup}i ) {
                    $testgroup = $items[$fixgroup];
                }
                next ASSEMBLY if( $testgroup !~ m{$matchGroup}i );
            }
            #print $items[$groupfield],"\n";
            my $assembly = $items[0];
            my $status
                = $items[$statusfield] =~ m{$statusMatch} ? $&
                : "none";
            next ASSEMBLY if( $status eq "none");
            $fullInfo{"$assembly"} = $_;
            push( @{ $assemblies{"$status"} }, $assembly );
        }
    }
    close($ASSEM);
    my $clines = keys %fullInfo;
    if( $clines > 0 ) {
        return($headInfo,\%fullInfo,\%assemblies);
    }
}

sub bringGenomes {
    my ( $status,$refIDs,$resultsdir ) = @_;
    my $maxTries = 5;
    my $list = $tempFolder . "/" . "$status.list";
    open( my $LS,">","$list" );
    for my $gnmID ( @{ $refIDs } ) {
        #$gnmID =~ s{\.\d+}{};
        print {$LS} $gnmID,"\n";
    }
    close($LS);
    my $zipfile    = "$tempFolder/$status.zip";
    my $tmpncbi    = "$tempFolder/$status";
    my $gotthemdir = "$tmpncbi/ncbi_dataset/data";
    my $metadata   = "$tmpncbi/$status.json.gz";
    ########### datasets commands:
    system qq(mkdir -p $tmpncbi);
    my $getMeta
        = qq(datasets summary genome accession --inputfile $list)
        . qq( --as-json-lines | )
        . qq( gzip --best > $metadata);
    #print "$getMeta\n";
    my $downloadCMD
        = qq(datasets download genome accession --inputfile $list)
        . qq( --include all --dehydrated --no-progressbar --filename $zipfile);
    #print $downloadCMD,"\n";
    my $unzipper = qq(unzip $zipfile -d $tmpncbi);
    #print $unzipper,"\n";
    my $rehydrater
        = qq(datasets rehydrate --gzip --no-progressbar --directory $tmpncbi);
    #print $rehydrater,"\n";
    if( $dry eq 'T' ) {
        print "running dry commands:\n";
        for my $cmd ( $getMeta, $downloadCMD, $unzipper, $rehydrater ) {
            print $cmd,"\n";
            if( $cmd =~ m{summary} ) {
                print qq(mv $metadata $resultsdir/ 2>&1 > /dev/null\n);
            }
        }
    }
    else {
        print "downloading:\n";
        for my $cmd ( $getMeta, $downloadCMD, $unzipper, $rehydrater ) {
            print $cmd,"\n";
            my $try    = 1;
            my $errors = 1;
            while( $errors > 0 && $try <= $maxTries ) {
                print "  download try: $try\n";
                my $log = qx($cmd 2>&1);
                $errors = 0;
                $errors += ( $log =~ s{error}{error}ig );
                print "errors: $errors\n";
                $try++;
                sleep 5;
            }
            if( $cmd =~ m{summary} ) {
                system qq(mv $metadata $resultsdir/ 2>&1 > /dev/null);
            }
            if( $cmd =~ m{dehydrated} ) {
                if( ! -f $zipfile ) {
                    die "no $zipfile to unzip and dehydrate\n";
                }
            }
        }
    }
}
