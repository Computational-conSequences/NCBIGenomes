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
######### make directories for results
########################################################################
my $tempFolder = tempdir("ncbi/tmp.XXXXXXXXXXXX");
unless( -d "$localGnms" ){
    system qq(mkdir -p $localGnms);
}
unless( -d "$logDir" ){
    system qq(mkdir -p $logDir);
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
my ($heading,$refInfo,$refStatus) = readRefSeq($assemblyfile);
my $total2get = 0;
print "the RefSeq genome database contains:\n";
for my $status ( @status ) {
    if( exists $refStatus->{"$status"} ) {
        my $count = @{ $refStatus->{"$status"} };
        print join(" ",
                   $status,"genomes",$count),"\n";
        $total2get += $count;
        open( my $STATUSF,"|-","bzip2 -9 > $localGnms/$status.info.bz2" );
        print {$STATUSF} $heading;
        for my $gcf ( @{ $refStatus->{"$status"} } ) {
            if( exists $refInfo->{"$gcf"} ) {
                print {$STATUSF} $refInfo->{"$gcf"},"\n";
            }
        }
        close($STATUSF);
    }
}
if( $total2get < 1 ) {
    die "no genomes to download\n";
}
else {
    print "will download $total2get genomes from RefSeq\n";
}

########################################################################
######## check for already present genomes/genomes to erase:
########################################################################
#### now let's check if all directories correspond to bichos in use
#### erase otherwise
print "   checking for genomes to erase\n";
my $toerase = 0;
my %keepers = 0;
my $erasefl = "$logDir/eraser-$group.log";
if( -f "$erasefl" ) {
    unlink("$erasefl");
}
open( my $BORRADOR,">","$erasefl.tmp" );
CHECKINGGNMS:
for my $status ( @status ) {
    my $statusDir = join("/",$localGnms,$status);
    if( ! -d $statusDir ) {
        next CHECKINGGNMS;
    }
    opendir( my $CHECKD,"$statusDir" );
    my @subdirs
        = grep{ m{^GC\S+_\d+} and ( -d "$statusDir/$_" ) } readdir($CHECKD);
    closedir($CHECKD);
    for my $subdir ( @subdirs ) {
        if( exists $refInfo->{"$subdir"} ) {
            $keepers{"$subdir"}++;
        }
        else {
            print {$BORRADOR} "rm -rf $statusDir/$subdir\n";
            $toerase++;
        }
    }
}
close($BORRADOR);
my $tokeep = keys %keepers;

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
        bringGenomes($status,\@{ $refStatus->{"$status"} },
                     $resultsdir,\%keepers);
    }
}

print "\n\tdone with $0\n\n";

########################################################################
########################################################################
########################################################################
###################### subroutines: ####################################
########################################################################
########################################################################
########################################################################
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
    my $assemblyFile = $_[0];
    #my($assemblyFile,$reftaxa) = @_;
    ### open assembly report to learn path to sequences/genome files
    my %fullInfo    = ();
    my %assemblies  = ();
    my $headInfo    = "";
    my $taxfield    = 0;
    my $statusfield = 0;
    my $groupfield  = 0;
    my $excluded    = 0;
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
                    elsif( $headings[$index] eq "excluded_from_refseq" ) {
                        $excluded = $index;
                        print "   excluded field is $excluded\n";
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
                    s{Biological Resource Center,\s*\t}{Biological Resource Center, };
                    @items = split(/\t/,$_);
                    print $items[$groupfield],"<---fixed?\n";
                }
                next ASSEMBLY if( $testgroup !~ m{$matchGroup}i );
            }
            next ASSEMBLY if( $items[$excluded] ne "na" );
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
    my ( $status,$refIDs,$resultsdir,$keepers ) = @_;
    my $toget = 0;
    my $list  = $tempFolder . "/" . "$status.list";
    my $subls = $tempFolder . "/" . "$status.sublist";
    open( my $LS,">","$list" );
    open( my $SUB,">","$subls" );
    for my $gnmID ( @{ $refIDs } ) {
        #$gnmID =~ s{\.\d+}{};
        print {$LS} $gnmID,"\n";
        #print "new is $new ($gnmID),",$keepers->{"$gnmID"},"\n";
        if( $new eq 'T' ) {
            if( ! exists $keepers->{"$gnmID"} ) {
                print {$SUB} $gnmID,"\n";
                $toget++;
            }
        }
        else {
            print {$SUB} $gnmID,"\n";
            $toget++;
        }
    }
    close($LS);
    close($SUB);
    if( $toget == 0 ) {
        print "   nothing to download\n";
    }
    else {
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
        my $downloadCMD
            = qq(datasets download genome accession --inputfile $subls)
            . qq( --include all --dehydrated --no-progressbar)
            . qq( --filename $zipfile);
        my $unzipper
            = qq(unzip $zipfile -d $tmpncbi);
        my $rehydrater
            = qq(datasets rehydrate --gzip --no-progressbar --max-workers 20)
            . qq( --directory $tmpncbi);
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
            my $maxTries = 10;
            for my $cmd ( $getMeta, $downloadCMD, $unzipper, $rehydrater ) {
                my $try    = 1;
                my $errors = 1;
                while( $errors > 0 && $try <= $maxTries ) {
                    print $cmd,": try $try\n";
                    my $log = qx($cmd 2>&1);
                    $errors = 0;
                    $errors += ( $log =~ s{error}{error}ig );
                    print "errors in try $try: $errors\n";
                    $try++;
                    sleep 60;
                }
                if( $cmd =~ m{summary} ) {
                    my $metatsv = cleanMeta($metadata);
                    print  qq(mv $metadata $resultsdir/ 2>&1 > /dev/null\n);
                    system qq(mv $metadata $resultsdir/ 2>&1 > /dev/null);
                    system qq(mv $metatsv  $resultsdir/ 2>&1 > /dev/null);
                }
                if( $cmd =~ m{dehydrated} ) {
                    if( ! -f $zipfile ) {
                        die "no $zipfile to unzip and dehydrate\n";
                    }
                }
            }
            ####### now move to proper directory
            print "now moving files to $resultsdir\n";
            opendir( my $NCBI,$gotthemdir );
            my @tomove
                = grep { m{GC\S+_\d+} and -d "$gotthemdir/$_" } readdir($NCBI);
            closedir($NCBI);
            for my $tomove ( @tomove ) {
                my $moveit
                    = qq(rsync -av --delete --remove-source-files)
                    . qq( $gotthemdir/$tomove/ $resultsdir/$tomove);
                #print $moveit,"\n";
                my $transfer = qx($moveit 2>&1);
                my $errors = 0;
                $errors += ( $transfer =~ s{}{}g );
                if ( $errors == 0 ) {
                    system qq(rmdir $gotthemdir/$tomove);
                }
            }
            print "   done moving files\n";
        }
    }
}

sub cleanMeta {
    my $metajson = $_[0];
    my $metatsv  = $metajson;
    $metatsv =~ s{json}{tsv};
    my %seen = ();
    print "cleaning $metajson to $metatsv\n";
    open( my $TSV,"|-","gzip --best > $metatsv.tmp" );
    for my $line ( qx(gzip -qdc $metajson | dataformat tsv genome) ) {
        my @test = split(/\s+/,$line);
        if( ! exists $seen{"$test[0]"} ) {
            print {$TSV} $line;
            $seen{"$test[0]"}++;
        }
    }
    close($TSV);
    rename( "$metatsv.tmp","$metatsv" );
    return($metatsv);
}
