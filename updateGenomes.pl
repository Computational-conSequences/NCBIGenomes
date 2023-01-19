#!/usr/bin/env perl

use strict;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;

my @acceptable = qw(
                       prokaryotes
                       animals
                       fungi
                       plants
                       protists
               );
my $acceptable = join("|",@acceptable);

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
    . qq(  This program downloads genomes from NCBI's RefSeq database\n\n)
    . qq(usage:\n)
    . qq(    $ownname -g [group] [options]\n\n)
    . qq(options:\n)
    . qq(   -g group to download [$acceptable],\n)
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
        . "        [$acceptable]\n"
        . $helpMsg;
}
if( $group !~ m{^($acceptable)$}i ) {
    die
        "    $group is not an acceptable option:\n"
        . "        [$acceptable]\n"
        . $helpMsg;
}
$group = lc($group);

my $groupMatch = ucfirst($group);
my $localDir   = "ncbi";
my $localLists = $localDir . "/genomeInfo";
my $localGnms  = $localDir . "/$groupMatch";
my $listFile
    = $group eq "prokaryotes" ? "$localLists/$group.txt"
    : "$localLists/eukaryotes.txt";
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
    print "will only enlist md5 download commands for $group\n";
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
print "reading full genome list:\n  $listFile\n";
my( $refCount,$reforigStatus ) = readGlist("$listFile","$statusMatch");
print "The whole genome files contain:\n";
for my $status ( @status ) {
    if( exists $refCount->{"$status"} ) {
        print "   ",join(" ",$refCount->{"$status"},$status,"genomes"),"\n";
    }
}

print "finding corresponding refSeq genomes:\n";
my ($heading,$refInfo,$refStatus)
    = readRefSeq($assemblyfile,$refCount,$reforigStatus);

########################################################################
######### make directories for results
########################################################################
unless( -d "$localDir" ){
    system "mkdir -p $localGnms" unless( -d "$localGnms");
}
unless( -d "$logDir" ){
    mkdir("$logDir") unless( -d "$logDir");
}

########################################################################
######## verify and bring genomes:
########################################################################
print "downloading from RefSeq:\n";
for my $status ( @status ) {
    my @ids = sort grep { $refStatus->{"$_"} eq "$status" } keys %{ $refInfo };
    my $count = @ids;
    if( $count > 0 ) {
        print "   ",join(" ","found",$count,$status,"genomes at RefSeq"),"\n";
        system("mkdir -p $localGnms/$status") unless( -d "$localGnms/$status" );
        bringGenomes($status,\@ids,$refInfo);
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

sub check_md5s {
    my $localPath = $_[0];
    if( -d "$localPath" ) {
        my $md5file = $localPath . "/md5checksums.txt";
        if( -f "$md5file" ) {
            my %md5sum = ();
            open( my $MD5F,"<","$md5file" );
            while(<$MD5F>) {
                my($md5sum,$file) = split;
                $file =~ s{.*/}{};
                $md5sum{"$file"} = $md5sum;
                #print join("\t",$file,$md5sum,$md5sum{"$file"}),"\n";
            }
            close($MD5F);
            ### learn file names
            opendir( my $LOCALD,"$localPath");
            my @files2check = grep { m{\.gz$} } readdir($LOCALD);
            closedir($LOCALD);
            my $count_f = @files2check;
            if( $count_f > 0 ) {
                my $failed = 0;
                for my $file2check ( @files2check ) {
                    my $full_file = $localPath . "/" . $file2check;
                    open( my $FL2CH,"<","$full_file" );
                    binmode($FL2CH);
                    my $md5sum = Digest::MD5->new->addfile($FL2CH)->hexdigest;
                    close($FL2CH);
                    if( $md5sum ne $md5sum{"$file2check"} ) {
                        $failed++;
                        #unlink("$full_file");
                    }
                }
                if( $failed > 0 ) {
                    return();
                }
                else {
                    return("Good2go");
                }
            }
            else {
                return();
            }
        }
        else {
            return();
        }
    }
    else {
        return();
    }
}

sub readGlist {
    my($listFile,$statusMatch) = @_;
    ### indexes for each necessary item
    my $iGroup    = '';
    my $iStatus   = '';
    my $iAssembly = '';
    ### to save the data:
    my %count     = ();
    my %status    = ();
    open( my $GNMS,"<","$listFile" )
        or die "I need a $listFile (run updateGenomeInfo.pl first)\n";
  GNMLINE:
    while(<$GNMS>) {
        chomp;
        my @items = split(/\t/,$_);
        if(  m{^#} ) {
            for my $index ( 0 .. $#items ) {
                if( $items[$index] =~ m{^Group} ) {
                    $iGroup = $index;
                }
                if( $items[$index] =~ m{^Status} ) {
                    $iStatus = $index;
                }
                if( $items[$index] =~ m{^Assembly} ) {
                    $iAssembly = $index;
                }
            }
        }
        else {
            if( $group ne "prokaryotes" ) {
                unless( $items[$iGroup] eq "$groupMatch" ) {
                    next GNMLINE ;
                }
            }
            my $status
                = $items[$iStatus] =~ m{$statusMatch}i ? ucfirst(lc($&))
                : 'none';
            if( $status eq 'none' ) {
                next GNMLINE;
            }
            my $assemblyID = $items[$iAssembly];
            $count{"$status"}++;
            $status{"$assemblyID"} = $status;
        }
    }
    close($GNMS);
    my $cstatus = keys %status;
    if( $cstatus > 0 ) {
        return(\%count,\%status);
    }
    else {
        return();
    }
}

sub readRefSeq {
    my($assemblyFile,$refCount,$refStatus) = @_;
    ### open assembly report to learn path to sequences/genome files
    my %genomeInfo = ();
    my %status     = ();
    my $headInfo   = "";
    open( my $ASSEM,"<","$assemblyFile" );
  ASSEMBLY:
    while(<$ASSEM>) {
        if( m{^#} ) {
            $headInfo .= $_;
            next ASSEMBLY;
        }
        chomp;
        my @items = split(/\t/,$_);
        ### item 10 is version status, we only want "latest"
        next ASSEMBLY if( $items[10] ne "latest");
        my $rsyncPath = $items[19] =~ m{(https|ftp)://} ? $items[19] : "none";
        next ASSEMBLY if( $rsyncPath eq "none");
        my $assembly_accession = $items[0];
        my $assembly_id = $items[17];
        ##### by checking the status I'm also checking that this is one of the
        ##### the genomes I want to download
        my $status
            = $refStatus->{"$assembly_id"} =~ m{$statusMatch} ? $&
            : "none";
        next ASSEMBLY if( $status eq "none" );
        ##### we want to use rsync, rather than ftp or wget
        $rsyncPath =~ s{https|ftp}{rsync}g;
        my $local_subdir = $assembly_accession;
        $genomeInfo{"$assembly_accession"} = $_;
        $status{"$assembly_accession"}     = $status;
    }
    close($ASSEM);
    my $clines = keys %genomeInfo;
    if( $clines > 0 ) {
        return($headInfo,\%genomeInfo,\%status);
    }
}

sub bringGenomes {
    my ($status,$refIDs,$refInfo) = @_;
    my $maxTries = 5;
    my $rsyncMD5
        = qq(rsync -aqL)
        . qq( --timeout=15 --contimeout=10 )
        . qq( );
    my $rsyncCmd
        = qq(rsync -aqL)
        . qq( --timeout=15 --contimeout=10)
        . qq( --exclude='*/')
        . qq( --delete --delete-excluded)
        . qq( --prune-empty-dirs)
        . qq( );
  ASSSEMBLYID:
    for my $gnmID ( @{ $refIDs } ) {
        my $info = $refInfo->{"$gnmID"};
        my @items = split(/\t/,$info);
        my $rsyncPath = $items[19];
        my $assembly_accession = $items[0];
        my $assembly_id = $items[17];
        ##### we want to use rsync, rather than ftp or wget
        $rsyncPath =~ s{https|ftp}{rsync}g;
        my $local_subdir = $assembly_accession;
        if( length("$local_subdir") > 1 ) {
            my $localPath = join("/",$localGnms,$status,$local_subdir);
            if( $new eq 'T' && -d "$localPath" ) {
                next ASSEMBLYID;
            }
            ##### added to make sure newed md5 file helps discover if the
            ##### remote files are the same as the local files
            my $returnStatus = 1;
            my $tries        = 1;
            my $md5command
                = "$rsyncMD5 $rsyncPath/md5checksums.txt $localPath/ 1>/dev/null";
            if( $dry eq 'T' ) {
                print $md5command,"\n";
            }
            else {
                while( $returnStatus != 0 && $tries < $maxTries ) {
                    print "   bringing md5checksums $local_subdir $status (try $tries)\n";
                    if( $tries > 1 ) {
                        sleep 60;
                    }
                    $returnStatus = qx($md5command);
                    $tries++;
                }
                ##### now let's check with new md5 file:
                if( check_md5s("$localPath") ) {
                    ## if it works, it means all files are fine
                    ## thus, do nothing here
                }
                else {
                    my $returnStatus2 = 1;
                    my $tries2        = 1;
                    while( $returnStatus2 != 0 && $tries2 < $maxTries ) {
                        print "      bringing genome files $local_subdir $status (try $tries2)\n";
                        if( $tries2 > 1 ) {
                            sleep 60;
                        }
                        $returnStatus2 =
                            qx($rsyncCmd $rsyncPath/ $localPath 1>/dev/null);
                        $tries2++;
                    }
                }
            }
        }
        else {
            print "problem with $assembly_id\n";
        }
    }
}
