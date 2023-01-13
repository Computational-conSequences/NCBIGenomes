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
$dry = $dry =~ m{^(T|F)$}i ? uc($1): $defDry;
$new = $new =~ m{^(T|F)$}i ? uc($1): $defNew;

### indexes for each necessary item
#my ($iGroup,$iAssembly,$iStatus)
#    = $group eq "prokaryotes" ? ("NA",18,15) : (4,8,16);
### open reports to learn genome status
### to hold the indexes for the different items we care about:
my $iGroup    = '';
my $iStatus   = '';
my $iAssembly = '';
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
        my $status     = $items[$iStatus];
        my $assemblyID = $items[$iAssembly];
        $count{"$status"}++;
        $status{"$assemblyID"} = $status;
    }
}
close($GNMS);

if( $dry eq "T" ) {
    print "will only enlist md5 download commands for $group\n";
}
else {
    print "will download $group\n";
}
my $matcher = join("|",@status);
print $matcher,"<--status to match\n";
if( $new eq 'T' ) {
    print "will download new files, won't update already present\n";
}
else {
    print "will run a complete update\n";
}
unless( -d "$localDir" ){
    system "mkdir -p $localGnms" unless( -d "$localGnms");
}
unless( -d "$logDir" ){
    mkdir("$logDir") unless( -d "$logDir");
}

for my $status ( @status ) {
    system("mkdir -p $localGnms/$status")
        unless( -d "$localGnms/$status" );
}

my $maxTries = 5;
$maxTries++;

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

for my $status ( sort keys %count ) {
    print join("\t",$status,$count{"$status"}),"\n";
}

### open assembly report to learn path to sequences/genome files
my %genome_info = ();
my $head_info   = "";
open( my $ASSEM,"<","$assemblyfile" );
ASSEMBLY:
while(<$ASSEM>) {
    if( m{^#} ) {
        $head_info .= $_;
        next ASSEMBLY;
    }
    chomp;
    my @items = split(/\t/,$_);
    ### item 10 is version status, we only want "latest"
    next ASSEMBLY if( $items[10] ne "latest");
    my $rsync_path = $items[19] =~ m{(https|ftp)://} ? $items[19] : "none";
    next ASSEMBLY if( $rsync_path eq "none");
    my $assembly_accession = $items[0];
    my $assembly_id = $items[17];
    ##### by checking the status I'm also checking that this is one of the
    ##### the genomes I want to download
    my $status
        = $status{"$assembly_id"} =~ m{$matcher} ? $&
        : "none";
    next ASSEMBLY if( $status eq "none" );
    ##### we want to use rsync, rather than ftp or wget
    $rsync_path =~ s{https|ftp}{rsync}g;
    my $local_subdir = $assembly_accession;
    $genome_info{"$assembly_accession"} = $_;
    if( length("$local_subdir") > 1 ) {
        my $local_path = join("/",$localGnms,$status,$local_subdir);
        if( $new eq 'T' && -d "$local_path" ) {
            next ASSEMBLY;
        }
        ##### added to make sure newed md5 file helps discover if the
        ##### remote files are the same as the local files
        my $returnStatus = 1;
        my $tries        = 1;
        my $md5command
            = "$rsyncMD5 $rsync_path/md5checksums.txt $local_path/ 1>/dev/null";
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
            if( check_md5s("$local_path") ) {
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
                        qx($rsyncCmd $rsync_path/ $local_path 1>/dev/null);
                    $tries2++;
                }
            }
        }
    }
    else {
        print "problem with $assembly_id\n";
    }
}
close($ASSEM);

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
    print {$STATUSF} $head_info;
    my $statusDir = join("/",$localGnms,$status);
    opendir( my $CHECKD,"$statusDir" );
    my @subdirs = grep{ m{^[A-Z]} } readdir($CHECKD);
    closedir($CHECKD);
    for my $subdir ( @subdirs ) {
        if( length($genome_info{"$subdir"}) > 0 ) {
            print {$STATUSF} $genome_info{"$subdir"},"\n";
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
    my $local_path = $_[0];
    if( -d "$local_path" ) {
        my $md5file = $local_path . "/md5checksums.txt";
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
            opendir( my $LOCALD,"$local_path");
            my @files2check = grep { m{\.gz$} } readdir($LOCALD);
            closedir($LOCALD);
            my $count_f = @files2check;
            if( $count_f > 0 ) {
                my $failed = 0;
                for my $file2check ( @files2check ) {
                    my $full_file = $local_path . "/" . $file2check;
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
