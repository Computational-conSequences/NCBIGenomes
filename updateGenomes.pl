#!/usr/bin/perl

use strict;
### group of genomes to bring
my $trygroup = shift @ARGV;

my @acceptable = qw(
                       prokaryotes
                       plants
                       fungi
                       animals
               );
my $acceptable = join("|",@acceptable);

my $group
    = $trygroup =~ m{^($acceptable)$}i ? lc($1) : "none";
if( $group eq "none" ) {
    die "    the first argument should be the group to download:\n"
        . "        [$acceptable]";
}
my $groupMatch = ucfirst($group);
my $localDir   = "ncbi";
my $localLists = $localDir . "/genomeInfo";
my $localGnms  = $localDir . "/$groupMatch";
my $listFile
    = $group eq "prokaryotes" ? "$localLists/$group.txt"
    : "$localLists/eukaryotes.txt";
my $logDir     = "ncbi/logs";
unless( -d "$localDir" ){
    system "mkdir -p $localGnms" unless( -d "$localGnms");
}
unless( -d "$logDir" ){
    mkdir("$logDir") unless( -d "$logDir");
}

### directories for each status
my @allStatus = qw(
                   Complete
                   Chromosome
                   Scaffold
                   Contig
              );

my @status = ();
my @preferred = @ARGV;
my $countpref = @preferred;
if( $countpref > 0 ) {
    for my $try ( @preferred ) {
        if( my @matches = grep { m{\b$try\b}i } @allStatus ) {
            print $matches[0], " matched\n";
            push(@status,$matches[0]);
        }
    }
}

my $countadded = @status;
if( $countadded == 0 ) {
    if( $countpref > 0 ) {
        die "  your arguments are no genome status:\n"
            . "       [" . join("|",@allStatus) . "]";
    }
    else {
        push(@status,@allStatus);
    }
}

my $matcher = join("|",@status);
print $matcher,"<--what we will match\n";

for my $status ( @status ) {
    system("mkdir -p $localGnms/$status")
        unless( -d "$localGnms/$status" );
}

my $maxTries = 5;
$maxTries++;

my $rsyncMD5 =  'rsync -aqzL'
    . ' --timeout=15 --contimeout=10 '
    . ' ';

my $rsyncCmd =  'rsync -aqzL'
    . ' --timeout=15 --contimeout=10'
    . ' --include="*_assembly_stats.txt"'
    . ' --include="*_assembly_regions.txt"'
    . ' --include="*_genomic.fna.gz"'
    . ' --include="*_genomic.gbff.gz"'
    . ' --include="*_genomic.gff.gz"'
    . ' --include="*_protein.faa.gz"'
    . ' --include="*_rna.fna.gz"'
    . ' --include="md5checksums.txt"'
    . ' --exclude="*/"'
    . ' --exclude="*"'
    . ' --delete --delete-excluded'
    . ' --prune-empty-dirs'
    . ' ';

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
open( my $GNMS,"<","$listFile" );
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

for my $status ( sort keys %count ) {
    print join("\t",$status,$count{"$status"}),"\n";
}

### open assembly report to learn path to sequences/genome files
my %genome_info = ();
my $head_info   = "";
open( my $ASSEM,"<","$localLists/assembly_summary_refseq.txt" );
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
    my $rsync_path = $items[19] =~ m{ftp://} ? $items[19] : "none";
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
    $rsync_path =~ s{ftp}{rsync}g;
    my $local_subdir = $assembly_accession;
    $genome_info{"$assembly_accession"} = $_;
    if( length("$local_subdir") > 1 ) {
        my $local_path = join("/",$localGnms,$status,$local_subdir);
        ##### added to make sure newed md5 file helps discover if the
        ##### remote files are the same as the local files
        my $returnStatus = 1;
        my $tries        = 1;
        while( $returnStatus != 0 && $tries < $maxTries ) {
            print "   bringing md5checksums $local_subdir (try $tries)\n";
            if( $tries > 1 ) {
                sleep 60;
            }
            $returnStatus =
                system "$rsyncMD5 $rsync_path/md5checksums.txt $local_path/ 1>/dev/null";
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
                print "      bringing genome files $local_subdir (try $tries2)\n";
                if( $tries2 > 1 ) {
                    sleep 60;
                }
                $returnStatus2 =
                    system "$rsyncCmd $rsync_path/ $local_path 1>/dev/null";
                $tries2++;
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
open( my $BORRADOR,">","$logDir/eraser-$group.log" );
for my $status ( @status ) {
    open( my $STATUSF,"|-","bzip2 -9 > $localGnms/$status.info.bz2" );
    print {$STATUSF} $head_info;
    my $status_dir = join("/",$localGnms,$status);
    opendir( my $CHECKD,"$status_dir" );
    my @subdirs = grep{ m{^[A-Z]} } readdir($CHECKD);
    closedir($CHECKD);
    for my $subdir ( @subdirs ) {
        if( length($genome_info{"$subdir"}) > 0 ) {
            print {$STATUSF} $genome_info{"$subdir"},"\n";
        }
        else {
            print {$BORRADOR} "rm -rf $status_dir/$subdir\n";
        }
    }
    close($STATUSF);
}
close($BORRADOR);
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
                    my $checkline = qx(md5sum $full_file);
                    my ($md5sum,$checkedfile) = split(/\s+/,$checkline);
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
