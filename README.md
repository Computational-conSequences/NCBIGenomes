# NCBIGenomes
These two Gabo-friendly scripts help retrieve genomes from NCBI. I wrote them for myself, but might improve them as time permits (and improvements from you will be appreciated).

First the user should run "updateGenomeInfo.pl" which will retrieve a few files from NCBI that contain the information necessary to figure out genome collections. The files will be saved under the directory:

ncbi/genomeInfo

Now, those files are used by the program: "updateGenomes.pl" which will retrieve genomes from NCBI's refseq by reading the necesssary info in those files. This program has options, which you can learn by running the program with no arguments.

1. It requires the user to type, as first argument, either of [prokaryotes|animals|fungi|plants|protists]

2. It can accept, after the first argument, a list of genome status [Complete|Chromosome|Scaffold|Contig]. In the absence of this argument, the program willl retreieve all of them.

So, for example, if the user wanted all the Complete prokaryotic genomes, this command would do (after running updateGenomeInfo.pl with no arguments):

updateGenomes.pl -g prokaryotes -s Complete -d F -n F

If you want to try a small set of genomes:

updateGenomes.pl -g bacteria -s Chromosome -d F -n F

The genomes are saved under the ncbi directory. For example, if you run:

updateGenomes.pl -g bacteria -d F -n F

The resulting genomes will be here:

ncbi/Bacteria

# NCBI has stopped rsync downloads (there is an rsync server, but it doesn't seem to have the right directory structure). So, here a few scripts that work with NCBI's command line programs (cli, from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)

1. installNCBIcli.zsh
- this zsh script will install the CLI in a mac, under an /opt/biotools directory in a mac, you must have the directory already, and it should be yours. You may modify it to install somewhere else if you have no power.

2. bringMetadata.zsh
- this zsh script will download genome info from NCBI, using NCBI's "datasets" command, and extract a sublist of refseq genomes. For example, zsh bringMetadata.zsh Klebsiella will bring info from most Klebsiella genomes at NCBI, and extract a list of the RefSeq ones.

3. bringGenomes.zsh
- this program will bring the genomes from the list, but check it, it expects the list to be in a "DataSets" directory. If you ran the example command above, this one, with the same argument, will use the resulting list and download the corresponding Klebsiella genomes.

Use with loads and loads of caution (run at your own risk and all of that, no warranties whatsoever).
