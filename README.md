# NCBIGenomes
These two Gabo-friendly scripts help retrieve genomes from NCBI. I wrote them for myself, but might improve them as time permits (and improvements will be appreciated).

First the user should run "updateGenomeInfo.pl" which will retrieve a few files from NCBI that contain the information necessary to figure out genome collections. The files will be saved under the directory:
ncbi/genomeInfo

Now, those files are used by the program: "updateGenomes.pl" which will retrieve genomes from NCBI's refseq by reading the necesssary info in those files. This program has options.

1. It requires the user to type, as first argument, either of [prokaryotes|plants|fungi|animals]
2. It can accept, after the first argument, a list of genome status [Complete|Chromosome|Scaffold|Contig]. In the absence of this argument, the program willl retreieve all of them.

So, for example, if the user wanted all the Complete prokaryotic genomes, this command would do (after running updateGenomeInfo.pl with no arguments):

updateGenomes.pl prokaryotes Complete

There's no dry option, there's no much in terms of warnings, etc.
