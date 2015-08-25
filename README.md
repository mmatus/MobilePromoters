# MobilePromoters
Perl pipeline to identify mobile promoters in prokaryotic genomes. 

Citation

## How to use it

1. Environment. The pipeline was written to run in a UNIX platform.

2. Make databases

	a. Prokaryotic genomes to be analyzed. You need the .fna and .ptt files provided by the NCBI ftp site. ftp://ftp.ncbi.nih.gov/genomes/Bacteria/
	
	b. RFAM database provided as a single file (.cm). ftp://selab.janelia.org/pub/Rfam

3. Install dependencies

	a. Perl (of course)

	b. Blastall

	c. Netclust. http://www.bioinformatics.nl/netclust/Download/nc/index.html
	
	d. CD-HIT. http://www.bioinformatics.org/project/filelist.php?group_id=350
	
	e. INFERNAL. http://infernal.janelia.org/
	
	f. IS Finder. Use it online with file *.ReprFinal. http://www-is.biotoul.fr/

4. Modify parameters.txt file located in pipeline/

To modify the settings of the pipeline, edit the input file ‘parameters.txt’. The mandatory parameters you need to establish are marked with a double asterisk (**), and they are the name of the run, the directory with the genomes, the directory with the RFAM database and the directory with the CD-HIT program. Modify only the text after the equal (=) sign.

```
\\Specify the run name**
run_name = MODIFY [Alphanumeric; e.g. New]
\\Directory of the genomes. Don't put a slash at the end.**
dir_genomes = MODIFY [Alphanumeric; e.g. GenomesNCBI2012]
\\Directory of the RFAM models. Don't put a slash at the end.**
dir_RFAM = MODIFY [Alphanumeric; e.g. RFAMdb]
\\Directory of the CD-HIT program. Don't put a slash at the end.**
dir_CDHIT = MODIFY [Alphanumeric; e.g. cd-hit-v4]
\\Run BLAST for all genomes. Only do it the first time you run the pipeline.**
run_blast = MODIFY [T/F; e.g. T]
\\Number of CPUs to use in the BLAST. Select one if working in your PC, or more if you're working in a cluster.**
cpus = MODIFY [1 for PC or >1 when working in a cluster; e.g. 8]
```

The rest of the parameters can be left as they are or modified according to the user’s preferences:


```
\\Minimum promoter size to be extracted
min_prom = 100
\\Minimum CDS size to be extracted
min_cds = 100
\\Max promoter region taken for doing BLAST
max_prom = 100
\\Max CDS region taken for doing BLAST
max_cds = 100
\\Set-off position for taking promoter and CDS regions for doing BLAST
set_off = 50
\\Minimum word size in BLAST
word_size = 11
\\Maximum e-value allowed in the BLAST
max_eval = 0.0001
\\BLAST score threshold for defining clusters in NetClust
score = 0
\\Minimum length of BLAST alignment for defining same promoter
min_length = 50
\\Minimum identity of BLAST alignment for defining same promoter
min_identity = 80
\\Percentage of allowed CDS in each cluster
max_cds_cluster = 0.2
\\Percentage of allowed paralogs in each cluster
max_paralogs = 0.2
\\Percentage of allowed upstream paralogs in each cluster
max_upparalogs = 0.2
\\Minimum length of alignment for defining same cluster in CD-HIT (use same as BLAST)
aL = 0.5
\\Identity level for defining same cluster in CD-HIT (use same as BLAST)
c = 0.8
```


5. Run master Perl script located in pipeline/

```
perl pipeline_PMP.pl parameters.txt
``

	a. _pipeline_PMP.pl_
	Script that coordinates the complete pipeline. To establish/modify the settings of the pipeline, you need to modify the file ‘parameteres.txt’ which is the input of this program (see below in 5). Next a brief explanation of the scripts called by pipeline_PMP.pl:
	
	b. _getCDSandPromRegions.pl_
	This script extracts the promoters and CDSs sequences from each genome.
	
	c. _blastClust.pl_	
	It calls for Netclust in order to produce inter-genome clusters of PMPs.
	
	d. _filterParalogs.pl_
	Filter for paralogy of the downstream CDSs of the clusters of PMPs.
	
	e. _filterUpstreamParalogs.pl_
	Filter for paralogy of the upstream CDSs of the clusters of PMPs.
	
	f. _getRep.pl_
	Get representatives and all members of each cluster of PMPs.
	
	g. _formatDatabase.pl_
	Add information to each cluster, such as functional annotation of the downstream CDS of the PMP.
	
	h. _formatcdhit.pl_
	Format file of representatives for using CD-HIT.
	
	i. _finalFormat.pl_
	Add a cluster ID to the non-redundant set of PMPs clusters.

