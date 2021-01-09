# mazzoni-begendiv
This repository contains scripts developed by members of Mazzoni group at IZW/BeGenDiv 

**Scripts are identified as   
GenNGS: scripts for handling NGS data in general;   
RAD: scripts for handling ddRAD/3RAD data;   
RefGenome: scripts for handling reference genome data and quality control)**   

### dereplicateFQ

```
usage: dereplicateFQ.py [-h] [--g] fastqIn out

dereplicates reads in a fastq file, takes the highest quality value for each
base

positional arguments:
  fastqIn     fastq input file
  out         fastq output file(dereplicated)

optional arguments:
  -h, --help  show this help message and exit
  --g, -gzip
```

### RAD_digestion

```
usage: RAD_digestion.py [-h] -g genomeFile -e enzymeFile [--dd]
                        [-o outputFile] [--q] [--min MIN] [--max MAX] [--gz]
                        [--v] [--rad]

In-silico RAD digestion.

optional arguments:
  -h, --help      show this help message and exit
  -g genomeFile   Sequences to digest in fasta or fastq format (f.e. a genome file).
  -e enzymeFile   List with enzymes/enzyme combinations including restriction pattern. Tab-separated.
                  Format:
                  enzyme1	pattern1	enzyme2	pattern2 -- the first enzyme is the forward enzyme, the second one in the list is reverse (format is as shown in the file downloaded by the get_rad_enzymes.py script)
  --dd            If argument is set, only fragments cut by two different enzymes will be returned.
  -o outputFile   Output file prefix. [Default: stdout]
  --q, -fastq     Genome file is in fastq format. Output will also be in fastq-format (qual. Phred33)
  --min MIN       Minimum size for fragments.
  --max MAX       Maximum size for fragments.
  --gz, -gzip     Enable if input file is gzipped.
  --v, -verbose   Enables progress information in console.
  --rad, -radseq  Allow fragments to have Start(=first/forward enzyme in list) and End(=second/reverse enzyme in list) of a scaffold as a restriction side. Enable if digesting already digested sequences.
```

### checkRestrictionSites

```
usage: checkRestrictionSites.py [-h] fqIn outCorrect info res1 res2

filters a fastq file to start and end with defined sequences (designed for
merged ddRAD data)

positional arguments:
  fqIn        fastq input file
  outCorrect  output file for sequences with correct restriction sites
  info        output file for numebr of reads with correct or incorrect
              restriction sites
  res1        overhang of 1. enzmye --> starting basepairs for all reads
  res2        overhang of 2. enzmye --> ending basepairs for all reads

optional arguments:
  -h, --help  show this help message and exit
```

### cluserFromPairs

Script to greedily construct single linkage clusters based on similarity values from pair-wise comparison with `vserach`. Clusters are identified by finding connect components in a similarity graph using the [NetworkX](https://networkx.github.io/) package.


```
usage: clusterFromPairs.py [-h] [-vs VSEARCH] [-fq FASTQ]

builds clusters based on pairwise identites from vsearch

optional arguments:
  -h, --help            show this help message and exit
  -vs VSEARCH, --vsearch VSEARCH
                        vsearch pairwise identity tabout file
  -fq FASTQ, --fastq FASTQ
                        fastq file used as an input for vsearch
```
