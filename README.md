# mazzoni-begendiv
This repository contains scripts developed by members of Mazzoni group at IZW/BeGenDiv 

**Scripts are identified as   
NGSdata: scripts for handling NGS data in general;   
RAD: scripts for handling ddRAD/3RAD data;   
RefGenome: scripts for handling reference genome data and quality control)**   

### dereplicateFQ (NGSdata)

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

### Filter_Reads.py (NGSdata)

```
usage = %prog [options] input.fasta [output.fasta]

filter fasta or fastq reads based on a list of IDs, minimum size or randomly

positional arguments:
  input.fasta   fasta input file
  output.fasta  fasta output file

optional arguments:
-q, --quiet  do not print status messages to the screen
-w, --supress-warnings  do not print warnings if a ID was found more or less than once
-u, --fastq  input file is fastq
-z, --gzip input file is gzipped
-l, --min-length  write only sequence with lengths at least X
-i, --id-list  write only sequence with an ID from this list. 
-r, --random   randomly sample X sequence from input file
-e, --regexp   use regular expression instead of exact matching for IDs
-a, --ignore-at  ignore the first letter of the query IDs if it is an @ 
-n, --negative  do exactly the opposite of what would normally be done
```

### RAD_digestion (RAD)

```
usage: RAD_digestion.py [-h] -g genomeFile -e enzymeFile [--dd]
                        [-o outputFile] [--q] [--min MIN] [--max MAX] [--gz]
                        [--v] [--rad]

In-silico RAD digestion of a fasta genome or from merged paired-end reads from ddRAD-like data

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

### checkRestrictionSites (RAD)

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

### clusterFromPairs (NGSdata or RAD)

Script to greedily construct single linkage clusters based on similarity values from pairwise comparison with `vsearch`. Clusters are identified by finding connect components in a similarity graph using the [NetworkX](https://networkx.github.io/) package.


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

### Create_Haplotype_Structure.py (RAD)

```
usage: Create_Haplotype_Structure.py popMap mappedCounts popOutDir outDir 

converts the output from Stacks populations into a Structure file based on haplotypes (rather than SNPs), adding new filtering options

positional arguments:
  popMap        population map used for stacks
  mappedCounts  tab separated file containing the readcounts per locus
  popOutDir     output directory of the stacks population module
  outDir        output directory

optional arguments:
  --readCOV  minimum read coverage needed for a locus [default=6]
  --popCOV   total number of populations that needs to cover a locus to include it in the structure output [default=4]
  --intraCOV  percentage of each population that needs to cover a locus to include it in the structure output [default=0.7]
  --intraCOVtotal  enable this option to have a strict int value as input for --intraCOV 
```
