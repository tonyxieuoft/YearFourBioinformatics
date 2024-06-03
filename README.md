# Requirements

The following must be installed:
- Selenium for Python with ChromeDriver (https://www.selenium.dev/downloads/, https://sites.google.com/chromium.org/driver/downloads?authuser=0)
- BioPython (https://biopython.org/wiki/Download)

More specifically, BioPython is required to access the NCBI Entrez API, and Selenium is required for the web-driver based automation of NCBI BLAST.

# Program Overview

The program offers automated functionality for two major use cases:
- Pulling exons for well-annotated genes and taxa from the NCBI Gene database
= Blasting gene exons of reference species against non-well annotated genomes to obtain full sets of sequences for given taxa.

When running the main program, the following will occur sequentially:
1) The program first asks the user for basic information (username, directory to download to)
2) Given user-inputted lists of genes and taxa, the program pulls out all available exon information from NCBI and directs it to an output file (skippable)
3) The program organizes seuqences pulled out in step 2 (or during a previous iteration of the program) into query files in preparation for BLAST.
4) Using the query files assembled in step 3, the program runs the NCBI BLAST program against genomes in the NCBI database for given taxa.
5) Finally, the program concatenates the results and outputs alignments by gene.

Details for steps 2-5 are provided below. 

# Pulling exons from NCBI

The user provides two files on hand for this part of the program: one that lists taxa to pull sequences for, and another that lists queries that the program will use to search for specific genes. 

### Taxa file

The taxa file should be in .txt format, and must be organized the following way:
```
taxon1
taxon2
taxon3
...
```
Note that there are no headers. Each line contains one taxon. 

### Gene query file

The structure of input for the gene query file is more involved, and the file should be in either .txt or .csv format. As a .txt file, it must be organized the following way:
```
marker: gene1query1, marker: gene1query2,   ...
marker: gene2query1, marker: gene2query2,   ...
marker: gene3query1, marker: gene3:query2,  ...
...
```
Again, there are no headers. All queries for a gene are on the same line, separated by commas similar to a comma-separated-file format (.csv). Before each query is a marker denoting the type of query. Markers are separated from the queries they denote via colons (':'). 

The following markers are available:
- `g`: indicates that a query is an abbreviated gene name (eg. 'RHO', 'GRK7').  
- `d`: indicates that a query is a gene description. (eg. 'rhodopsin', 'G protein-coupled receptor kinase 7')

 For genes are known under multiple possible abbreviated names or description, multiple instances of the same marker can be used in a given line.

 ### Output directory

 The results pulled out from the program are outputted in a layer of nested folders with the following structure: `General Folder -> Gene -> Taxon -> Species`. Each species folder contains fasta files corresponding to different transcript versions. 

 
