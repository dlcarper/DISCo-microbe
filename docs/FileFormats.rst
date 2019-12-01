==========================
File Format Descriptions
==========================

Create module
==============

**Input files**

--i-alignment: This is an alignment of all sequences the user would like to evaluate in FASTA format

:Example: 
  >S003715306
  
  GCGGTA-AT-ACGTAG-GGAGCAAGCGTTGTC-CGG-ATTTATTGG-GCGTAAA-GAGCTCGTAG-G-CGGCTT-GGCAAGT-CGGATGTGAAA-CC-CCCAGG-CTTAACC-TGGGG-C--C- 	GCCATTCGA-TAC-TGC-TATGG-C-TT-GAGTTCGGTA-GGGGAT-TG-TGGA-ATT-CC-C-GGTGTAGCGGTGAAATGCGCAG-ATATCG-GG-AGGA-ACACC-AATG-GCGAAGGCAG- 	CAAT-CTGGGC-CGACACT-GA-CGCTGA-GG--A-GCGAAA-GCGTGGG-G-AGCAAA-CAGGATTAGATA
  
  >S003614093
  
  GCGGTA-AT-ACGTAG-GGAGCAAGCGTTGTC-CGG-AATTATTGG-GCGTAAA-GAGCTCGTAG-G-CGGTTC-GGTAAGT-CGGGTGTGAAA-AC-TCAAGG-CTCAACC-TTGAG-A--C-	GCCACTCGA-TAC-TGC-CGTGA-C-TT-GAGTCCGGTA-GAGGAG-TG-TGGA-ATT-CC-T-GGTGTAGCGGTGAAATGCGCAG-ATATCA-GG-AGGA-ACACC-AGCG-GCGAAGGCGG-	CACT-CTGGGC-CGGTACT-GA-CGCTGA-GG--A-GCGAAA-GCATGGG-G-AGCAAA-CAGGATTAGATA
  
  >S001611178
  
  GCGGTA-AT-ACGTAG-GGCGCGAGCGTTGTC-CGG-AATTATTGG-GCGTAAA-GGGCTCGTAG-G-CGGCTT-GTTGCGC-CTGCTGTGAAA-AC-GCGGGG-CTTAACT-CCGCG-C-GT-	G-CAGTGGG-TAC-GGG-CA-GG-C-TT-GAGTGTGGTA-GGGGTG-AC-TGGA-ATT-CC-A-GGTGTAGCGGTGGAATGCGCAG-ATATCT-GG-AGGA-ACACCGAT-G-GCGAAGGCAG-	GTCA-CTGGGC-CATTACT-GA-CGCTGA-GG--A-GCGAAA-GCGTGGG-T-AGCGAA-CAGGATTAGATA

--p-include-strains: A list of community members the community the user would like added. Each sequence identifier (must match what is in alignment)is on its own line

:Example:
  S003715306
  
  S003614093
  
  S001611178

--i-metadata: Information to combine with the community output. File must contain a header, be tab-delimited, and contain the identifiers in the first column

:Example:
  ID	Phylum	Class
  
  S003715306	Actinobacteria	Actinobacteria
  
  S003614093	Actinobacteria	Actinobacteria
  
  S001611178	Actinobacteria	Actinobacteria
  
  S000014419	Actinobacteria	Actinobacteria

--i-distance-database: Pre-calculated distance database of sequences, this is created when a previous create command has been run with the same strains. It is a tab delimited file with each line comparing two sequences. The sequence identifiers are in the first two columns and the edit distance between them is in the third

:Example:
  S003715306	S003614093	28
  
  S003715306	S001611178	50
  
  S003715306	S000014419	48
  
  S003715306	S000015295	49
  
  S003715306	S000022350	42
  
  S003715306	S000129061	44

**Output files**

 --o-community-list: A tab delimited list of strains, with each strain on its own line with a header line. If metadata is supplied it will be combined with this output

:Example:
  ID	Phylum	Class
  
  S003715306	Actinobacteria	Actinobacteria
  
  S003614093	Actinobacteria	Actinobacteria
  
  S001611178	Actinobacteria	Actinobacteria
  
  S000014419	Actinobacteria	Actinobacteria
  
--o-fasta: A FASTA file containing only the strains in the constructed community

Subsample module
================

**Input files**

 --i-input-community: Tab seperated file with taxa ids in the first column with metadata in additional columns, output of create module
 
:Example:
  ID	Phylum	Class
  
  S003715306	Actinobacteria	Actinobacteria
  
  S003614093	Actinobacteria	Actinobacteria
  
  S001611178	Actinobacteria	Actinobacteria
  
  S000014419	Actinobacteria	Actinobacteria

--p-proportion: File of the relative proportions of each taxonomic rank desired in final community. Each rank is contained on its own line. The rank and the proportion are seperated by a tab.

:Example:
  Actinobacteria	0.1
  
  Aquificae	0.001
  
  Bacteroidia	0.05
  
  Flavobacteriia	0.001
  
  Sphingobacteriia	0.003

**Output files**

A file with each sequence identifier on its own lines for the subsampled community

