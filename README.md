# assembly_stats
The script provides some statistics for assemblies. It will ask for the following files:

File containing PE reads (FASTA format)
File containing scaffolds  (FASTA format)
File containing scaffolds that are longer than 1kb (FASTA format)

The output is placed in a file named assembly_stats.csv and contains the following values:
Fields 1-7 refer to all assembled scaffolds
Fields 8-14 refer to scaffodls longer than 1kb

1. Number of bp in assembly
2. Number of reads used for the assembly
3. Number of bp in all assembled scaffolds
4. Percent of bp that were assembled (calculated from items 1 and 3)
5. Number of scaffolds that were assembled
6. Average scaffold length
7. The length of the longest scaffold
8. N50
9. Number of bp found in scaffolds that are longer than 1kb
10. Percent of bp that were assembled into contigs longer than 1kb (calculated from items 1 and 9)
11. Number of scaffolds that are longer than 1kb
12. The average length of scaffolds that are longer than 1kb
13. N50 of scaffolds that are longer than 1kb
14. Percent of scaffolds that are longer than 1kb (calculated from 11 and 5)


requires:
Python 3.x, Biopython, numpy 
