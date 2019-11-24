# README

## Description

This pipeline can be used to roughly detect structral variations, such as insertions, deletions, copy Number Variations (i.e. mostly tandem duplications),inversions,and those with breakpoints unambiguously being defined ones, using high-quality chromosome-level assemblies between the genomes of closely related species or within species. It is designed mainly based on the Lastz/Chain/Net tools. Currently, It works well for rice (Oryza sativa), fruitfly (Drosophila melanogaster) and Brassica genomes. 

## Dependency 

LASTZ: You can download from http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html. 

RepeatMasker: You can download from http://www.repeatmasker.org.

KentUtils: You can download from http://hgdownload.soe.ucsc.edu/admin/exe/.

## How to use  

### Firstly, you need the soft-masked genome sequences for both reference and query. This can significantly reduce the runing time for Lastz alignment. 

RepeatMasker -e crossmatch -s -pa 64 -xsmall -gff -lib $TE $reference

RepeatMasker -e crossmatch -s -pa 64 -xsmall -gff -lib $TE $query
 

### Second, performing the pairwise alignment of reference and query sequences.


### Third, running the SV discovery steps with the bash script.



## Further efforts

Further efforts will include: 1) Inter/Intra-Chromosomal duplications and other complex SVs; 2) genotyping SVs with short or long reads methods; 3) multiple genome comparsion.

## Contact and citations

If you have any question, please feel free to contact yiliao1022@gmail.com
