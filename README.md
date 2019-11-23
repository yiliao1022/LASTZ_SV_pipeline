# README

## Description

This pipeline can be used to roughly detect structral variations, such as insertions, deletions, copy Number Variations (i.e. mostly tandem duplications),inversions,and those with breakpoints unambiguously being defined ones, using high-quality chromosome-level assemblies between the genomes of closely related species or within species. It is designed mainly based on the Lastz/Chain/Net tools. Currently, It works well for rice (Oryza sativa), fruitfly (Drosophila melanogaster) and Brassica genomes. Further efforts needed to address the following targets: 1) Inter/Intra-Chromosomal duplications and other complex SVs; 2) genotyping SVs with short or long reads methods; 3) multiple genome comparsion.

## Dependency 

LASTZ: You can download from http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html;
RepeatMasker: You can download from http://www.repeatmasker.org;
KentUtils: You can download from http://hgdownload.soe.ucsc.edu/admin/exe/.

## HOW TO USE  

Firstly, you need the soft-masked genome sequences for both reference and query. This can significantly reduce the runing time for Lastz alignment. 
