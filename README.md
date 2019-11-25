# README

## Description

This pipeline can be used to roughly detect structral variations, such as insertions, deletions, copy Number Variations (i.e. mostly tandem duplications),inversions,and those with breakpoints unambiguously being defined, using high-quality chromosome-level assemblies between the genomes of closely related species or within species. It is designed mainly based on the Lastz/Chain/Net tools. Currently, It works well for rice (Oryza sativa), fruitfly (Drosophila melanogaster) and Brassica genomes. 

## Dependency 

LASTZ: You can download from http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html. 

RepeatMasker: You can download from http://www.repeatmasker.org.

KentUtils: You can download from http://hgdownload.soe.ucsc.edu/admin/exe/.

## How to use  

### Firstly, you need the soft-masked genome sequences for both reference and query. This can significantly reduce the runing time for Lastz alignment. 

export REF=/your/paths/to/reference
export QUERY=/your/paths/to/query

RepeatMasker -e crossmatch -s -pa 64 -xsmall -gff -lib $TE $reference
RepeatMasker -e crossmatch -s -pa 64 -xsmall -gff -lib $TE $query
 
### Second, performing the pairwise alignment of reference and query sequences.

Lastz $reference $query --format=axt --output=$reference.lz

### Third, running the SV discovery steps with the bash script.

export TARGET_name=target
export Query_name=query

sh LastzSV.sh

## Output format

Target "\t" Start "\t" End "\t" SV_type "\t" Query "\t" Start "\t" End "\t" strand "\t" Forward:Backward bases "\t" sequence in target

Dere.utg000001l 32      36      Deletion        ORE.2L  21422524        21422524        4       +       29:15   AAAC
Dere.utg000001l 51      53      Deletion        ORE.2L  21422539        21422539        2       +       15:120  TA
Dere.utg000001l 173     175     Deletion        ORE.2L  21422663        21422663        2       +       120:165 GG
Dere.utg000001l 340     344     Deletion        ORE.2L  21422830        21422830        4       +       165:47  aaaa

## Further effort

Further efforts will include: 1) Inter/Intra-Chromosomal duplications and other complex SVs; 2) genotyping SVs with short or long reads methods; 3) multiple genome comparsion.

## Contact and citation

If you have any question, please feel free to contact yiliao1022@gmail.com
