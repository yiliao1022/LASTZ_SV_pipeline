#!/bin/bash

export REF=/data/users/liaoy12/liaoy12/TAD/UsingNaoassemblies/ref/Dere.fa
export Sam=/data/users/liaoy12/liaoy12/TAD/genomes/A1.fa
export Rname="Target"
export Qname="Query"

faSize -detailed $REF > $REF.sizes
faSize -detailed $Sam > $Sam.sizes
faToTwoBit $REF $REF.2bit
faToTwoBit $Sam $Sam.2bit

for i in ./*.lz; do axtChain -linearGap=/data/users/liaoy12/Pipeplines/linux.x86_64/medium $i $REF.2bit $Sam.2bit $i.chain; done
chainMergeSort *chain > all.chain
chainPreNet all.chain $REF.sizes $Sam.sizes all.chain.filter
chainNet -minSpace=1 all.chain.filter $REF.sizes $Sam.sizes all.chain.filter.tnet all.chain.filter.qnet
netSyntenic all.chain.filter.tnet all.chain.filter.tnet.synnet
netSyntenic all.chain.filter.qnet all.chain.filter.qnet.synnet
netToAxt all.chain.filter.tnet.synnet all.chain.filter $REF.2bit $Sam.2bit all.chain.filter.tnet.synnet.axt
axtToMaf all.chain.filter.tnet.synnet.axt $REF.sizes $Sam.sizes -tPrefix=$Rname. -qPrefix=$Qname. all.chain.filter.tnet.synnet.axt.maf

cat all.chain.filter.tnet.synnet | grep 'net\|^ fill\|^  gap' > synnet1.txt
cat all.chain.filter.tnet.synnet | grep 'net\|^  gap\|^   fill' > synnet2.txt
cat all.chain.filter.tnet.synnet | grep 'net\|^   fill\|^    gap' > synnet3.txt
cat all.chain.filter.tnet.synnet | grep 'net\|^    gap\|^     fill' > synnet4.txt
cat all.chain.filter.tnet.synnet | grep 'net\|top' > Net_top.txt

sed -i 's/^ fill/fill/g' Net_top.txt
sed -i 's/net/fillnet/g' synnet1.txt
sed -i 's/^ fill/fill/g' synnet1.txt
sed -i 's/^  gap/gap/g' synnet1.txt
sed -i 's/net/gapnet/g' synnet2.txt
sed -i 's/^  gap/gap/g' synnet2.txt
sed -i 's/^   fill/fill/g' synnet2.txt
sed -i 's/net/fillnet/g' synnet3.txt
sed -i 's/^   fill/fill/g' synnet3.txt
sed -i 's/^    gap/gap/g' synnet3.txt
sed -i 's/net/gapnet/g' synnet4.txt
sed -i 's/^    gap/gap/g' synnet4.txt
sed -i 's/^     fill/fill/g' synnet4.txt

perl ~/Pipeplines/Lastz_SV/SVcallerl2.pl -input synnet2.txt -output level2 -refname $Rname -queryname $Qname -refseq $REF
perl ~/Pipeplines/Lastz_SV/SVcallerl1.pl -input synnet3.txt -output level3 -refname $Rname -queryname $Qname -refseq $REF
perl ~/Pipeplines/Lastz_SV/SVcallerl2.pl -input synnet4.txt -output level4 -refname $Rname -queryname $Qname -refseq $REF
perl ~/Pipeplines/Lastz_SV/SVcallerl1.pl -input synnet1.txt -output level1 -refname $Rname -queryname $Qname -refseq $REF

perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level1.Complex.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level1.Deletion.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level2.CNV.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level2.INV.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level3.Complex.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level3.Deletion.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level4.CNV.txt
perl ~/Pipeplines/Lastz_SV/filter.pl Net_top.txt level4.INV.txt

cat level1.Complex.txt.synteny.out level3.Complex.txt.synteny.out > $Rname.$Qname.complex.sv
cat level1.Deletion.txt.synteny.out level3.Deletion.txt.synteny.out > $Rname.$Qname.deletion.sv
cat level2.CNV.txt.synteny.out level4.CNV.txt.synteny.out > $Rname.$Qname.CNV.sv
cat level2.INV.txt.synteny.out level4.INV.txt.synteny.out > $Rname.$Qname.INV.sv

rm *.txt *.out
mkdir SVnew
mv *.sv SVnew



