## Colin D., 2018 and 2020
## Splits a big multifasta into individual fasta files so they can be recombined
## Only extracts the files in the fa.fai (so you can delete some lines out of this,
## and have the genomes excluded from the out directory
## Uses a fasta index and samtools


#gen=refSeqs_allKingdoms_2020_03.fa
#gen=selectedgenomes.fa
#gen=all_1_cln_1strainperspecies.fa
gen=all_bact.fa.masked.fa

mkdir out

for i in `cut -f1 $gen.fai`
	do
	samtools faidx $gen $i > out/$i.fa
done

