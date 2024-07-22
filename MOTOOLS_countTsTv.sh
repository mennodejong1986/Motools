#!/bin/bash
# This script calculates proportion of single nucleotide variants and Ts/Tv ratio

# An alternative could be to use vcftools --TsTv-summary or --FILTER-summary
# However, vcftools applies some additional filters, resulted in a (slightly) different estimate

####################################################
MYVCFGZ=input.vcf.gz
####################################################

echo "Unzipping vcf file..."
zcat $MYVCFGZ > mytemp.vcf

echo "Counting total number of sites..."
grep -v '#' mytemp.vcf | wc -l > myoverallcounts.txt

echo "Selecting SNVs..."
# 0/0 = homozygous reference
# 0/1 = heterozygous (one allele as reference, one allele alternative)
# 1/1 = homozygous alternative
grep '1/1' mytemp.vcf > myhomo.txt
grep '0/1' mytemp.vcf > myhetero.txt

echo "Extracting transitions..."
# transition is A <=> G and C <=> T:
cat myhomo.txt | awk '$4 == "A" && $5 == "G" || $4 == "G" && $5 == "A" || $4 == "C" && $5 == "T" || $4 == "T" && $5 == "C"' > myhomo_transitions.txt
cat myhetero.txt | awk '$4 == "A" && $5 == "G" || $4 == "G" && $5 == "A" || $4 == "C" && $5 == "T" || $4 == "T" && $5 == "C"' > myhetero_transitions.txt

echo "Extracting ambiguous..."
cut -f5 myhomo.txt | grep ',' | sed 's/,/\t/g' > myhomo_ambiguous.txt
cut -f5 myhetero.txt | grep ',' | sed 's/,/\t/g' > myhetero_ambiguous.txt
cat myhomo_ambiguous.txt | awk '$1 == "A" && $2 == "G" || $1 == "G" && $2 == "A" || $1 == "C" && $2 == "T" || $1 == "T" && $2 == "C"' > myhomo_ambiguous_Ts.txt
cat myhetero_ambiguous.txt | awk '$1 == "A" && $2 == "G" || $1 == "G" && $2 == "A" || $1 == "C" && $2 == "T" || $1 == "T" && $2 == "C"' > myhetero_ambiguous_Ts.txt
wc -l myhetero_ambiguous_Ts.txt myhomo_ambiguous_Ts.txt | grep -v 'total' | sed 's/^ *//g' | sed 's/myhetero_ambiguous_Ts.txt/hetero.amb.Ts/g' | sed 's/myhomo_ambiguous_Ts.txt/fixed.amb.Ts/g' | sed 's/ /\t/g' > myambiguous.txt

echo "Combining results..."
wc -l myhomo.txt myhetero.txt | grep 'total' | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
wc -l myhetero.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
echo "NA" >> myoverallcounts.txt
wc -l myhetero_transitions.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
wc -l myhetero_ambiguous.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
wc -l myhomo.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
echo "NA" >> myoverallcounts.txt
wc -l myhomo_transitions.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt
wc -l myhomo_ambiguous.txt | sed 's/^ *//g' | cut -f1 -d ' ' >> myoverallcounts.txt

tr '\n' '\t' < myoverallcounts.txt | awk '{ $4 = $3 - $5 - $6 } { $8 = $7 - $9 - $10 } { $11 = $5/$4 } { $12 = $9/$8 } { $13 = ($5 + $9)/($4 + $8) } { $14 = ($5 + $9)/($5 + $4 + $9 + $8) }   1 ' | tr ' ' '\n' > mytemp.txt
echo "sites_SNVs_hetero_hetero.Tv_hetero.Ts_hetero.amb_fixed_fixed.Tv_fixed.Ts_fixed.amb_Ts/Tv.hetero_Ts/Tv.fixed_Ts/Tv.all_Ts.prop" | tr '_' '\n' > myrownames.txt
paste mytemp.txt myrownames.txt > tstv_counts.txt
cat tstv_counts.txt myambiguous.txt > mytemp.txt
mv mytemp.txt tstv_counts.txt

rm mytemp.vcf myambiguous.txt myhomo_ambiguous_Ts.txt myhetero_ambiguous_Ts.txt myhomo_transitions.txt myhomo.txt myhomo_ambiguous.txt myhetero.txt myhetero_transitions.txt myhetero_ambiguous.txt myoverallcounts.txt myrownames.txt

echo " "
echo "Results have been stored in file 'tstv_counts.txt'"
echo "Output table:"
cat tstv_counts.txt