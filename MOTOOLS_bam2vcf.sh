#!/bin/bash
# This script is for snp calling from bam file

########################  USER-DEFINED SECTION ##################################

NRSETS=20                                        # Number of subdivisions (this is basically the number of commands you want to run parallel)
SAMTOOLS=/path/to/samtools                       # Path to samtools executable
BCFTOOLS=/path/to/bcftools                       # Path to bcftools executable
REFERENCE=/path/to/referencegenome.fasta         # Path to reference genome
PREFIX=allsites                                  # Prefix of output bed file.
BAMFILES=mybamfiles.txt                          # a txt file with (path and) names of input bam file(s), one per line

createbed=TRUE
runmpileup=TRUE
combinebcf=TRUE
callgenotypes=TRUE
filtervcf=TRUE

#################################################################################

if [[ "$createbed" = TRUE ]]
then
echo "Start creating bed files"
${SAMTOOLS} faidx ${REFERENCE}                                                  # find out total number of bp in your genome
cut -f1-2 ${REFERENCE}.fai > contigs.bed
awk '{ total += $2 } END { print total }' contigs.bed > nrsites.txt             # I find: 2607875777. We are going to divide this over NRSETS equal batches (by generating bedfiles, which serve as input for samtools mpileup)
NRSITES=$(cat nrsites.txt)
cut -f2 contigs.bed > column2.txt
perl -lne 'print $sum+=$_' column2.txt > cum.txt
paste -d '\t' contigs.bed cum.txt > contigs.cum.bed

echo "Subdivisions:"
seq 1 $NRSETS > mynumbers.txt
for digit in $(cat mynumbers.txt)
do
num=$(($digit-1))
part=$((${NRSITES}/${NRSETS}))
var1=$((${num} * ${part}))
var2=$((${digit} * ${part}))
echo $digit
echo $var1
echo $var2
awk -v mystart="${var1}" -v myend="${var2}" -v mydigit="${digit}" '{ if ( $3 >= mystart && $3 <= myend ) print > "mybed"mydigit".txt" }' contigs.cum.bed
awk -v mystart="${var1}" -v myend="${var2}" '{ if ( $3 >= mystart && $3 <= myend ) print > "mybed.txt" }' contigs.cum.bed
cut -f1 mybed.txt > contignames.txt
cut -f2 mybed.txt > contiglengths.txt
sed -i 's/^/0___/' contiglengths.txt
sed -i 's/___/\t/g' contiglengths.txt
paste -d '\t' contignames.txt contiglengths.txt > mybed.txt
awk '{ if ($2 == "0") $2=1; print $0 }' mybed.txt | sed 's/ /\t/g' > mybed.temp.txt     # bcftools expects 1 based reference system
mv mybed.temp.txt mybed.txt                                                             # bcftools expects 1 based reference system
mv mybed.txt mybed${digit}.txt
done
rm mynumbers.txt
rm nrsites.txt
rm contig*
rm cum.txt
rm column2.txt
echo "Finished creating bed files."
fi

if [[ "$runmpileup" = TRUE ]]
then
echo "Starting bcftools mpileup now."
for mybedfile in mybed*txt
do
echo ${mybedfile}
${BCFTOOLS} mpileup -A -C50 -R ${mybedfile} --output-type b -f ${REFERENCE} --bam-list ${BAMFILES} > ${PREFIX}.${mybedfile}.bcf &
done
wait
echo "Finished mpileup."
fi
wait

if [[ "$combinebcf" = TRUE ]]
then
echo "Start concatenating bcf files to one single file."
${BCFTOOLS} concat --output-type b ${PREFIX}.mybed*.bcf > ${PREFIX}.combined.bcf
echo "Finished combining."
rm ./*mybed*bcf
rm ./mybed*.txt
echo "Removed intermediate bcf- and bed-files."
fi

if [[ "$callgenotypes" = TRUE ]]
then
echo "Genotype calling using bcftools call..."
echo "Crucially, excluding -v (--variant-sites) option, because info needed for both homozygous and heterozygous sites."
${BCFTOOLS} call -m --threads 10 --output-type z ${PREFIX}.combined.bcf > ${PREFIX}.vcf.gz
fi
wait

if [[ "$filtervcf" = TRUE ]]
then
echo "Filtering vcf file..."
$BCFTOOLS filter -i "QUAL>=30 && DP>=12 && DP<=40" ${PREFIX}.vcf.gz -O z | $BCFTOOLS view --exclude-types indels -O z -o allsites.vcf.gz
echo "Finished filtering vcf file."
fi
wait

echo "Finished analysis :-)"


############### EXPLANATION ##########

# You might think an easy way to split the files is simply by dividing the number of contigs evenly.
# However, some contigs can be bigger than others. If the big contigs are clumped, you can get a very uneven distribution.
# For example, his is what the distribution of sites I got when dividing the number of contigs by 4:
# As a result, some mpileup run would take might much longer than others, and time gained is low.
# Therefore, better to divide by number of sites, instead of by contig.

# the awk loop above execute these commands:
# awk '{ if ($4 >= 0*2607875777/16 && $4 <= 1*2607875777/16 ) print > "mybed1.txt" }' contigs.cum.bed
# awk '{ if ($4 >= 1*2607875777/16 && $4 <= 2*2607875777/16 ) print > "mybed2.txt" }' contigs.cum.bed
# etc
# until
# awk '{ if ($4 >= 15*2607875777/16 && $4 <= 16*2607875777/16 ) print > "mybed16.txt" }' contigs.cum.bed







