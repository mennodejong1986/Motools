# Motools

Motools, named after Motoo Kimura, is a small tool for investigating the accuracy of K2P-corrected distance estimates. Motools calculates the fit between observed and predicted estimates of transition and transversion-type difference between pairs of individuals (of from different species). 

The Motools pipeline expects as input a filtered bam-file, which contains the mapping information of sequencing reads relative to the reference genome of a relatively closely related species. Minimum required average read depth is roughly 8x, because Motools will return imprecise estimates if heterozygous sites cannot be determined accurately.

Users of Motools can call genotypes using their own preferred caller, or make use of the MOTOOLS_bam2vcf.sh script, which runs bcftools mpileup, bcftools call and bcftools filter. If using another method/caller, make sure to use a software which outputs a gVCF-file (i.e., a vcf file containing genotype calls for all sites, including monomorphic (homozygous) sites).

If using the MOTOOLS_bam2vcf.sh script, users need to edit the user-defined section (e.g. define the paths to bcftools and samtools, and to the bam file and reference genome). To run, make the script an executable (chmod +x MCART_bam2vcf.sh script), and afterwards execute:

./MOTOOLS_bam2vcf.sh

To subsequently count the numbers and proportions of single nucleotide variants (SNVs), divided over transition-type and transversion-type differences, execute:

./MOTOOLS_countTsTv.sh

This script will generate a file called ‘tstv_counts.txt’.

To estimate TMRCA and generate output plots, execute in R:

source(“MOTOOLS.txt”)

run_motools(tstv_file=‘tstv_counts.txt’,s_0=0.9987,mut_rate=2.5*10^-8, ngen=400000)
