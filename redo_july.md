Keep the same:
```
echo 'export PATH=/home/aygera/tools/array-analysis-cli-linux-x64-v2.1.0/array-analysis-cli/:$PATH' >> ~/.bashrc
source ~/.bashrc
array-analysis-cli genotype call \
    --num-threads 32 \
    --idat-sample-sheet ./sample_sheet/dauren_gwas_sample_sheet.csv \
    --bpm-manifest ./manifest_bpm/InfiniumImmunoArray-24v2-0_A2.bpm    \
    --cluster-file ./cluster/InfiniumImmunoArray-24v2-0_A_ClusterFile.egt \    
    --output-folder ./gtc

# for some reason it doesn't save into gtc but rather in dauren_gwas; i manually moved resulting gtc files into gtc folder

array-analysis-cli genotype gtc-to-vcf \
    --bpm-manifest ./manifest_bpm/InfiniumImmunoArray-24v2-0_A2.bpm \
    --genome-fasta-file /home/aygera/biostar/dauren_gwas/GRCh38_genome/GRCh38_genome.fa \
    --gtc-sample-sheet ./sample_sheet/dauren_gwas_sample_sheet.csv \
    --csv-manifest ./manifest_csv/InfiniumImmunoArray-24v2-0_A2.csv \
    --output-folder ./vcf

bgzip vcf/*vcf*
for f in vcf/*.vcf.gz; do tabix -p vcf -f $f;done
bgzip dauren_gwas.vcf
bcftools merge /home/aygera/biostar/dauren_gwas/vcf/*.vcf.gz -o dauren_gwas.vcf
```

Starting anew with dauren_gwas.vcf and phenotypes.txt (i fixed it so 2023 has cases and 2024 controls and made sure the namings are consistent)

1) checking key stats
conda create -n bcftools_env -c bioconda -c conda-forge bcftools -y
conda activate bcftools_env
bcftools stats ./dauren_gwas.vcf.gz
```
# This file was produced by bcftools stats (1.24+htslib-1.24) and can be plotted using plot-vcfstats.
# The command line was:	bcftools stats  ./dauren_gwas.vcf.gz
#
# Definition of sets:
# ID	[2]id	[3]tab-separated file names
ID	0	./dauren_gwas.vcf.gz
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	192
SN	0	number of records:	247814
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	246554
SN	0	number of MNPs:	0
SN	0	number of indels:	1260
SN	0	number of others:	0
SN	0	number of multiallelic sites:	515
SN	0	number of multiallelic SNP sites:	515
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV	[2]id	[3]ts	[4]tv	[5]ts/tv	[6]ts (1st ALT)	[7]tv (1st ALT)	[8]ts/tv (1st ALT)
TSTV	0	177742	69347	2.56	177673	68881	2.58
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS	[2]id	[3]allele count	[4]number of SNPs	[5]number of transitions	[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
SiS	0	1	247089	177742	69347	1260	0	0	1260
# AF, Stats by non-reference allele frequency:
# AF	[2]id	[3]allele frequency	[4]number of SNPs	[5]number of transitions	[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
AF	0	0.000000	247089	177742	69347	1260	0	0	1260
# QUAL, Stats by quality
# QUAL	[2]id	[3]Quality	[4]number of SNPs	[5]number of transitions (1st ALT)	[6]number of transversions (1st ALT)	[7]number of indels
QUAL	0	.	246554	177673	68881	1260
# IDD, InDel distribution:
# IDD	[2]id	[3]length (deletions negative)	[4]number of sites	[5]number of genotypes	[6]mean VAF
IDD	0	-38	1	0	.
IDD	0	-32	1	0	.
IDD	0	-23	1	0	.
IDD	0	-20	1	0	.
IDD	0	-15	1	0	.
IDD	0	-14	4	0	.
IDD	0	-13	3	0	.
IDD	0	-12	4	0	.
IDD	0	-11	5	0	.
IDD	0	-10	5	0	.
IDD	0	-9	3	0	.
IDD	0	-8	7	0	.
IDD	0	-7	3	0	.
IDD	0	-6	11	0	.
IDD	0	-5	22	0	.
IDD	0	-4	87	0	.
IDD	0	-3	63	0	.
IDD	0	-2	119	0	.
IDD	0	-1	490	0	.
IDD	0	1	310	0	.
IDD	0	2	51	0	.
IDD	0	3	23	0	.
IDD	0	4	31	0	.
IDD	0	5	5	0	.
IDD	0	6	1	0	.
IDD	0	7	2	0	.
IDD	0	8	1	0	.
IDD	0	9	1	0	.
IDD	0	14	1	0	.
IDD	0	15	1	0	.
IDD	0	17	1	0	.
IDD	0	24	1	0	.
# ST, Substitution types:
# ST	[2]id	[3]type	[4]count
ST	0	A>C	9685
ST	0	A>G	40685
ST	0	A>T	5782
ST	0	C>A	11149
ST	0	C>G	7934
ST	0	C>T	48153
ST	0	G>A	48299
ST	0	G>C	7931
ST	0	G>T	10982
ST	0	T>A	5941
ST	0	T>C	40605
ST	0	T>G	9943
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP	[2]id	[3]bin	[4]number of genotypes	[5]fraction of genotypes (%)
```

*** so i need to split multiallelics (515) and keep an eye on indels (1260)

2) check ref/alt consistency
```
bcftools norm -c w -f ../GRCh38_genome/GRCh38_genome.fa ./dauren_gwas.vcf.gz -o /dev/null
# Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:	247814/0/0/85/0/0/0
```

3) I will realign those 85 snps and split multiallelics
Problem: I also need to convert my chip names into normal names and if i just split my multiallelics it might get broken. Need to find a solution
