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
I will convert my chip names into normal names now while it's still vcf using dbsnp:

making plots of my vcf (2 not generated because its microarray)
```
bcftools stats -f ../GRCh38_genome/GRCh38_genome.fa \
  ./dauren_gwas.vcf.gz > my.stats.log
mkdir plots
plot-vcfstats -T "Variants" \
  -P -p ./plots/ my.stats.log
```

Annotation:
```
conda install bioconda::snpsift
bcftools index dauren_gwas.vcf.gz

wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi

bcftools index common_all_20180418.vcf.gz                         # might not be needed
SnpSift annotate common_all_20180418.vcf.gz dauren_gwas.vcf.gz > dauren_annotated.vcf
```

235193/247814 variants now have rsIDs, as well as gene names    

The file has this info now:
```
##INFO=<ID=WTD,Number=0,Type=Flag,Description="Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory">
##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description="First dbSNP Build for RS">
##INFO=<ID=COMMON,Number=1,Type=Integer,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">
##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
##INFO=<ID=TPA,Number=0,Type=Flag,Description="Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)">
##INFO=<ID=NOV,Number=0,Type=Flag,Description="Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common.">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">

1	1243468	1_1168711;rs115005664	G	A	.	PASS	ASP;CAF=0.9968,0.003195;COMMON=1;GENEINFO=C1QTNF12:388581|FAM132A:388581;HD;KGPhase1;KGPhase3;NSN;REF;RS=115005664;RSPOS=1243468;SAO=0;SSR=0;TOPMED=0.99568361365953109,0.00431638634046890;VC=SNV;VLD;VP=0x050000000605040436000100;WGT=1;dbSNPBuildID=132	GT:GS:BAF:LRR	0/0:0.8352:0.00521636:-0.121801	0/0:0.8352:0:0.0350399
```

Checking how many annotated:
```
bcftools view -H dauren_annotated.vcf | grep -c "RS="
#235193
bcftools view -H dauren_annotated.vcf | wc -l
#247814
```

Are the unannotated indels of SNVs?
```
bcftools view -v indels dauren_annotated.vcf -H | grep -vc "RS="
#308
bcftools view -v snps dauren_annotated.vcf -H | grep -vc "RS="
#12313
```

I have 1260 indels and 246554 snps so 308/1260 = 0.25 and 12313/246554 = 0.05


I might want to split my multiallelics before annotating:
```
bcftools norm -m -any -f /home/aygera/biostar/dauren_gwas/GRCh38_genome/GRCh38_genome.fa \
    dauren_gwas.vcf.gz -Oz -o dauren_gwas.split.vcf.gz
tabix -p vcf dauren_gwas.split.vcf.gz
```

Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:	247814/515/0/85/0/0/0


Re-annotating (might want to redo it with the entire dbsnp, not just common variants):
```
SnpSift annotate common_all_20180418.vcf.gz dauren_gwas.split.vcf.gz > dauren_annotated.vcf
bgzip dauren_annotated.vcf
tabix -p vcf dauren_annotated.vcf.gz
```

Re-check annotation rate:
```
bcftools view -H dauren_annotated.vcf.gz | wc -l
# 248349
bcftools view -H dauren_annotated.vcf.gz | grep -c "RS="
# 235283
bcftools view -v indels dauren_annotated.vcf.gz -H | grep -vc "RS="
# 229
bcftools view -v snps dauren_annotated.vcf.gz -H | grep -vc "RS="
# 12837
bcftools view -v indels dauren_annotated.vcf.gz -H | wc -l
1260
bcftools view -v snps dauren_annotated.vcf.gz -H | wc -l
247089

echo "scale=4; 229/1260" | bc
# .1817
echo "scale=4; 12837/247089" | bc
# .0519
```

Now ID cleanup and deduplication:
```
mkdir -p ./tmp_sort

# 1. Set ID to RS number where matched (chr_pos will be applied separately below for unmatched)
bcftools annotate --set-id '%INFO/RS' dauren_annotated.vcf.gz -Oz -o step1.vcf.gz
tabix -p vcf step1.vcf.gz

# Check how many still have ID="." (RS empty) vs proper rs numbers
bcftools query -f '%ID\n' step1.vcf.gz | grep -c "^\.$"

# 2. Split matched/unmatched, set chr_pos ID only on the unmatched subset, recombine.
#    (bcftools annotate -i/-e filters the WHOLE output, not just which rows get touched —
#     so we split first instead of trying to scope the --set-id call in place.)
bcftools view -e 'ID="."' step1.vcf.gz -Oz -o matched.vcf.gz
bcftools view -i 'ID="."' step1.vcf.gz -Oz -o unmatched.vcf.gz
tabix -p vcf matched.vcf.gz

bcftools annotate --set-id '%CHROM\_%POS' unmatched.vcf.gz -Oz -o unmatched.id.vcf.gz
tabix -p vcf unmatched.id.vcf.gz

bcftools concat matched.vcf.gz unmatched.id.vcf.gz -Oz -o step2_unsorted.vcf.gz
bcftools sort -T ./tmp_sort/ step2_unsorted.vcf.gz -Oz -o step2.vcf.gz
tabix -p vcf step2.vcf.gz

# sanity check: should equal step1's total record count
bcftools view -H step1.vcf.gz | wc -l
bcftools view -H step2.vcf.gz | wc -l

# 3. Drop withdrawn / non-unique-mapping / suspect variants
#    (WTD is Type=Flag, Number=0 — test presence directly, don't compare to =1)
bcftools view -e 'WTD || INFO/WGT>1 || INFO/SSR!=0' step2.vcf.gz -Oz -o step3.vcf.gz

# 4. Remove exact duplicate records (same CHROM:POS:REF:ALT)
bcftools norm -d exact step3.vcf.gz -Oz -o step4.vcf.gz

# 5. Check for duplicate IDs (different sites, same rsID — dbSNP merge artifacts)
bcftools query -f '%ID\n' step4.vcf.gz | sort | uniq -d > dup_ids.txt
wc -l dup_ids.txt

# 6. Remove ambiguous strand SNPs (A/T, C/G)
bcftools view -e '(REF="A" & ALT="T") | (REF="T" & ALT="A") | (REF="C" & ALT="G") | (REF="G" & ALT="C")' \
    step4.vcf.gz -Oz -o dauren_final.vcf.gz
tabix -p vcf dauren_final.vcf.gz

# 7. Final count log
echo "Split+annotated:      $(bcftools view -H dauren_annotated.vcf.gz | wc -l)"
echo "After ID fix (step2): $(bcftools view -H step2.vcf.gz | wc -l)"
echo "After suspect filter:  $(bcftools view -H step3.vcf.gz | wc -l)"
echo "After dedup:            $(bcftools view -H step4.vcf.gz | wc -l)"
echo "Final (no ambig strand):$(bcftools view -H dauren_final.vcf.gz | wc -l)"
```

Output:
```
13066
Checking the headers and starting positions of 2 files
Concatenating matched.vcf.gz	28.113457 seconds
Concatenating unmatched.id.vcf.gz
The chromosome block 1 is not contiguous, consider running with -a.
Writing to ./tmp_sort/SwOwp3
[W::bgzf_read_block] EOF marker is absent. The input may be truncated
[E::vcf_parse_format_check7] Number of columns at Y:22255793 does not match the number of samples (169 vs 192)
Error encountered while parsing the input
Cleaning
[W::bgzf_read_block] EOF marker is absent. The input may be truncated
248349
[W::bcf_sr_add_hreader] No BGZF EOF marker; file 'step2.vcf.gz' may be truncated
[E::vcf_parse_format_check7] Number of columns at Y:22255793 does not match the number of samples (169 vs 192)
Error: VCF parse error
235278
[W::bcf_sr_add_hreader] No BGZF EOF marker; file 'step2.vcf.gz' may be truncated
[E::vcf_parse_format_check7] Number of columns at Y:22255793 does not match the number of samples (169 vs 192)
Error: VCF parse error
Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:	235097/0/0/0/0/6/0
11 dup_ids.txt
Split+annotated:      248349
[W::bcf_sr_add_hreader] No BGZF EOF marker; file 'step2.vcf.gz' may be truncated
[E::vcf_parse_format_check7] Number of columns at Y:22255793 does not match the number of samples (169 vs 192)
Error: VCF parse error
After ID fix (step2): 235278
After suspect filter:  235097
After dedup:            235091
Final (no ambig strand):209976
```

Still not right. consul
