Versions of everything:
```
bcftools 1.24
plink1.9 and plink2
array analysis cli linux x64 2.1.0
bgzip 1.24
SnpSift 5.4c (build 2026-02-23)
tabix 1.24
```

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

bgzip -c dauren_gwas.vcf > dauren_gwas.vcf.gz
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

Installing needed tools and files for annotation:
```
conda install bioconda::snpsift
bcftools index dauren_gwas.vcf.gz

wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi
```


Split my multiallelics before annotating:
```
bcftools norm -m -any -f /home/aygera/biostar/dauren_gwas/GRCh38_genome/GRCh38_genome.fa \
    dauren_gwas.vcf.gz -Oz -o dauren_gwas.split.vcf.gz
tabix -p vcf dauren_gwas.split.vcf.gz
```

Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:	247814/515/0/85/0/0/0


Annotating (might want to redo it with the entire dbsnp, not just common variants):
```
SnpSift annotate common_all_20180418.vcf.gz dauren_gwas.split.vcf.gz > dauren_annotated.vcf
bgzip dauren_annotated.vcf
tabix -p vcf dauren_annotated.vcf.gz
```

Check annotation rate:
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
# 1260
bcftools view -v snps dauren_annotated.vcf.gz -H | wc -l
# 247089

echo "scale=4; 229/1260" | bc
# .1817
echo "scale=4; 12837/247089" | bc
# .0519
```

For the next step, I consistently had an issue with my drop out rates. Turns out, bcftools concat is the problem. Will use awk instead:
```
# 1. Set ID to RS number where matched (chr_pos will be applied separately below for unmatched)
bcftools annotate --set-id '%INFO/RS' dauren_annotated.vcf.gz -Oz -o step1.vcf.gz
tabix -p vcf step1.vcf.gz

# Adding rs to rsIDs:
bcftools view step1.vcf.gz | \
  awk 'BEGIN{OFS="\t"} /^#/{print; next} {if($3=="."){$3=$1"_"$2} else if($3 !~ /^rs/){$3="rs"$3}; print}' | \
  bgzip > step2.vcf.gz
tabix -p vcf step2.vcf.gz
bgzip -t step2.vcf.gz && echo "bgzip OK"
bcftools view -H step2.vcf.gz | wc -l   # expect 248349

bcftools view -e 'WTD || (INFO/WGT!="." && INFO/WGT>1) || (INFO/SSR!="." && INFO/SSR!=0)' \
    step2.vcf.gz -Oz -o step3.vcf.gz

bcftools view -H step3.vcf.gz | wc -l

bcftools norm -d exact step3.vcf.gz -Oz -o step4.vcf.gz
bcftools view -H step4.vcf.gz | wc -l

bcftools query -f '%ID\n' step4.vcf.gz | sort | uniq -d > dup_ids.txt
wc -l dup_ids.txt

bcftools view -e '(REF="A" & ALT="T") | (REF="T" & ALT="A") | (REF="C" & ALT="G") | (REF="G" & ALT="C")' \
    step4.vcf.gz -Oz -o dauren_final3.vcf.gz
tabix -p vcf dauren_final3.vcf.gz

bgzip -t dauren_final3.vcf.gz && echo "bgzip OK"
bcftools view -H dauren_final3.vcf.gz | wc -l

echo "rs-annotated:  $(bcftools query -f '%ID\n' dauren_final3.vcf.gz | grep -c '^rs')"
echo "chr_pos fallback: $(bcftools query -f '%ID\n' dauren_final3.vcf.gz | grep -vc '^rs')"
```

I need to de duplicate my values again after re-introducing my un-annotated variants:
```
# Pass 1: identify which CHROM:POS:ID combos are duplicated
bcftools query -f '%CHROM\t%POS\t%ID\n' dauren_final3.vcf.gz | sort | uniq -d | awk '{print $1"\t"$2"\t"$3}' > dup_keys.txt
cat dup_keys.txt   # should show all 6: the rs729172 one plus these 5

# Pass 2: rewrite ID for any row matching a dup key, appending ALT
bcftools view dauren_final3.vcf.gz | \
  awk -v dupfile=dup_keys.txt '
    BEGIN{
      OFS="\t"
      while((getline line < dupfile) > 0){
        split(line, a, "\t")
        dup[a[1]"\t"a[2]"\t"a[3]]=1
      }
    }
    /^#/{print; next}
    {
      key=$1"\t"$2"\t"$3
      if(key in dup){ $3=$3"_"$5 }
      print
    }
  ' | bgzip > dauren_final4.vcf.gz

tabix -p vcf dauren_final4.vcf.gz
```

Checking specific numbers - how many were filtered out:
```
# Withdrawn = 0
bcftools view -i 'WTD' step2.vcf.gz -H | wc -l

# Non-unique mapping (WGT > 1) = 0
bcftools view -i 'INFO/WGT>1' step2.vcf.gz -H | wc -l

# Suspect (SSR != 0) = 181
bcftools view -i 'INFO/SSR!="." && INFO/SSR!=0' step2.vcf.gz -H | wc -l

# 7. Final count log
echo "Split+annotated:         248349"
echo "After ID fix (step2):    $(bcftools view -H step2.vcf.gz | wc -l)"
echo "After suspect filter:    $(bcftools view -H step3.vcf.gz | wc -l)"
echo "After dedup:             $(bcftools view -H step4.vcf.gz | wc -l)"
echo "Final (no ambig strand): $(bcftools view -H dauren_final4.vcf.gz | wc -l)"
```

Given how many times this exact -i-scoping mistake has bitten this pipeline, it's worth adopting as a standing rule for the rest of your workflow: never use bcftools annotate -i/-e when you want to edit a subset of records — it always filters the whole output. Use awk (single pass, safe) or the split/tabix-and-recombine-with-concat approach only when the two subsets are truly non-overlapping in genomic coordinates (which yours weren't, hence the earlier concat corruption too).

Checking duplicates:
```
bcftools query -f '%CHROM\t%POS\t%ID\n' dauren_final4.vcf.gz | sort | uniq -d
```

All good. Moving on to PLINK conversion:
```
mkdir plink
cd plink
plink -vcf ../dauren_final4.vcf.gz --double-id --pheno ../phenotypes.tsv \
    --make-bed --out d1 --allow-extra-chr
echo "d1: $(wc -l < d1.fam) samples, $(wc -l < d1.bim) variants"
```

