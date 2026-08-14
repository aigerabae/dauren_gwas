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
```

Ref alt consistency of each vcf file:
```
mkdir -p ~/biostar/dauren_gwas/redo_july/before_merging
cd ~/biostar/dauren_gwas/redo_july/before_merging

REF=/home/aygera/biostar/dauren_gwas/GRCh38_genome/GRCh38_genome.fa
VCFDIR=~/biostar/dauren_gwas/vcf

# make sure the fasta is indexed (needed by bcftools norm)
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

# summary log
echo -e "sample\ttotal_records\tref_mismatches" > refcheck_summary.tsv

for vcf in ${VCFDIR}/*.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)
    echo "=== Checking $sample ==="

    # -w = warn only, don't modify; just want a diagnostic count per sample
    bcftools norm --check-ref w -f "$REF" "$vcf" -Oz -o /dev/null 2> "${sample}.refcheck.log"

    total=$(bcftools view -H "$vcf" | wc -l)
    mismatches=$(grep -c "REF_MISMATCH" "${sample}.refcheck.log")

    echo -e "${sample}\t${total}\t${mismatches}" >> refcheck_summary.tsv
done

echo ""
echo "=== Summary ==="
column -t refcheck_summary.tsv
```

Checking the results:
```
cd ~/biostar/dauren_gwas/redo_july/before_merging
cat refcheck_summary.tsv

# total mismatches across all samples
awk -F'\t' 'NR>1{sum+=$3} END{print "Total mismatches across all samples:", sum}' refcheck_summary.tsv

# any sample with a notably higher rate than others is worth a closer look —
# could indicate a specific batch/plate issue
awk -F'\t' 'NR>1{print $1, $3/$2*100"%"}' refcheck_summary.tsv | sort -k2 -rn | head
```

No strand errors are present. Proceeding to merging:
```
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
    --make-bed --out d1 --allow-extra-chr --allow-no-sex
echo "d1: $(wc -l < d1.fam) samples, $(wc -l < d1.bim) variants"
```

# PLINK QC:
# 0. Setup and initial conversion
```
echo "d1: $(wc -l < d1.fam) samples, $(wc -l < d1.bim) variants" >> ../qc_log/counts.txt
cat ../qc_log/counts.txt

plink --bfile d1 --list-duplicate-vars ids-only suppress-first --out dup_check --allow-no-sex
wc -l dup_check.dupvar
```

Getting a log of how many variants are retained/lost:
```
mkdir -p ../qc_log
echo "d1: $(wc -l < d1.fam) samples, $(wc -l < d1.bim) variants" >> ../qc_log/counts.txt
cat ../qc_log/counts.txt
```

# 1. Duplicate variant check (sanity check post-VCF-conversion)
```
# if nonzero, remove
if [ -s dup_check.dupvar ]; then
    plink --bfile d1 --exclude dup_check.dupvar --make-bed --out d1b --allow-no-sex
else
    cp d1.bed d1b.bed; cp d1.bim d1b.bim; cp d1.fam d1b.fam
fi
echo "d1b (dedup): $(wc -l < d1b.fam) samples, $(wc -l < d1b.bim) variants" >> ../qc_log/counts.txt
```

# 2. Sample-level QC — call rate, heterozygosity, sex
```
# 2a. Initial call rate filters (liberal, just to stabilize downstream stats)
plink --bfile d1b --geno 0.05 --make-bed --out d2 --allow-no-sex
plink --bfile d2 --mind 0.05 --make-bed --out d3 --allow-no-sex
echo "d3 (geno/mind pre-filter): $(wc -l < d3.fam) samples, $(wc -l < d3.bim) variants" >> ../qc_log/counts.txt

# 2b. Heterozygosity outliers
plink --bfile d3 --het --out het_check --allow-no-sex
```

```python
import pandas as pd
het = pd.read_csv('het_check.het', sep=r'\s+')
het['HET_RATE'] = (het['N(NM)'] - het['O(HOM)']) / het['N(NM)']
mean, sd = het['HET_RATE'].mean(), het['HET_RATE'].std()
outliers = het[(het['HET_RATE'] < mean - 3*sd) | (het['HET_RATE'] > mean + 3*sd)]
outliers[['FID','IID']].to_csv('het_outliers.txt', sep='\t', index=False, header=False)
print(f"Mean het rate: {mean:.4f}, SD: {sd:.4f}, outliers: {len(outliers)}")
```

Mean het rate: 0.2824, SD: 0.0059, outliers: 1

```bash
plink --bfile d3 --remove het_outliers.txt --make-bed --out d4 --allow-no-sex
echo "d4 (het-outlier removed): $(wc -l < d4.fam) samples, $(wc -l < d4.bim) variants" >> ../qc_log/counts.txt
```

2c. Sex check — documented threshold, not eyeballed
```
plink --bfile d4 --impute-sex ycount --make-bed --out d5 --allow-no-sex
```

```python
import pandas as pd
sc = pd.read_csv('d5.sexcheck', sep=r'\s+')
sc.to_csv('sexcheck_full.tsv', sep='\t', index=False)
print(sc['STATUS'].value_counts())
# PROBLEM rows are where imputed sex disagrees with reported/pedigree sex
problems = sc[sc['STATUS'] == 'PROBLEM']
problems[['FID','IID']].to_csv('sex_discordant.txt', sep='\t', index=False, header=False)
print(f"Sex-discordant samples: {len(problems)}")
```

STATUS
PROBLEM    184
Name: count, dtype: int64
Sex-discordant samples: 184

Making histogram:
```
import pandas as pd
import matplotlib.pyplot as plt

sc = pd.read_csv('d5.sexcheck', sep=r'\s+')

plt.figure()
plt.scatter(sc['F'], sc['YCOUNT'] if 'YCOUNT' in sc.columns else sc.index, s=10)
plt.title("F vs. YCOUNT")
plt.xlabel("F (X-heterozygosity)")
plt.ylabel("YCOUNT")
plt.savefig('f_vs_y_d5.png')

plt.figure()
sc['F'].hist(bins=30)
plt.title('F Distribution')
plt.xlabel('F')
plt.ylabel('Frequency')
plt.savefig('f_hist_d5.png')

print(sc[['PEDSEX','SNPSEX','F']].describe())
print(sc['PEDSEX'].value_counts())
print(sc['SNPSEX'].value_counts())

plt.figure(figsize=(10,5))
plt.hist(sc['F'], bins=60)
plt.axvline(0.30, color='green', linestyle='--', label='current female-max')
plt.axvline(0.45, color='red', linestyle='--', label='current male-min')
plt.xlim(0.0, 0.45)   # wider view showing the whole female cluster
plt.legend()
plt.savefig('f_hist_female_side.png')

import matplotlib.pyplot as plt
sc = pd.read_csv('d5.sexcheck', sep=r'\s+')
plt.figure(figsize=(10,5))
plt.hist(sc['F'], bins=60)
plt.axvline(0.30, color='green', linestyle='--', label='current female-max')
plt.axvline(0.45, color='red', linestyle='--', label='current male-min')
plt.xlim(0.25, 0.55)
plt.legend()
plt.savefig('f_hist_zoomed.png')
```

Inspect f_vs_y.png/histograms (as you did earlier) to confirm the F/YCOUNT threshold PLINK used makes sense for your cohort before trusting STATUS; adjust with --sex-threshold-female/--sex-threshold-male if your data's cluster separation warrants a non-default cutoff, and document whichever threshold you land on.

```
plink --bfile d4 --impute-sex 0.31 0.43 --make-bed --out d5 --allow-no-sex
# d6 is skipped because i am not removing any samples from sex disordancy
```

# 3. Relatedness — KING instead of PI_HAT
```
plink2 --bfile d5 --king-cutoff 0.177 --make-bed --out d7
# 0.177 = standard KING cutoff for 2nd-degree-or-closer relatives; adjust if you want a different degree threshold
echo "d7 (KING-pruned): $(wc -l < d7.fam) samples, $(wc -l < d7.bim) variants" >> ../qc_log/counts.txt
```

plink2 --king-cutoff writes d7.king.cutoff.out.id / .in.id listing who was removed/kept — keep these files for your methods reporting (replaces the old manual related_remove.txt).

# 4. Variant-level QC — call rate, HWE, MAF
```
# differential missingness between cases/controls (catches genotyping artifacts)
plink --bfile d7 --test-missing --out diffmiss_check --allow-no-sex
awk '$5 < 1e-5' diffmiss_check.missing > diffmiss_fail.txt
wc -l diffmiss_fail.txt
```
```
plink --bfile d7 --exclude diffmiss_fail.txt --make-bed --out d8 --allow-no-sex

# HWE filter in controls only
plink --bfile d8 --filter-controls --hwe 1e-6 --write-snplist --out hwe_pass_controls
plink --bfile d8 --extract hwe_pass_controls.snplist --make-bed --out d9 --allow-no-sex

echo "d9 (HWE filtered): $(wc -l < d9.fam) samples, $(wc -l < d9.bim) variants" >> ../qc_log/counts.txt

# MAF — real threshold this time, not 0.001 (which only removed monomorphic sites at your N)
plink --bfile d9 --maf 0.01 --make-bed --out d10 --allow-no-sex
echo "d10 (MAF>=0.01): $(wc -l < d10.fam) samples, $(wc -l < d10.bim) variants" >> ../qc_log/counts.txt
```
Adjust --maf up (0.05) if power calculations suggest 0.01 still leaves you underpowered — worth checking sample size vs. detectable MAF explicitly given your N.

# 5. LD pruning + internal PCA + scree plot 
```
plink --bfile d10 --indep-pairwise 50 5 0.2 --out prune --allow-no-sex
plink --bfile d10 --extract prune.prune.in --make-bed --out d10_pruned --allow-no-sex

plink2 --bfile d10_pruned --pca 10 --out pca
```

Inspect pca.eigenvec / a scree plot of pca.eigenval before deciding how many PCs to retain as covariates:
making a scree plot
```python
import pandas as pd
import matplotlib.pyplot as plt
eigenval = pd.read_csv('pca.eigenval', header=None)
plt.figure()
plt.plot(range(1, len(eigenval)+1), eigenval[0], 'o-')
plt.xlabel('PC'); plt.ylabel('Eigenvalue'); plt.title('PCA Scree Plot')
plt.savefig('scree_plot.png')
```

# 6 Saving ancestry outliers according to the internal PCA:
```python
# Flag ancestry/genotyping outliers from PC1/PC2 — not case/control related, general outliers
import pandas as pd
pca = pd.read_csv('pca.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]

mean1, sd1 = pca['PC1'].mean(), pca['PC1'].std()
mean2, sd2 = pca['PC2'].mean(), pca['PC2'].std()

outliers = pca[
    (abs(pca['PC1'] - mean1) > 4*sd1) |
    (abs(pca['PC2'] - mean2) > 4*sd2)
]
print(f"Ancestry/PCA outliers: {len(outliers)}")
outliers[['FID','IID']].to_csv('pca_outliers.txt', sep='\t', index=False, header=False)
print(outliers[['FID','IID','PC1','PC2']])
```

The outliers:
Ancestry/PCA outliers: 6
                     FID                  IID       PC1       PC2
12   206667660001_R08C02  206667660001_R08C02 -0.469707 -0.445571 - yes
57   206767120002_R08C01  206767120002_R08C01  0.050645  0.464647 - no
61   206767120002_R10C01  206767120002_R10C01 -0.444747  0.453026 - yes
131  207859430005_R10C01  207859430005_R10C01 -0.422022  0.045590 - yes
133  207859430005_R11C01  207859430005_R11C01  0.061089  0.464488 - yes
139  207859430006_R03C01  207859430006_R03C01 -0.355511  0.127528 - yes


# 7 PCA with other referent populations to see any outliers:
```
I used worldwide populations with tatars and kazakhs for WGS and infile
cd ref_pops
plink --bfile ../d10 --bmerge merged6.bed merged6.bim merged6.fam --out together
plink --bfile ../d10 --exclude together.missnp --make-bed --out d11
plink --bfile d11 --bmerge merged6.bed merged6.bim merged6.fam --out together
plink --bfile together --geno 0.05 --make-bed --out together2 --allow-no-sex
plink --bfile together2 --mind 0.05 --make-bed --out together3 --allow-no-sex
awk '{print $1, $2, $2, $2}' together3.fam > update_ids.txt
plink2 --bfile together3 --update-ids update_ids.txt --make-bed --out together4

plink2 --bfile together4 --pca 10 --out pca

# metadata:
cat infile.txt | awk '{print $1"\t"$2"\tprevious"}' > m1.txt
cat ../d10.fam | awk '{print $1 "\tcurrent" "\tcurrent" }' > m2.txt
cat m1.txt m2.txt > m3.txt
cat m3.txt | awk '{print $1 "\t" tolower($2) "\t" tolower($3)}' > m4.txt

awk '$2 ~ /^(azeri|buryat|hazara|polish|spanish|turkmen|ukrainian|uygur|tuvan|yakut|han_china|russian_central|altaian|french|german|kazakh|koryak|lezgin|pathan|tatar|current)$/' m4.txt > m5.txt
take from m4 only those rows where column 2 is either: azeri, buryat, hazara, polish, spanish, turkmen, ukrainian, uygur, tuvan, yakut, han_china, russian_central, altaian, french, german, kazakh, koryak, lezgin, pathan, tatar, current

python plot_eigenvec.py pca.eigenvec m5.txt
```

Samples: 
206767120002_R10C01    206767120002_R10C01
206667660001_R08C02    206667660001_R08C02
207859430005_R10C01    207859430005_R10C01
207859430006_R03C01    207859430006_R03C01
206767120003_R10C02    206767120003_R10C02

are deviating from the general kazakh population of the PCA.

I will add all 6 into ref_pops/deviating_samples.txt and delete all 6 shown in the internal PCA since 5/6 are obviously deviating and the other one is also not quite as close to the center of the cluster

```
plink --bfile d10 --remove ref_pops/deviating_samples.txt --make-bed --out d11
```

# 8 - PC1 vs PC2 for case and control 0 
before removing 6 outliers:
```bash
plink2 --bfile d10 --pca 10 --out pca_179
```

```python
import pandas as pd
import matplotlib.pyplot as plt

# check the raw structure first — don't skip this
with open('pca_179.eigenvec') as f:
    print(f.readline())

pca = pd.read_csv('pca.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]   # strip leading '#' from '#FID' if present
print(pca.columns.tolist())
print(pca.head())

fam = pd.read_csv('d10.fam', sep=r'\s+', header=None,
                   names=['FID','IID','PAT','MAT','SEX','PHENO'])

merged = pca.merge(fam[['FID','IID','PHENO']], on=['FID','IID'])
print(merged.shape)                       # sanity check row count — should be close to your sample N
print(merged[['PC1','PC2','PHENO']].describe())   # confirm PC1/PC2 are small floats, PHENO is 1/2

plt.figure(figsize=(7,6))
for pheno, label, color in [(2,'Case','red'), (1,'Control','blue')]:
    sub = merged[merged['PHENO']==pheno]
    plt.scatter(sub['PC1'], sub['PC2'], label=label, alpha=0.6, c=color)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend()
plt.title('PCA colored by case/control status (179)')
plt.savefig('pca_case_control_before.png', dpi=150)
```

After removing 6 outliers:
```bash
plink2 --bfile d11 --pca 10 --out pca_173
```

```python
import pandas as pd
import matplotlib.pyplot as plt

# check the raw structure first — don't skip this
with open('pca.eigenvec') as f:
    print(f.readline())

pca = pd.read_csv('pca_173.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]   # strip leading '#' from '#FID' if present
print(pca.columns.tolist())
print(pca.head())

fam = pd.read_csv('d11.fam', sep=r'\s+', header=None,
                   names=['FID','IID','PAT','MAT','SEX','PHENO'])

merged = pca.merge(fam[['FID','IID','PHENO']], on=['FID','IID'])
print(merged.shape)                       # sanity check row count — should be close to your sample N
print(merged[['PC1','PC2','PHENO']].describe())   # confirm PC1/PC2 are small floats, PHENO is 1/2

plt.figure(figsize=(7,6))
for pheno, label, color in [(2,'Case','red'), (1,'Control','blue')]:
    sub = merged[merged['PHENO']==pheno]
    plt.scatter(sub['PC1'], sub['PC2'], label=label, alpha=0.6, c=color)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend()
plt.title('PCA colored by case/control status (173)')
plt.savefig('pca_case_control_after.png', dpi=150)
```

Checking if any more outliers want removal:
```
import pandas as pd

pca = pd.read_csv('pca_173.eigenvec', sep=r'\s+')  # this is the N=173 recompute
pca.columns = [c.lstrip('#') for c in pca.columns]

mean1, sd1 = pca['PC1'].mean(), pca['PC1'].std()
mean2, sd2 = pca['PC2'].mean(), pca['PC2'].std()
outliers = pca[
    (abs(pca['PC1'] - mean1) > 4*sd1) |
    (abs(pca['PC2'] - mean2) > 4*sd2)
]
print(f"Round 2 outliers: {len(outliers)}")
print(outliers[['FID','IID','PC1','PC2']])
outliers[['FID','IID']].to_csv('pca_outliers_round2.txt', sep='\t', index=False, header=False)
```

Round 2 outliers: 1
                     FID                  IID       PC1       PC2
100  207859430002_R08C01  207859430002_R08C01  0.414859 -0.285921

Another round:
```
echo "207859430002_R08C01    207859430002_R08C01" > remove_round2.txt
plink --bfile d11 --remove remove_round2.txt --make-bed --out d12
plink2 --bfile d12 --pca 10 --out pca_172
```

And checking if more outliers come out:
```python
import pandas as pd
pca = pd.read_csv('pca_172.eigenvec', sep=r'\s+')  # this is the N=173 recompute
pca.columns = [c.lstrip('#') for c in pca.columns]

mean1, sd1 = pca['PC1'].mean(), pca['PC1'].std()
mean2, sd2 = pca['PC2'].mean(), pca['PC2'].std()
outliers = pca[
    (abs(pca['PC1'] - mean1) > 4*sd1) |
    (abs(pca['PC2'] - mean2) > 4*sd2)
]
print(f"Round 3 outliers: {len(outliers)}")
print(outliers[['FID','IID','PC1','PC2']])
outliers[['FID','IID']].to_csv('pca_outliers_round2.txt', sep='\t', index=False, header=False)
```

No more outliers. I will graph the PCA as well:
```
import pandas as pd
import matplotlib.pyplot as plt

pca = pd.read_csv('pca_172.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]   # strip leading '#' from '#FID' if present
print(pca.columns.tolist())
print(pca.head())

fam = pd.read_csv('d11.fam', sep=r'\s+', header=None,
                   names=['FID','IID','PAT','MAT','SEX','PHENO'])

merged = pca.merge(fam[['FID','IID','PHENO']], on=['FID','IID'])
print(merged.shape)                       # sanity check row count — should be close to your sample N
print(merged[['PC1','PC2','PHENO']].describe())   # confirm PC1/PC2 are small floats, PHENO is 1/2

plt.figure(figsize=(7,6))
for pheno, label, color in [(2,'Case','red'), (1,'Control','blue')]:
    sub = merged[merged['PHENO']==pheno]
    plt.scatter(sub['PC1'], sub['PC2'], label=label, alpha=0.6, c=color)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend()
plt.title('PCA colored by case/control status (172)')
plt.savefig('pca_case_control_after2.png', dpi=150)
```

Might want to double check whether the same samples that have sex imputation issues are the same ones that cluster are outliers in pca_case_control

# 9. Association testing — PCA-adjusted logistic regression, not --model
All 10 PCs:
```
plink2 --bfile d12 --glm firth-fallback --covar pca_172.eigenvec --covar-name PC1-PC10 \
    --ci 0.95 --out gwas_firth

# extract clean additive results
grep -w "ADD" gwas_firth.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_sorted.txt
head -20 gwas_firth_sorted.txt
```
    
```results
# (without outliers ethnic): -
# has the lower p value
5	175833369	rs1560036	G	A	Y	A	G	0.380814	N	ADD	172	2.9215	0.261219	1.75088	4.87479	4.10421	4.057e-05	.
4	4367673	rs4330304	G	A	Y	A	G	0.267442	N	ADD	172	3.28557	0.291881	1.85422	5.82185	4.07543	4.59299e-05	.
11	80251941	rs1265425	T	C	Y	C	T	0.25	N	ADD	172	0.302031	0.296069	0.169059	0.539593	-4.04374	5.26062e-05	.
19	18052445	rs62121092	G	A	Y	A	G	0.119186	N	ADD	172	0.187453	0.432414	0.0803196	0.437487	-3.87181	0.000108031	.
10	6178762	rs11257090	A	G	Y	G	A	0.406433	N	ADD	171	0.371281	0.25796	0.223937	0.615572	-3.84088	0.000122592	.
3	161256680	rs4470535	A	G	Y	G	A	0.318713	N	ADD	171	0.332545	0.287584	0.18926	0.584309	-3.82838	0.000128991	.
4	87741379	rs2732204	G	A	Y	A	G	0.415698	N	ADD	172	2.47418	0.237696	1.55276	3.94238	3.81122	0.000138285	.
21	41885483	rs1045071	T	G	Y	G	T	0.252907	N	ADD	172	3.39515	0.321484	1.80805	6.37538	3.8022	0.000143414	.
9	18164547	rs10810914	A	G	Y	G	A	0.40407	N	ADD	172	0.346687	0.279786	0.200348	0.599918	-3.78623	0.000152953	.
3	161296488	rs1478568	T	C	Y	C	T	0.319767	N	ADD	172	0.336694	0.287847	0.191522	0.591904	-3.7818	0.000155697	.
11	2729473	rs170786	C	T	Y	T	C	0.238372	N	ADD	172	3.00231	0.290877	1.6977	5.30946	3.77954	0.000157118	.
11	76661008	rs3740779	G	A	Y	A	G	0.427326	N	ADD	172	0.370141	0.263132	0.220998	0.619934	-3.77708	0.000158678	.
3	161275966	rs6441370	T	C	Y	C	T	0.22619	N	ADD	168	0.276072	0.341244	0.141434	0.538877	-3.77177	0.000162092	.
19	51532449	rs7259990	T	C	Y	C	T	0.424419	N	ADD	172	2.67529	0.261664	1.60193	4.46786	3.76077	0.00016939	.
14	82139260	rs12883722	C	T	Y	T	C	0.125	N	ADD	172	0.203832	0.423409	0.0888928	0.46739	-3.75631	0.000172435	.
15	38901041	rs4923807	T	C	Y	C	T	0.366279	N	ADD	172	2.65547	0.26034	1.59419	4.42327	3.75133	0.0001759	.
13	99399080	rs12586091	T	C	Y	C	T	0.360465	N	ADD	172	2.80447	0.275012	1.63591	4.80774	3.74971	0.000177038	.
10	122120697	rs11200411	C	T	Y	T	C	0.18314	N	ADD	172	3.69231	0.348791	1.86383	7.31458	3.74509	0.00018033	.
3	161348506	rs336577	T	C	Y	C	T	0.215116	N	ADD	172	0.28244	0.340025	0.145043	0.549992	-3.71822	0.000200632	.
3	161354266	rs336570	A	G	Y	G	A	0.215116	N	ADD	172	0.28244	0.340025	0.145043	0.549992	-3.71822	0.000200632	.
```



Need to rerun everything after this stage because now i use 172 samples. P values weren't significant with any number of samples but at least there are clearly no outliers in PCA
Current significant snp list is the final one so far.



Only 1 PC:
```
plink2 --bfile d11 --glm firth-fallback --covar pca.eigenvec --covar-name PC1 \
    --ci 0.95 --out gwas_firth_pc1
grep -w "ADD" gwas_firth_pc1.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_pc1_sorted.txt
head -20 gwas_firth_pc1_sorted.txt
```

"Association results were consistent across sensitivity analyses using either all 10 principal components or PC1 alone as covariates, with 11 of the top 20 loci overlapping between models (Table SX), supporting the robustness of these findings to the choice of population structure adjustment."

I will use PC1-PC10 for the main results.
If you have age or other covariates in phenotypes.tsv, add them: --covar combined_covars.txt --covar-name PC1-PC10,AGE,SEX.


# 10 Compute λ_GC and generate QQ/Manhattan plots on this corrected file
```
# Load once, clean immediately
import pandas as pd
from scipy import stats
import numpy as np

df = pd.read_csv('gwas_firth.PHENO1.glm.logistic.hybrid', sep='\t', low_memory=False)
df = df[df['TEST']=='ADD'].dropna(subset=['P'])
df['#CHROM'] = df['#CHROM'].astype(str).str.strip()

chisq = stats.chi2.isf(df['P'], 1)
lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)
print(f"Lambda GC: {lambda_gc:.4f}")
print(df['#CHROM'].unique())   # confirm clean, 23 unique values, no duplicates

import matplotlib.pyplot as plt

# QQ plot
observed = -np.log10(np.sort(df['P']))
expected = -np.log10(np.linspace(1/len(df), 1, len(df)))

plt.figure(figsize=(6,6))
plt.scatter(expected, observed, s=3, alpha=0.5)
plt.plot([0, max(expected)], [0, max(expected)], 'r--')
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.title(f'QQ Plot (λ_GC = {lambda_gc:.3f})')
plt.savefig('qq_plot_firth.png', dpi=150)

# Manhattan plot
df['CHR_NUM'] = df['#CHROM'].replace({'X': '23', 'Y': '24', 'MT': '25'})
df['CHR_NUM'] = pd.to_numeric(df['CHR_NUM'], errors='coerce')
df = df.dropna(subset=['CHR_NUM'])
df['-logP'] = -np.log10(df['P'])

plt.figure(figsize=(16,5))
colors = ['#1f77b4','#ff7f0e']
x_offset = 0
xticks, xlabels = [], []
for i, chrom in enumerate(sorted(df['CHR_NUM'].unique())):
    sub = df[df['CHR_NUM']==chrom]
    plt.scatter(sub['POS']+x_offset, sub['-logP'], c=colors[i%2], s=4)
    xticks.append(x_offset + sub['POS'].median())
    xlabels.append(str(int(chrom)) if chrom < 23 else {23:'X',24:'Y',25:'MT'}[chrom])
    x_offset += sub['POS'].max()

plt.axhline(-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide (5e-8)')
plt.axhline(-np.log10(1e-5), color='blue', linestyle=':', label='Suggestive (1e-5)')
plt.xticks(xticks, xlabels, rotation=90)
plt.ylabel('-log10(p)')
plt.legend()
plt.title('Manhattan Plot (Firth-corrected, PC1-10 adjusted)')
plt.tight_layout()
plt.savefig('manhattan_plot_firth.png', dpi=150)
```

# 11 calculating effective number of tests
```bash
plink --bfile d10_pruned --indep-pairwise 50 5 0.2 --out prune_eff
wc -l prune_eff.prune.in
```

```python
# Simple LD-pruned-count Bonferroni as a first-pass effective threshold
n_pruned = sum(1 for _ in open('prune_eff.prune.in'))
alpha_naive = 0.05 / n_pruned
print(f"Pruned variant count: {n_pruned}")
print(f"Effective-test-adjusted threshold: {alpha_naive:.3e}")

# Pruned variant count: 50923
# Effective-test-adjusted threshold: 9.819e-07
```

# 12 MAF table sanity check
```
plink --bfile d10 --freq --out d10_freq
awk 'NR==FNR{ids[$3]; next} FNR==1 || ($3 in ids)' \
    <(head -20 gwas_firth_sorted.txt) d10_freq.frq
```

# 13 formal power calculation
```R
install.packages("remotes")
remotes::install_version("genpwr", version = "1.0.4", repos = "https://cloud.r-project.org")

library(genpwr)
# adjust N/Case.Rate to your actual final d10 (or d10_noanc) sample counts
power_result <- genpwr.calc(
  calc = "power", model = "logistic",
  N = 179, Case.Rate = 91/179,
  MAF = seq(0.01, 0.5, 0.01),
  OR = c(1.5, 2, 3, 5),
  Alpha = 5e-8,   # or your effective-test threshold from above
  True.Model = "Additive", Test.Model = "Additive"
)
write.csv(power_result, "power_calculation.csv", row.names = FALSE)
```

Visualizing it:
```python
import pandas as pd
pw = pd.read_csv('power_calculation.csv')
print(pw.columns.tolist())
print(pw.head(10))

import matplotlib.pyplot as plt

power_col = [c for c in pw.columns if c.startswith('Power_at_Alpha_5e-08')][0]  # grab it dynamically

plt.figure(figsize=(8,6))
for or_val in sorted(pw['OR'].unique()):
    sub = pw[pw['OR']==or_val].sort_values('MAF')
    plt.plot(sub['MAF'], sub[power_col], marker='o', label=f'OR={or_val}')

plt.axhline(0.8, color='red', linestyle='--', label='80% power threshold')
plt.xlabel('MAF')
plt.ylabel('Power')
plt.title(f'Power by MAF and OR (N={pw["N_total"].iloc[0]}, α=5e-8)')
plt.legend()
plt.savefig('power_curve.png', dpi=150)
```

"Post-hoc power calculations indicated this study was adequately powered (≥80%) only to detect large-effect common variants (OR≥5, MAF≥0.33) at genome-wide significance (α=5×10⁻⁸); power to detect more modest effects (OR≤3), which are more typical of complex-trait susceptibility loci, remained below 30% across the tested MAF range (Figure SX). This limits our ability to draw firm conclusions from the absence of genome-wide significant findings and supports interpreting the suggestive associations reported here (Table X) as hypothesis-generating, warranting replication in larger, independent cohorts."

Power calculation for suggestive variants:
```R
power_suggestive <- genpwr.calc(
  calc = "power", model = "logistic",
  N = 179, Case.Rate = 91/179,
  MAF = seq(0.01, 0.5, 0.01),
  OR = c(1.5, 2, 3, 5),
  Alpha = 1e-5,   # your suggestive/effective threshold
  True.Model = "Additive", Test.Model = "Additive"
)
write.csv(power_suggestive, "power_calculation_suggestive.csv", row.names = FALSE)
```


Visualizing:
```python
import pandas as pd
pw = pd.read_csv('power_calculation_suggestive.csv')
print(pw.columns.tolist())
print(pw.head(10))

import matplotlib.pyplot as plt

power_col = [c for c in pw.columns if c.startswith('Power')][0]  # grab it dynamically

plt.figure(figsize=(8,6))
for or_val in sorted(pw['OR'].unique()):
    sub = pw[pw['OR']==or_val].sort_values('MAF')
    plt.plot(sub['MAF'], sub[power_col], marker='o', label=f'OR={or_val}')

plt.axhline(0.8, color='red', linestyle='--', label='80% power threshold')
plt.xlabel('MAF')
plt.ylabel('Power')
plt.title(f'Power by MAF and OR (N={pw["N_total"].iloc[0]}, α=1e-05)')
plt.legend()
plt.savefig('power_curve_suggestive.png', dpi=150)
```

"The P values of Bonferroni corrected thresholds for suggestive, 5 and 1% genome-wide significant levels were 1, 0.05 and 0.01, respectively, divided by the number of SNPs used in the GWAS. The suggestive level was first proposed by Lander and Kruglyak [17] and represents the threshold where, under the null hypothesis, one false positive is expected per genome scan."

So the suggestive line is calculated as the -log10( 1 / number of variants) = 5.2552725051

Reference: Guo, Y., Huang, Y., Hou, L., Ma, J., Chen, C., Ai, H., ... & Ren, J. (2017). Genome-wide detection of genetic markers associated with growth and fatness in four pig populations using four approaches. Genetics Selection Evolution, 49(1), 21.

Calculating power:
```
plink --bfile d10 --indep-pairwise 50 5 0.2 --out prune
```
Sum of "leaving" counts across chromosomes = 54,327 pruned-in variants (out of 181,595; matches 181595 - 127268 = 54327).
Bonferroni threshold = 0.05 / 54327 = 9.20 × 10⁻⁷

Method 2:
```
#!/bin/bash
# export_dosages.sh
# Run from ~/biostar/dauren_gwas/redo_july/plink/
mkdir -p simpleM/dosages
cd simpleM/dosages

for chr in $(seq 1 22); do
  plink2 --bfile ../../d10 \
    --chr ${chr} \
    --export A \
    --out chr${chr}_dosage
done

echo "Done. Exported chr1-22, excluded chr23 (X) — handled separately."
```

```R
#!/usr/bin/env Rscript
# simpleM_autosomes.R
# Effective number of independent tests (Li & Ji 2005 / Gao et al. 2008 simpleM),
# windowed for tractability. Autosomes only (chr1-22); chrX excluded and
# handled/reported separately due to het-haploid QC warning on chr23.
# running from ~/biostar/dauren_gwas/redo_july/plink/simpleM

suppressMessages(library(data.table))

WINDOW_SIZE <- 2000
OVERLAP     <- 500
THRESHOLD   <- 0.995
DOSAGE_DIR  <- "dosages"
OUT_LOG     <- "simpleM_per_chr.csv"

simpleM_one_chr <- function(dosage_file, window_size, overlap, threshold) {
  geno <- fread(dosage_file, header = TRUE)
  # plink2 --export A columns: FID IID PAT MAT SEX PHENOTYPE <SNP1_A> <SNP2_A> ...
  meta_cols <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE")
  snp_cols <- setdiff(colnames(geno), meta_cols)
  mat <- as.matrix(geno[, ..snp_cols])
  mode(mat) <- "numeric"

  n_snp <- ncol(mat)
  if (n_snp < 2) return(list(m_eff = n_snp, n_snp = n_snp))

  step <- window_size - overlap
  starts <- seq(1, n_snp, by = step)

  m_eff_total <- 0
  for (s in starts) {
    e <- min(s + window_size - 1, n_snp)
    block <- mat[, s:e, drop = FALSE]

    # drop monomorphic / all-NA columns within this window
    keep <- apply(block, 2, function(x) {
      v <- var(x, na.rm = TRUE)
      !is.na(v) && v > 0
    })
    block <- block[, keep, drop = FALSE]
    k_block <- ncol(block)
    if (k_block < 2) {
      m_eff_total <- m_eff_total + k_block * (step / window_size)
      next
    }

    R <- suppressWarnings(cor(block, use = "pairwise.complete.obs"))
    R[is.na(R)] <- 0
    diag(R) <- 1

    eig <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    eig[eig < 0] <- 0  # guard tiny negative eigenvalues from numerical noise

    cum_var <- cumsum(eig) / sum(eig)
    k <- which(cum_var >= threshold)[1]
    if (is.na(k)) k <- length(eig)

    # scale by non-overlapping fraction of the window to avoid double-counting
    # the boundary SNPs shared with the next window
    weight <- if (e == n_snp) 1 else (step / window_size)
    m_eff_total <- m_eff_total + k * weight
  }

  list(m_eff = m_eff_total, n_snp = n_snp)
}

results <- data.frame(chr = integer(), n_snp = integer(), m_eff = numeric())

for (chr in 1:22) {
  f <- file.path(DOSAGE_DIR, paste0("chr", chr, "_dosage.raw"))
  if (!file.exists(f)) {
    warning(paste("Missing file, skipping:", f))
    next
  }
  cat("Processing chr", chr, "...\n")
  r <- simpleM_one_chr(f, WINDOW_SIZE, OVERLAP, THRESHOLD)
  results <- rbind(results, data.frame(chr = chr, n_snp = r$n_snp, m_eff = r$m_eff))
  cat(sprintf("  chr%d: n_snp=%d, M_eff=%.1f\n", chr, r$n_snp, r$m_eff))
}

fwrite(results, OUT_LOG)

total_snp  <- sum(results$n_snp)
total_meff <- sum(results$m_eff)
threshold_corrected <- 0.05 / total_meff

cat("\n==================== SUMMARY (AUTOSOMES ONLY, chr1-22) ====================\n")
cat(sprintf("Total SNPs (autosomal, post-QC):        %d\n", total_snp))
cat(sprintf("Total effective independent tests:      %.1f\n", total_meff))
cat(sprintf("simpleM-corrected significance threshold: %.3e\n", threshold_corrected))
cat("chrX (23) excluded from this calculation — see README note.\n")
cat("=============================================================================\n")

summary_out <- data.frame(
  metric = c("total_autosomal_snp", "total_M_eff", "corrected_threshold",
             "window_size", "overlap", "cum_var_threshold"),
  value  = c(total_snp, total_meff, threshold_corrected,
             WINDOW_SIZE, OVERLAP, THRESHOLD)
)
fwrite(summary_out, "simpleM_summary.csv")
cat("\nWrote simpleM_per_chr.csv and simpleM_summary.csv\n")
```

simpleM-corrected threshold = 0.05 / 17,699.8 = 2.825 × 10⁻⁶

The genome-wide significance threshold of 5×10⁻⁸ was derived for arrays with genome-wide backbone coverage approximating ~1,000,000 independent tests and is likely overly conservative for a targeted array with 181,595 variants. We therefore additionally computed an array-specific effective number of independent tests using the eigenvalue-based method of Li and Ji (2005), implemented following the windowed approach of Gao et al. (2008) (simpleM), applied to autosomal variants (n = 180,388; chromosome X excluded, see below). Pairwise LD correlation matrices were computed within sliding 2,000-SNP windows (500-SNP overlap), and the number of eigenvalues required to explain 99.5% of cumulative variance was summed across windows and chromosomes, yielding M_eff = 17,699.8 and a corrected significance threshold of 0.05/17,699.8 = 2.83×10⁻⁶. As a complementary, more conservative estimate, we also report a naive Bonferroni correction on the LD-pruned variant set (r² < 0.2; 54,327 variants; threshold = 9.20×10⁻⁷). No variant in this study reached significance under the genome-wide standard, the naive array-wide Bonferroni threshold, the LD-pruned threshold, or the simpleM-corrected threshold; the lowest observed p-value (2.6×10⁻⁵, rs[XXX]) falls short of even the least conservative of these. Chromosome X was excluded from the M_eff calculation due to [an unresolved excess of heterozygous-haploid genotype calls in male samples flagged during QC / pending pseudoautosomal-region splitting]; X-linked association results are reported [descriptively only / with a separate note], consistent with standard practice for excluding sex chromosomes from genome-wide M_eff estimation.
