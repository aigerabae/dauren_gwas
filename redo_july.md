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

# pruning and PCA
plink --bfile together4 --indep-pairwise 50 5 0.2 --out prune --allow-no-sex
plink --bfile together4 --extract prune.prune.in --make-bed --out together4_pruned --allow-no-sex
plink2 --bfile together4_pruned --pca 10 --out pca

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

# pruning for PCA
plink --bfile d11 --indep-pairwise 50 5 0.2 --out prune --allow-no-sex
plink --bfile d11 --extract prune.prune.in --make-bed --out d11_pruned --allow-no-sex
```

# 8 - PC1 vs PC2 for case and control 0 
before removing 6 outliers:
```bash
plink2 --bfile d10_pruned --pca 10 --out pca_179
```

```python
import pandas as pd
import matplotlib.pyplot as plt

# check the raw structure first — don't skip this
with open('pca_179.eigenvec') as f:
    print(f.readline())

pca = pd.read_csv('pca_179.eigenvec', sep=r'\s+')
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
plink2 --bfile d11_pruned --pca 10 --out pca_173
```

```python
import pandas as pd
import matplotlib.pyplot as plt

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

Round 2 outliers: 5
                     FID                  IID       PC1       PC2
85   206767120003_R12C01  206767120003_R12C01 -0.081173 -0.428174
100  207859430002_R08C01  207859430002_R08C01  0.464879  0.373393
119  207859430005_R05C02  207859430005_R05C02  0.420347 -0.384214
149  207859430006_R12C02  207859430006_R12C02 -0.097316 -0.421918
166  207859430008_R09C02  207859430008_R09C02  0.436728  0.123302


Another round:
```
echo "206767120003_R12C01  206767120003_R12C01
207859430002_R08C01  207859430002_R08C01
207859430005_R05C02  207859430005_R05C02
207859430006_R12C02  207859430006_R12C02
207859430008_R09C02  207859430008_R09C02" > remove_round2.txt
plink --bfile d11 --remove remove_round2.txt --make-bed --out d12

# pruning for PCA
plink --bfile d12 --indep-pairwise 50 5 0.2 --out prune --allow-no-sex
plink --bfile d12 --extract prune.prune.in --make-bed --out d12_pruned --allow-no-sex

# pca
plink2 --bfile d12_pruned --pca 10 --out pca_168
```

And checking if more outliers come out:
```python
import pandas as pd
pca = pd.read_csv('pca_168.eigenvec', sep=r'\s+')  # this is the N=173 recompute
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

Round 3 outliers: 4
                     FID                  IID       PC1       PC2
54   206767120002_R07C01  206767120002_R07C01 -0.553443 -0.133591
70   206767120003_R04C01  206767120003_R04C01 -0.319814 -0.066294
115  207859430005_R04C02  207859430005_R04C02 -0.131545 -0.362335
164  207859430008_R11C01  207859430008_R11C01 -0.276053  0.339133

I will remove them again:
```
echo "206767120002_R07C01    206767120002_R07C01
206767120003_R04C01    206767120003_R04C01
207859430005_R04C02    207859430005_R04C02
207859430008_R11C01    207859430008_R11C01" > remove_round3.txt
plink --bfile d12 --remove remove_round3.txt --make-bed --out d13
```

Visualizing after round 2 (before round 3):
```
import pandas as pd
import matplotlib.pyplot as plt

pca = pd.read_csv('pca_168.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]   # strip leading '#' from '#FID' if present
print(pca.columns.tolist())
print(pca.head())

fam = pd.read_csv('d12.fam', sep=r'\s+', header=None,
                   names=['FID','IID','PAT','MAT','SEX','PHENO'])

merged = pca.merge(fam[['FID','IID','PHENO']], on=['FID','IID'])
print(merged.shape)                       # sanity check row count — should be close to your sample N
print(merged[['PC1','PC2','PHENO']].describe())   # confirm PC1/PC2 are small floats, PHENO is 1/2

plt.figure(figsize=(7,6))
for pheno, label, color in [(2,'Case','red'), (1,'Control','blue')]:
    sub = merged[merged['PHENO']==pheno]
    plt.scatter(sub['PC1'], sub['PC2'], label=label, alpha=0.6, c=color)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend()
plt.title('PCA colored by case/control status (168)')
plt.savefig('pca_case_control_after2.png', dpi=150)
```

Running PCA on round3:
```
# pruning for PCA
plink --bfile d13 --indep-pairwise 50 5 0.2 --out prune --allow-no-sex
plink --bfile d13 --extract prune.prune.in --make-bed --out d13_pruned --allow-no-sex

# pca
plink2 --bfile d13_pruned --pca 10 --out pca_164
```

Calculating if any more outliers are left after round 3:
```
import pandas as pd
pca = pd.read_csv('pca_164.eigenvec', sep=r'\s+')  # this is the N=173 recompute
pca.columns = [c.lstrip('#') for c in pca.columns]

mean1, sd1 = pca['PC1'].mean(), pca['PC1'].std()
mean2, sd2 = pca['PC2'].mean(), pca['PC2'].std()
outliers = pca[
    (abs(pca['PC1'] - mean1) > 4*sd1) |
    (abs(pca['PC2'] - mean2) > 4*sd2)
]
print(f"Round 4 outliers: {len(outliers)}")
print(outliers[['FID','IID','PC1','PC2']])
outliers[['FID','IID']].to_csv('pca_outliers_round3.txt', sep='\t', index=False, header=False)
```

No more outliers. Visualizing:
```
import pandas as pd
import matplotlib.pyplot as plt

pca = pd.read_csv('pca_164.eigenvec', sep=r'\s+')
pca.columns = [c.lstrip('#') for c in pca.columns]   # strip leading '#' from '#FID' if present
print(pca.columns.tolist())
print(pca.head())

fam = pd.read_csv('d13.fam', sep=r'\s+', header=None,
                   names=['FID','IID','PAT','MAT','SEX','PHENO'])

merged = pca.merge(fam[['FID','IID','PHENO']], on=['FID','IID'])
print(merged.shape)                       # sanity check row count — should be close to your sample N
print(merged[['PC1','PC2','PHENO']].describe())   # confirm PC1/PC2 are small floats, PHENO is 1/2

plt.figure(figsize=(7,6))
for pheno, label, color in [(2,'Case','red'), (1,'Control','blue')]:
    sub = merged[merged['PHENO']==pheno]
    plt.scatter(sub['PC1'], sub['PC2'], label=label, alpha=0.6, c=color)
plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend()
plt.title('PCA colored by case/control status (168)')
plt.savefig('pca_case_control_after3.png', dpi=150)
```

Might want to double check whether the same samples that have sex imputation issues are the same ones that cluster are outliers in pca_case_control

re-making scree plot with 164 samples:
```
Inspect pca.eigenvec / a scree plot of pca.eigenval before deciding how many PCs to retain as covariates:
making a scree plot
```python
import pandas as pd
import matplotlib.pyplot as plt
eigenval = pd.read_csv('pca_164.eigenval', header=None)
plt.figure()
plt.plot(range(1, len(eigenval)+1), eigenval[0], 'o-')
plt.xlabel('PC'); plt.ylabel('Eigenvalue'); plt.title('PCA Scree Plot')
plt.savefig('scree_plot2.png')
```

Following outlier removal, the PCA eigenvalue spectrum showed no discrete elbow (Supplementary Figure X), consistent with a homogeneous residual population; association analyses were therefore performed adjusting for [N] PCs, with results robust to inclusion of up to 10 PCs (Supplementary Table Y)."

# 9. Association testing — PCA-adjusted logistic regression, not --model
All 10 PCs:
```
plink2 --bfile d13 --glm firth-fallback --covar pca_164.eigenvec --covar-name PC1-PC10 \
    --ci 0.95 --out gwas_firth

# extract clean additive results
grep -w "ADD" gwas_firth.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_sorted.txt
head -20 gwas_firth_sorted.txt
```
    
```results
5	175833369	rs1560036	G	A	Y	A	G	0.384146	N	ADD	164	3.35783	0.277933	1.94752	5.78943	4.358221.31122e-05	.
3	95829911	rs13081814	A	G	Y	G	A	0.219512	N	ADD	164	4.30402	0.346498	2.1824	8.48817	4.212292.52794e-05	.
X	2623614	rs3795179	A	G	Y	G	A	0.20122NADD	164	0.210716	0.382377	0.0995906	0.445839	-4.07253	4.65057e-05	.
11	80251941	rs1265425	T	C	Y	C	T	0.25	N	ADD	164	0.273695	0.318485	0.146613	0.510931	-4.06845	4.73273e-05	.
4	4367673	rs4330304	G	A	Y	A	G	0.27439NADD	164	3.64393	0.320846	1.94297	6.83401	4.03016	5.57379e-05	.
14	102101959	rs10873531	A	G	Y	G	A	0.198171	N	ADD	164	0.256903	0.339516	0.13206	0.499765-4.00292	6.25657e-05	.
19	18052445	rs62121092	G	A	Y	A	G	0.121951	N	ADD	164	0.174404	0.44975	0.0722317	0.421098-3.88301	0.000103171	.
15	38901041	rs4923807	T	C	Y	C	T	0.375	N	ADD	164	2.93181	0.285189	1.67642	5.1273	3.771610.000162199	.
9	86824392	rs4878016	G	A	Y	A	G	0.289634	N	ADD	164	0.305911	0.314244	0.165238	0.566343	-3.76924	0.000163744	.
10	93022978	rs835259	G	T	Y	T	G	0.16358	N	ADD	162	4.35175	0.39721	1.99783	9.47915	3.70227	0.000213681	.
3	190574239	rs3773971	A	G	Y	G	A	0.204268	N	ADD	164	3.99861	0.374876	1.91785	8.33688	3.697080.000218092	.
14	102099124	rs7145597	G	A	Y	A	G	0.131098	N	ADD	164	0.206063	0.42752	0.0891445	0.476328-3.69473	0.000220119	.
15	57316925	rs11071311	C	T	Y	T	C	0.439024	N	ADD	164	0.36001	0.278481	0.208579	0.621381-3.66855	0.000243928	.
2	136105973	2_136105973	C	A	Y	A	C	0.0864198	N	ADD	162	7.71144	0.557101	2.5878	22.97953.66667	0.000245728	.
3	70963429	rs3773494	A	G	Y	G	A	0.310976	N	ADD	164	3.00408	0.301434	1.66391	5.42367	3.649130.000263132	.
3	161275966	rs6441370	T	C	Y	C	T	0.21875	N	ADD	160	0.275161	0.354826	0.137265	0.551588	-3.63671	0.000276138	.
2	118692587	rs512681	C	A	Y	A	C	0.134969	N	ADD	163	4.66735	0.423727	2.0342	10.709	3.635810.000277105	.
14	95858869	rs1885152	C	T	Y	T	C	0.310976	N	ADD	164	2.80459	0.28452	1.60578	4.89839	3.62455	0.00028946	.
2	129887241	rs4662668	A	G	Y	G	A	0.307927	N	ADD	164	2.94994	0.300613	1.63656	5.31736	3.598590.000319945	.
13	99399080	rs12586091	T	C	Y	C	T	0.362805	N	ADD	164	2.68696	0.275131	1.567	4.60737	3.592510.000327514	.
```

Only 1 PC:
```
plink2 --bfile d13 --glm firth-fallback --covar pca_164.eigenvec --covar-name PC1 \
    --ci 0.95 --out gwas_firth_pc1
grep -w "ADD" gwas_firth_pc1.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_pc1_sorted.txt
head -20 gwas_firth_pc1_sorted.txt
```

```results
5	175833369	rs1560036	G	A	Y	A	G	0.384146	N	ADD	164	2.96919	0.258143	1.79022	4.92459	4.215842.48845e-05	.
3	95829911	rs13081814	A	G	Y	G	A	0.219512	N	ADD	164	3.32509	0.311544	1.80558	6.12337	3.856590.000114981	.
19	18052445	rs62121092	G	A	Y	A	G	0.121951	N	ADD	164	0.203762	0.419376	0.0895671	0.463549	-3.79327	0.000148679	.
15	38901041	rs4923807	T	C	Y	C	T	0.375	N	ADD	164	2.6218	0.254792	1.59119	4.31996	3.782940.000154988	.
10	15029619	rs11818063	C	T	Y	T	C	0.35061	N	ADD	164	2.7835	0.272245	1.63251	4.74599	3.760250.000169742	.
11	76661008	rs3740779	G	A	Y	A	G	0.426829	N	ADD	164	0.3888	0.254252	0.236215	0.639951-3.71556	0.000202756	.
3	161275966	rs6441370	T	C	Y	C	T	0.21875	N	ADD	160	0.297926	0.328902	0.156368	0.567637	-3.68167	0.000231707	.
4	4367673	rs4330304	G	A	Y	A	G	0.27439NADD	164	2.73778	0.273682	1.60118	4.68119	3.67999	0.000233243	.
22	37129691	rs84460	A	G	Y	G	A	0.359756N	ADD	164	2.52254	0.252571	1.53762	4.13835	3.6634	0.00024889	.
X	2623614	rs3795179	A	G	Y	G	A	0.20122NADD	164	0.302967	0.32616	0.15987	0.574148	-3.66118	0.000251056	.
14	95858869	rs1885152	C	T	Y	T	C	0.310976	N	ADD	164	2.64824	0.266099	1.572	4.46129	3.659890.000252325	.
11	80251941	rs1265425	T	C	Y	C	T	0.25	N	ADD	164	0.354534	0.283351	0.203455	0.617798	-3.6596	0.000252608	.
10	93022978	rs835259	G	T	Y	T	G	0.16358	N	ADD	162	3.66841	0.358321	1.8175	7.40425	3.627360.000286329	.
14	102101959	rs10873531	A	G	Y	G	A	0.198171	N	ADD	164	0.333806	0.305711	0.183347	0.607737	-3.589	0.00033195	.
10	122120697	rs11200411	C	T	Y	T	C	0.185976	N	ADD	164	3.24122	0.328348	1.70301	6.16877	3.581410.000341739	.
16	8402882	rs12597409	A	G	Y	G	A	0.265244N	ADD	164	0.370747	0.277858	0.215062	0.639133-3.57101	0.000355602	.
3	161256680	rs4470535	A	G	Y	G	A	0.315951	N	ADD	163	0.376216	0.274727	0.219578	0.644593	-3.55841	0.000373107	.
3	54271643	rs7615792	G	A	Y	A	G	0.237805	N	ADD	164	0.35973	0.289553	0.203943	0.634518-3.53097	0.000414034	.
11	61000572	rs2905517	A	G	Y	G	A	0.243902	N	ADD	164	2.79362	0.291386	1.57812	4.94534	3.525690.000422375	.
3	161296488	rs1478568	T	C	Y	C	T	0.317073	N	ADD	164	0.3818	0.274546	0.222916	0.653927-3.5071	0.000453017	.
```

Ancestry outliers were identified by iterative PCA: at each round, PCA was performed on an LD-pruned variant set (r²<0.2, 50-SNP window, 5-SNP step), and samples deviating more than 4 SD from the cohort mean on PC1 or PC2 were flagged and removed. This process was repeated, re-pruning and recomputing PCA at each round, until no further outliers were detected (three rounds; N=179→173→168→164). The resulting eigenvalue spectrum showed no discrete elbow (Supplementary Figure X), consistent with a genetically homogeneous residual cohort. 

Calculating for PC1-3 + sex adjusted:
```
plink2 --bfile d13 --glm firth-fallback sex --covar pca_164.eigenvec \
    --covar-name PC1-PC3 --ci 0.95 --out gwas_firth_pc1-3

grep -w "ADD" gwas_firth_pc1-3.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_pc1-3_sorted.txt
head -100 gwas_firth_pc1-3_sorted.txt
head -39 gwas_firth_pc1-3_sorted.txt > top39.txt
```

Top 39 hits because they are suggestive under 5×10⁻⁴ threshold
```results:
5	175833369	rs1560036	G	A	Y	A	G	0.384146	N	ADD	164	2.97915	0.258558	1.79477	4.94512	4.22203	2.42112e-05	.
3	95829911	rs13081814	A	G	Y	G	A	0.219512	N	ADD	164	3.40698	0.314977	1.83764	6.31654	3.8918	9.95023e-05	.
15	38901041	rs4923807	T	C	Y	C	T	0.375	N	ADD	164	2.7146	0.261194	1.62696	4.52934	3.82339	0.000131632	.
10	15029619	rs11818063	C	T	Y	T	C	0.35061	N	ADD	164	2.86019	0.275403	1.66714	4.90703	3.81582	0.00013573	.
19	18052445	rs62121092	G	A	Y	A	G	0.121951	N	ADD	164	0.199424	0.423012	0.0870382	0.456927	-3.81153	0.000138111	.
3	161275966	rs6441370	T	C	Y	C	T	0.21875	N	ADD	160	0.277438	0.336564	0.143444	0.5366	-3.80955	0.000139221	.
11	76661008	rs3740779	G	A	Y	A	G	0.426829	N	ADD	164	0.379303	0.257935	0.228787	0.628842	-3.75838	0.000171016	.
10	93022978	rs835259	G	T	Y	T	G	0.16358	N	ADD	162	3.94683	0.368882	1.91538	8.13281	3.72182	0.000197793	.
3	161256680	rs4470535	A	G	Y	G	A	0.315951	N	ADD	163	0.34899	0.283275	0.200303	0.608047	-3.71622	0.000202228	.
14	95858869	rs1885152	C	T	Y	T	C	0.310976	N	ADD	164	2.6996	0.268288	1.59563	4.56738	3.70164	0.000214209	.
11	80251941	rs1265425	T	C	Y	C	T	0.25	N	ADD	164	0.341822	0.291724	0.192968	0.605503	-3.67973	0.000233482	.
4	4367673	rs4330304	G	A	Y	A	G	0.27439	N	ADD	164	2.744530.274811	1.60158	4.70314	3.67383	0.000238942	.
X	2623614	rs3795179	A	G	Y	G	A	0.20122	N	ADD	164	0.2992460.329751	0.156799	0.571101	-3.65879	0.00025341	.
20	11812314	rs6134366	G	A	Y	A	G	0.141104	N	ADD	163	0.247386	0.383952	0.116562	0.525045	-3.63797	0.000274799	.
3	161296488	rs1478568	T	C	Y	C	T	0.317073	N	ADD	164	0.359195	0.281559	0.206855	0.623726	-3.6365	0.000276369	.
3	161232242	rs4518145	C	T	Y	T	C	0.22561	N	ADD	164	0.306937	0.326271	0.16193	0.581798	-3.62004	0.000294562	.
22	37129691	rs84460	A	G	Y	G	A	0.359756	N	ADD	164	2.54843	0.258765	1.53466	4.23188	3.61516	0.000300158	.
16	8402882	rs12597409	A	G	Y	G	A	0.265244	N	ADD	164	0.364697	0.279682	0.210798	0.630954	-3.60656	0.000310289	.
3	161348506	rs336577	T	C	Y	C	T	0.207317	N	ADD	164	0.303097	0.331553	0.158257	0.580497	-3.60034	0.000317803	.
3	161354266	rs336570	A	G	Y	G	A	0.207317	N	ADD	164	0.303097	0.331553	0.158257	0.580497	-3.60034	0.000317803	.
3	161354558	rs336569	G	T	Y	T	G	0.207317	N	ADD	164	0.303097	0.331553	0.158257	0.580497	-3.60034	0.000317803	.
3	161226077	rs4370045	C	T	Y	T	C	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161232526	rs4273380	G	A	Y	A	G	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161235821	rs6775627	G	A	Y	A	G	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161245322	rs4597724	C	T	Y	T	C	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161297801	rs3849524	T	C	Y	C	T	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161308929	rs6441373	G	T	Y	T	G	0.207317	N	ADD	164	0.306601	0.329492	0.160735	0.58484	-3.58798	0.000333253	.
3	161298633	rs7647327	C	T	Y	T	C	0.20679	N	ADD	162	0.304756	0.331433	0.159161	0.583537	-3.58517	0.000336857	.
14	102101959	rs10873531	A	G	Y	G	A	0.198171	N	ADD	164	0.334048	0.309004	0.182299	0.612116	-3.5484	0.000387581	.
10	122120697	rs11200411	C	T	Y	T	C	0.185976	N	ADD	164	3.25604	0.333122	1.69486	6.25523	3.54378	0.000394441	.
5	170054447	rs17072023	C	T	Y	T	C	0.22561	N	ADD	164	0.342202	0.303252	0.188866	0.620029	-3.53617	0.000405969	.
2	118692587	rs512681	C	A	Y	A	C	0.134969	N	ADD	163	4.1236	0.401296	1.87799	9.05441	3.53038	0.000414968	.
9	89661201	rs7848364	C	T	Y	T	C	0.137195	N	ADD	164	3.99335	0.393685	1.846	8.63857	3.5171	0.000436284	.
11	61000572	rs2905517	A	G	Y	G	A	0.243902	N	ADD	164	2.8222	0.294998	1.58301	5.03141	3.51702	0.000436413	.
6	31103868	rs9263600	C	T	Y	T	C	0.146341	N	ADD	164	3.78526	0.378849	1.80143	7.95377	3.51357	0.000442127	.
15	57316925	rs11071311	C	T	Y	T	C	0.439024	N	ADD	164	0.408742	0.255521	0.247714	0.674449	-3.50136	0.000462882	.
13	99399080	rs12586091	T	C	Y	C	T	0.362805	N	ADD	164	2.45007	0.256749	1.48127	4.0525	3.49024	0.000482586	.
19	51532449	rs7259990	T	C	Y	C	T	0.435976	N	ADD	164	2.39044	0.24971	1.46529	3.8997	3.48996	0.000483085	.
3	182217390	rs9853234	T	C	Y	C	T	0.317073	N	ADD	164	2.50539	0.263613	1.49447	4.20014	3.48407	0.000493845	.
```

I will stick with PC-PC3 because it has a bir more variation captured than PC1 but holds more degrees of freedom than PC1-PC10

# 10 Compute λ_GC and generate QQ/Manhattan plots on this corrected file
```
# Load once, clean immediately
import pandas as pd
from scipy import stats
import numpy as np

df = pd.read_csv('gwas_firth_pc1-3.PHENO1.glm.logistic.hybrid', sep='\t', low_memory=False)
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
plt.axhline(-np.log10(5e-4), color='yellow', linestyle='-.', label='Weak suggestive (5e-4)')
plt.xticks(xticks, xlabels, rotation=90)
plt.ylabel('-log10(p)')
plt.legend()
plt.title('Manhattan Plot (Firth-corrected, PC1-3 + Sex adjusted)')
plt.tight_layout()
plt.savefig('manhattan_plot_firth.png', dpi=150)
```

Lambda GC: 1.0619


# 11 calculating effective number of tests
```bash
plink --bfile d13_pruned --indep-pairwise 50 5 0.2 --out prune_eff
wc -l prune_eff.prune.in
```

```python
# Simple LD-pruned-count Bonferroni as a first-pass effective threshold
n_pruned = sum(1 for _ in open('prune_eff.prune.in'))
alpha_naive = 0.05 / n_pruned
print(f"Pruned variant count: {n_pruned}")
print(f"Effective-test-adjusted threshold: {alpha_naive:.3e}")

# Pruned variant count: 50863
# Effective-test-adjusted threshold: 9.830e-07
```

# 12 MAF table sanity check
```
# for all
plink --bfile d13 --freq --out d13_freq
awk 'NR==FNR{ids[$3]; next} FNR==1 || ($3 in ids)' \
    <(head -20 gwas_firth_pc1-3_sorted.txt) d13_freq.frq

# for cases
awk '$6 == 2 {print $1, $2}' d13.fam > cases.txt
plink --bfile d13 --keep cases.txt --make-bed --out d13_case
plink --bfile d13_case --freq --out d13_case

# for controls
awk '$6 == 1 {print $1, $2}' d13.fam > controls.txt
plink --bfile d13 --keep controls.txt --make-bed --out d13_control
plink --bfile d13_control --freq --out d13_control
```

# 13 formal power calculation
```R
install.packages("remotes")
remotes::install_version("genpwr", version = "1.0.4", repos = "https://cloud.r-project.org")

library(genpwr)
# adjust N/Case.Rate to your actual final d13 (or d13_noanc) sample counts
power_result <- genpwr.calc(
  calc = "power", model = "logistic",
  N = 164, Case.Rate = 84/164,
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

Post-hoc power calculations at the final analytic sample size (N=164) indicated this study was adequately powered (≥80%) only to detect large-effect common variants (OR≥5, MAF≥0.45) at genome-wide significance (α=5×10⁻⁸); power to detect more modest effects (OR≤3) remained below 21% across the tested MAF range, and power for OR≤2 was negligible (<2%) regardless of allele frequency (Figure SX). This limits our ability to draw firm conclusions from the absence of genome-wide significant findings and supports interpreting the suggestive associations reported here (Table X) as hypothesis-generating, warranting replication in larger, independent cohorts.


Power calculation for suggestive variants:
```R
power_suggestive <- genpwr.calc(
  calc = "power", model = "logistic",
  N = 164, Case.Rate = 84/179,
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

At the conventional suggestive-significance threshold (α=1×10⁻⁵), power improved modestly but remained limited for anything short of large effects: OR=5 crossed 80% power at MAF≥0.20, OR=3 reached a maximum of 58% power even at MAF=0.5, and power for OR=2 remained below 9% and for OR=1.5 below 1% across the entire tested MAF range (Figure SX). This confirms that even under a substantially relaxed significance threshold, the study remains underpowered to detect small-to-moderate effect sizes typical of complex-trait susceptibility loci.

A suggestion:
```
"The P values of Bonferroni corrected thresholds for suggestive, 5 and 1% genome-wide significant levels were 1, 0.05 and 0.01, respectively, divided by the number of SNPs used in the GWAS. 
The suggestive level was first proposed by Lander and Kruglyak [17] and represents the threshold where, under the null hypothesis, one false positive is expected per genome scan."
So the suggestive line is calculated as the -log10( 1 / number of variants) = 5.2552725051

Reference: Guo, Y., Huang, Y., Hou, L., Ma, J., Chen, C., Ai, H., ... & Ren, J. (2017). Genome-wide detection of genetic markers associated with growth and fatness in four pig populations using four approaches. Genetics Selection Evolution, 49(1), 21.
```

Calculating power:
```
plink --bfile d13 --indep-pairwise 50 5 0.2 --out prune
```
Sum of "leaving" counts across chromosomes = 54,350 pruned-in variants (out of 181,595; matches 181595 - 127245 = 54350).
Bonferroni threshold = 0.05 / 54350 = 9.20 × 10⁻⁷

Method 2:
```
#!/bin/bash
# export_dosages.sh
# Run from ~/biostar/dauren_gwas/redo_july/plink/
mkdir -p simpleM/dosages
cd simpleM/dosages

for chr in $(seq 1 22); do
  plink2 --bfile ../../d13 \
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

The genome-wide significance threshold of 5×10⁻⁸ was derived for arrays with genome-wide backbone coverage approximating ~1,000,000 independent tests and is likely overly conservative for a targeted array with 181,595 variants. We therefore additionally computed an array-specific effective number of independent tests using the eigenvalue-based method of Li and Ji (2005), implemented following the windowed approach of Gao et al. (2008) (simpleM), applied to autosomal variants (n = 180,388; chromosome X excluded, see below). Pairwise LD correlation matrices were computed within sliding 2,000-SNP windows (500-SNP overlap), and the number of eigenvalues required to explain 99.5% of cumulative variance was summed across windows and chromosomes, yielding M_eff = 17,699.8 and a corrected significance threshold of 0.05/17,699.8 = 2.83×10⁻⁶. As a complementary, more conservative estimate, we also report a naive Bonferroni correction on the LD-pruned variant set (r² < 0.2; 54,327 variants; threshold = 9.20×10⁻⁷). No variant in this study reached significance under the genome-wide standard, the naive array-wide Bonferroni threshold, the LD-pruned threshold, or the simpleM-corrected threshold; the lowest observed p-value (2.6×10⁻⁵, rs[XXX]) falls short of even the least conservative of these. Chromosome X was excluded from the M_eff calculation due to [an unresolved excess of heterozygous-haploid genotype calls in male samples flagged during QC / pending pseudoautosomal-region splitting]; X-linked association results are reported [descriptively only / with a separate note], consistent with standard practice for excluding sex chromosomes from genome-wide M_eff estimation.

The effective number of independent tests, estimated via windowed eigenvalue decomposition (Li & Ji, 2005; Gao et al., 2008) on the final analytic sample (N=164, autosomal variants, n=180,388), was M_eff=16,319.2, yielding a study-specific significance threshold of 0.05/16,319.2 = 3.06×10⁻⁶. The top association (rs1560036, p=2.42×10⁻⁵) did not reach this threshold.

```
Result:
==================== SUMMARY (AUTOSOMES ONLY, chr1-22) ====================
Total SNPs (autosomal, post-QC):        180388
Total effective independent tests:      16319.2
simpleM-corrected significance threshold: 3.064e-06
chrX (23) excluded from this calculation — see README note.
=============================================================================
```

# 14: I will try other models:
1. Sex not adjusted:
```
plink2 --bfile d13 --glm firth-fallback --covar pca_164.eigenvec --covar-name PC1-PC3 \
    --ci 0.95 --out gwas_firth_pc1-3_nosex
grep -w "ADD" gwas_firth_pc1-3_nosex.PHENO1.glm.logistic.hybrid | awk '$18!="NA"' | sort -gk 18,18 > gwas_firth_pc1-3_nosex_sorted.txt
head -20 gwas_firth_pc1-3_nosex_sorted.txt
```

Since age at recruitment was unavailable for this cohort, the primary association model adjusted for genetically-imputed sex (derived from X-chromosome heterozygosity; see Methods) in addition to the first three principal components. Adjustment for sex had minimal impact on association estimates relative to the PC-only model: the top-associated variant (rs1560036) showed nearly identical effect size and significance with sex included (OR=2.98, p=2.42×10⁻⁵) compared to without (OR=2.98, p=2.42×10⁻⁵), and the full top-20 list of associated loci was unchanged in composition and rank order between the two models. This stability suggests that residual sex-related confounding, if present, does not meaningfully influence the reported associations, and supports the inclusion of sex as a conservative but low-impact covariate in the final model.

```

2. Full Firth vs. Firth-fallback (quick, worth doing)
```
plink2 --bfile d13 --glm firth sex --covar pca_164.eigenvec \
    --covar-name PC1-PC3 --ci 0.95 --out gwas_full_firth
```

Checking result:
```
awk -F'\t' '$10=="ADD" && $17!="NA"' gwas_full_firth.PHENO1.glm.firth | sort -gk17,17 | head -20
```
As a robustness check on the primary Firth-fallback logistic regression, we additionally fit full Firth regression (bias-corrected penalized likelihood applied uniformly across all variants) under the same covariate specification (genetically-imputed sex, PC1–3). Effect estimates and association rankings were highly concordant between the two approaches: the top-associated variant (rs1560036) showed OR=2.98 (p=2.42×10⁻⁵) under Firth-fallback and OR=2.84 (p=3.25×10⁻⁵) under full Firth, with the full top-20 loci list identical in composition and nearly identical in rank order between models. As expected, full Firth produced systematically more conservative estimates (effect sizes shrunk slightly toward the null), consistent with its uniform application of bias correction rather than selective correction at variants showing evidence of separation. This concordance indicates the reported associations are not artifacts of the fallback regression strategy.

3. Permutation
```
plink --bfile d13 --model --mperm 10000 --out gwas_perm
```

or for genuinely adaptive (variable permutation count, stops early for clearly-null variants, spends more permutations on promising ones):
```
plink --bfile d13 --assoc --perm --out gwas_perm
```

Checking result:
```
sort -gk4,4 gwas_perm.model.best.mperm | head -20        # the first
sort -gk3,3 gwas_perm.assoc.perm | head -20              # the second
```

The pointwise empirical p-value for rs1560036 from adaptive permutation (1,000,000 permutations) was 1.1×10⁻⁵, closely matching the Firth-regression analytical p-value (2.42×10⁻⁵–3.25×10⁻⁵ across model specifications), confirming the analytical approximation was reliable at this sample size. However, the genome-wide max(T) permutation-corrected p-value (accounting for simultaneous testing of all 181,595 variants) was 0.656, indicating this association does not exceed chance expectation once multiple testing across the full variant set is properly accounted for.

4. Genotyping model
```
plink2 --bfile d13 --glm firth-fallback sex genotypic --covar pca_164.eigenvec \
    --covar-name PC1-PC3 --ci 0.95 --out gwas_firth_genotypic
```

Checking result:
```
awk -F'\t' '$11=="ADD" && $18!="NA"' gwas_firth_genotypic.PHENO1.glm.logistic.hybrid | sort -gk18,18 | head -20
awk -F'\t' '$11=="DOMDEV"' gwas_firth_genotypic.PHENO1.glm.logistic.hybrid | sort -gk18,18 | head -20
```

As an additional sensitivity analysis, we fit a genotypic model allowing for non-additive (dominance-deviation) effects at each variant, using the same covariate specification (genetically-imputed sex, PC1–3). The additive-term (ADD) results from this model showed reduced statistical evidence relative to the primary additive-only model: the top hit remained rs1560036, but its significance attenuated (OR=3.00, p=1.28×10⁻⁴ under the genotypic model vs. p=2.42×10⁻⁵ under the additive-only model), and several loci from the primary top-20 list (e.g., rs13081814, rs4923807, rs62121092, rs11818063) fell outside the top ranks entirely. This reduction in power is expected: the genotypic model estimates an additional parameter (dominance deviation) per variant, which increases the effective number of parameters relative to the available sample size (N=164). Dominance-deviation terms could not be reliably estimated for a substantial number of low-frequency variants (minor allele frequency approximately 0.01–0.10 in the affected set), returning undefined estimates due to near-perfect collinearity between the additive and dominance terms (CORR_TOO_HIGH) — a well-recognized limitation of fitting non-additive genetic models at low allele frequencies in modestly-sized cohorts. Given these constraints, the genotypic model was not adopted for the primary analysis; the additive model, which is both statistically better-supported at this sample size and the standard convention for GWAS discovery analyses, was retained as primary. This sensitivity analysis does not provide evidence for dominance effects at the top loci and instead illustrates that any deviation from additivity in this dataset cannot be reliably distinguished from noise given current sample size.

15. Annotation
```
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' top39.txt > top39_for_annotation.tsv
conda activate bcftols_env
pip install pandas requests tqdm
```

I ran the rest in ipynb annotating but i also added eqtl data via script
```
python add_gtex_eqtl.py --input top39_FINAL_annotation.tsv --output your_annotated_table_with_gtex.tsv 
```

# 16. PRS
```
# getting bed file for liftover
awk 'BEGIN{OFS="\t"} {print "chr"$1, $4-1, $4, $2}' d13.bim > d13_liftover.bed
# i downloaded liftovered file as my_d13_to37.bed and 1284 fails as liftover_fail.txt; 1206 of them are rsids and 78 aren't
```

for those 1206 i did this:
```python
import pandas as pd
import requests
import time

def lookup_rsids_grch37(rsids):
    url = "https://grch37.rest.ensembl.org/variation/human"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    results = []
    chunk_size = 200  # Ensembl POST limit
    for i in range(0, len(rsids), chunk_size):
        chunk = rsids[i:i+chunk_size]
        payload = {"ids": chunk}
        
        response = requests.post(url, headers=headers, json=payload)
        if response.status_code == 200:
            data = response.json()
            for rsid, info in data.items():
                mappings = info.get('mappings', [])
                for mapping in mappings:
                    results.append({
                        "rsID": rsid,
                        "Chromosome": mapping.get('seq_region_name'),
                        "Start": mapping.get('start'),
                        "End": mapping.get('end'),
                        "Allele_String": mapping.get('allele_string'),
                        "Strand": mapping.get('strand')
                    })
        else:
            print(f"Error fetching chunk starting at index {i}: Status code {response.status_code}")
            print(f"Response: {response.text}")  # <-- see the actual reason
            print(f"First few IDs in chunk: {chunk[:5]}")
        
        time.sleep(0.2)  # be polite to the API, avoid rate-limit 429s later
            
    return pd.DataFrame(results)

# Read and clean the rsid list properly
df = pd.read_csv('for_additional_anno.txt', header=None)  # set header=0 if there IS a header row
my_rsids = df[0].dropna().astype(str).str.strip().unique().tolist()

# sanity check before hitting the API
print(f"Total rsids: {len(my_rsids)}")
print(f"Sample: {my_rsids[:5]}")
bad = [r for r in my_rsids if not r.startswith('rs')]
print(f"Malformed entries (should be empty): {bad[:10]}")

df_37 = lookup_rsids_grch37(my_rsids)
df_37.to_csv("rsids_grch37.csv", index=False)
```

```bash
# had to remove all non usual chromosomes
awk -F, '$2 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/' rsids_grch37.csv  > rsids_grch37_clean.csv

# combining them together:
# From the liftOver BED: strip "chr" prefix, position = end column (1-based)
awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print $4, $1, $3}' my_d13_to37.bed > from_liftover.txt

# From the Ensembl lookup CSV: rsid, chr, pos (already 1-based)
awk -F, 'BEGIN{OFS="\t"} {print $1, $2, $3}' rsids_grch37_clean.csv > from_ensembl.txt

# Combine, check for any rsid appearing in both (shouldn't happen, but verify)
cat from_liftover.txt from_ensembl.txt | sort -k1,1 > combined_grch37.txt
awk '{print $1}' combined_grch37.txt | uniq -d > overlap_check.txt
wc -l overlap_check.txt

# use as dictionary
awk 'BEGIN{OFS="\t"} {print $1, $2}' combined_grch37.txt > update_chr.txt   # rsid, new chromosome
awk 'BEGIN{OFS="\t"} {print $1, $3}' combined_grch37.txt > update_pos.txt  # rsid, new position

plink --bfile d13 \
  --update-chr update_chr.txt \
  --update-map update_pos.txt \
  --make-bed --out d13_37
```

181595 snps in d13 became 181517 snps after liftover.

Reformatting GWAMA to PLINK:
```
# Build chr:pos key from your target data (using YOUR SNP id)
awk '{print $1":"$4, $2}' d13_37.bim | sort -k1,1 > bim_key.txt

# Build chr:pos key from GWAMA (carry effect_allele, other_allele, odds_ratio, p_value)
awk 'NR>1{print $2":"$3, $4, $5, $7, $12}' GWAMA.Asia.Final.gcFEchbp.out | sort -k1,1 > gwama_key.txt

# Join on chr:pos
join -1 1 -2 1 bim_key.txt gwama_key.txt > merged.txt
wc -l merged.txt

awk 'BEGIN{print "SNP","A1","BETA"} {print $2, $3, log($5)}' merged.txt > asia_score_input.txt

awk 'BEGIN{print "SNP","P"} {print $2, $6}' merged.txt > asia_clump_pvals.txt
```

Clumping:
```
plink --bfile d13_37 \
  --clump asia_clump_pvals.txt \
  --clump-p1 0.05 --clump-p2 0.05 \
  --clump-r2 0.1 --clump-kb 250 \
  --clump-snp-field SNP --clump-field P \
  --out asia_clumped
```

3602 clumps formed from 8488 top variants.

Extract snps:
```
awk 'NR>1{print $3}' asia_clumped.clumped > asia_clumped.clumped.snplist
wc -l asia_clumped.clumped.snplist
```

Run actual scoring:
```
plink --bfile d13_37 \
  --score asia_score_input.txt 1 2 3 header \
  --extract asia_clumped.clumped.snplist \
  --out prs_asia_base
```

Sanity check the output:
```
head prs_asia_base.profile
wc -l prs_asia_base.profile
```


Getting insights:
```
plink --bfile d13_37 --indep-pairwise 200 50 0.25 --out pruned
plink --bfile d13_37 --extract pruned.prune.in --pca 10 --out d13_pca
```

```R
install.packages(c("dplyr", "pROC"))

library(dplyr)
library(pROC)

prs <- read.table("prs_asia_base.profile", header=TRUE)
pca <- read.table("d13_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID","IID", paste0("PC",1:10))

df <- merge(prs, pca, by=c("FID","IID"))

# Standardize PRS
df$PRS_z <- scale(df$SCORE)

# PHENO in .fam/.profile is coded 1=control, 2=case by default -- convert to 0/1
df$TB <- df$PHENO - 1

# Logistic regression: PRS alone
m1 <- glm(TB ~ PRS_z, data=df, family=binomial)
summary(m1)

# Logistic regression: PRS + ancestry PCs (adjust N of PCs as needed, check scree/variance explained)
m2 <- glm(TB ~ PRS_z + PC1 + PC2 + PC3, data=df, family=binomial)
summary(m2)

# AUC
pred <- predict(m2, type="response")
roc_obj <- roc(df$TB, pred)
auc(roc_obj)

# Nagelkerke pseudo-R^2
nagelkerke_r2 <- function(model) {
  n <- nobs(model)
  ll_full <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(update(model, . ~ 1)))
  
  cox_snell <- 1 - exp((2/n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2/n) * ll_null)
  nagelkerke <- cox_snell / max_r2
  
  cat("Cox & Snell R2:", round(cox_snell, 4), "\n")
  cat("Nagelkerke R2:", round(nagelkerke, 4), "\n")
  return(nagelkerke)
}

nagelkerke_r2(m2)
```

Output:
```
Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  0.04934    0.15684   0.315    0.753
PRS_z        0.17925    0.15929   1.125    0.260

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 227.25  on 163  degrees of freedom
Residual deviance: 225.97  on 162  degrees of freedom
AIC: 229.97

Number of Fisher Scoring iterations: 4


Call:
glm(formula = TB ~ PRS_z + PC1 + PC2 + PC3, family = binomial, 
    data = df)

Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  0.05089    0.15786   0.322    0.747
PRS_z        0.17697    0.15987   1.107    0.268
PC1          1.67725    2.04896   0.819    0.413
PC2          2.18933    2.06289   1.061    0.289
PC3          1.11049    2.05102   0.541    0.588

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 227.25  on 163  degrees of freedom
Residual deviance: 223.88  on 159  degrees of freedom
AIC: 233.88

Number of Fisher Scoring iterations: 4

Setting levels: control = 0, case = 1
Setting direction: controls < cases
Area under the curve: 0.5854
Cox & Snell R2: 0.0204 
Nagelkerke R2: 0.0272 
[1] 0.02719539
```

Looking at results:
```
> df$prs_decile <- ntile(df$PRS_z, 10)
table(df$prs_decile, df$TB)

# odds ratio: top decile vs bottom decile
top_vs_bottom <- df %>% filter(prs_decile %in% c(1,10))
m3 <- glm(TB ~ factor(prs_decile), data=top_vs_bottom, family=binomial)
summary(m3)
    
      0  1
  1  10  7
  2   8  9
  3   9  8
  4  10  7
  5   7  9
  6   8  8
  7  10  6
  8   4 12
  9   7  9
  10  7  9

Call:
glm(formula = TB ~ factor(prs_decile), family = binomial, data = top_vs_bottom)

Coefficients:
                     Estimate Std. Error z value Pr(>|z|)
(Intercept)           -0.3567     0.4928  -0.724    0.469
factor(prs_decile)10   0.6080     0.7049   0.863    0.388

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 45.717  on 32  degrees of freedom
Residual deviance: 44.965  on 31  degrees of freedom
AIC: 48.965

Number of Fisher Scoring iterations: 4
```


Same for europe:
```
# Join on chr:pos
awk 'NR>1{print $2":"$3, $4, $5, $7, $12}' GWAMA.Europe.Final.gcFEchbp.out | sort -k1,1 > europe_key.txt
join -1 1 -2 1 bim_key.txt europe_key.txt > merged_europe.txt
wc -l merged_europe.txt

# Check for duplicate SNP IDs from multi-allelic matches
awk '{print $2}' merged_europe.txt | sort | uniq -d | wc -l

# Build score input and clump p-value files
awk 'BEGIN{print "SNP","A1","BETA"} {print $2, $3, log($5)}' merged_europe.txt > europe_score_input.txt
awk 'BEGIN{print "SNP","P"} {print $2, $6}' merged_europe.txt > europe_clump_pvals.txt

# Clump
plink --bfile d13_37 \
  --clump europe_clump_pvals.txt \
  --clump-p1 0.05 --clump-p2 0.05 \
  --clump-r2 0.1 --clump-kb 250 \
  --clump-snp-field SNP --clump-field P \
  --out europe_clumped

# Extract clumped SNP list
awk 'NR>1{print $3}' europe_clumped.clumped > europe_clumped.clumped.snplist
wc -l europe_clumped.clumped.snplist

# Score
plink --bfile d13_37 \
  --score europe_score_input.txt 1 2 3 header \
  --extract europe_clumped.clumped.snplist \
  --out prs_europe_base
```

output:
```
Among remaining phenotypes, 84 are cases and 80 are controls.
Warning: 172350 lines skipped in --score file (172349 due to variant ID
mismatch, 1 due to allele code mismatch); see prs_europe_base.nopred for
details.
--score: 3849 valid predictors loaded.
--score: Results written to prs_europe_base.profile .
```




Getting insights for europe:
```R
library(dplyr)
library(pROC)

prs <- read.table("prs_europe_base.profile", header=TRUE)
pca <- read.table("d13_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID","IID", paste0("PC",1:10))

df <- merge(prs, pca, by=c("FID","IID"))

# Standardize PRS
df$PRS_z <- scale(df$SCORE)

# PHENO in .fam/.profile is coded 1=control, 2=case by default -- convert to 0/1
df$TB <- df$PHENO - 1

# Logistic regression: PRS alone
m1 <- glm(TB ~ PRS_z, data=df, family=binomial)
summary(m1)

# Logistic regression: PRS + ancestry PCs (adjust N of PCs as needed, check scree/variance explained)
m2 <- glm(TB ~ PRS_z + PC1 + PC2 + PC3, data=df, family=binomial)
summary(m2)

# AUC
pred <- predict(m2, type="response")
roc_obj <- roc(df$TB, pred)
auc(roc_obj)

# Nagelkerke pseudo-R^2
nagelkerke_r2 <- function(model) {
  n <- nobs(model)
  ll_full <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(update(model, . ~ 1)))
  
  cox_snell <- 1 - exp((2/n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2/n) * ll_null)
  nagelkerke <- cox_snell / max_r2
  
  cat("Cox & Snell R2:", round(cox_snell, 4), "\n")
  cat("Nagelkerke R2:", round(nagelkerke, 4), "\n")
  return(nagelkerke)
}

nagelkerke_r2(m2)
```

Comparing asia and europe:
```
library(dplyr)
library(pROC)

pca <- read.table("d13_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID","IID", paste0("PC",1:10))

prs_asia <- read.table("prs_asia_base.profile", header=TRUE)
prs_eur  <- read.table("prs_europe_base.profile", header=TRUE)

df_combined <- merge(prs_asia[,c("FID","IID","PHENO","SCORE")], 
                      prs_eur[,c("FID","IID","SCORE")], 
                      by=c("FID","IID"), 
                      suffixes=c("_asia","_europe"))

df_combined <- merge(df_combined, pca, by=c("FID","IID"))

df_combined$TB <- df_combined$PHENO - 1
df_combined$PRS_asia_z <- as.numeric(scale(df_combined$SCORE_asia))
df_combined$PRS_eur_z  <- as.numeric(scale(df_combined$SCORE_europe))

# sanity check these are actually different
cor(df_combined$PRS_asia_z, df_combined$PRS_eur_z)
head(df_combined[,c("PRS_asia_z","PRS_eur_z")])

m_asia <- glm(TB ~ PRS_asia_z + PC1 + PC2 + PC3, data=df_combined, family=binomial)
m_eur  <- glm(TB ~ PRS_eur_z + PC1 + PC2 + PC3, data=df_combined, family=binomial)

roc_asia <- roc(df_combined$TB, predict(m_asia, type="response"))
roc_eur  <- roc(df_combined$TB, predict(m_eur, type="response"))

auc(roc_asia)
auc(roc_eur)
roc.test(roc_asia, roc_eur)
```

For a writeup:  
To evaluate whether polygenic risk scores (PRS) for tuberculosis (TB) susceptibility, derived from published multi-ancestry GWAS data, could meaningfully stratify genetic risk in a Kazakh cohort, we constructed and validated PRS using summary statistics from the International Tuberculosis Host Genetics Consortium (ITHGC) multi-ancestry meta-analysis (Schurz et al., 2024). Because Kazakh ancestry is genetically intermediate between East Asian and European populations and was not itself represented in the ITHGC discovery cohorts, we tested two ancestry-specific base datasets in parallel — the Asian and European fixed-effects meta-analysis summary statistics (GWAMA) — to assess which showed better transferability to our target sample. After harmonizing genomic coordinates between datasets (lifting our target genotype data to GRCh37 and resolving variant identifiers via chromosome:position matching, since a substantial fraction of GWAMA variants and our array-genotyped SNPs used inconsistent ID formats), we performed linkage disequilibrium clumping (p<0.05, r²<0.1, 250kb windows) on each ancestry-specific base independently, retaining 3,602 independent variants from the Asian base (from 8,488 nominally significant SNPs prior to clumping) and 3,850 independent variants from the European base (from 10,312 nominally significant SNPs prior to clumping). PRS were calculated in PLINK as the weighted sum of risk allele dosages using effect sizes (log odds ratios) from each respective base, then standardized and evaluated in logistic regression models adjusted for ancestry principal components, with model performance assessed via AUC and Nagelkerke pseudo-R². Both ancestry-specific PRS showed a weak, non-significant association with TB case status in the expected direction (Asia-base: OR direction positive, p=0.26, AUC=0.585, Nagelkerke R²=0.027; Europe-base: p=0.50, AUC=0.572, Nagelkerke R²=0.021), and a DeLong's test comparing the two ROC curves found no statistically significant difference in predictive performance between the ancestry-specific bases (Z=0.34, p=0.74). These results suggest a modest, directionally consistent genetic signal that neither ancestry-specific PRS captures with statistical confidence in this cohort (n=164), most plausibly reflecting a combination of limited target sample size and the well-documented challenge of cross-ancestry PRS transferability for TB — a disease for which, as of this analysis, no validated PRS model exists in the PGS Catalog for any population.  

***need to rephrase to avoid plagiarism***
	                                                      |   Asia	 |   Europe  |  
SNPs shared between my data and the base (total overlap)  |  167980  |   176199  |  
SNPs passing p<0.05 ("top variants," pre-LD-pruning)	  |  8,488   | 	 10,312  |  
Final independent SNPs after clumping                     |  3,602   |   3,850   |  
