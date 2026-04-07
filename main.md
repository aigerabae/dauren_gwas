```bash
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

```bash
bgzip vcf/*vcf*
for f in vcf/*.vcf.gz; do tabix -p vcf -f $f;done
bcftools merge /home/aygera/biostar/dauren_gwas/vcf/*.vcf.gz -o dauren_gwas.vcf
```

```bash
echo 'export PATH=/home/aygera/tools/plink1.9/:$PATH' >> ~/.bashrc
echo 'export PATH=/home/aygera/tools/plink2/:$PATH' >> ~/.bashrc
source ~/.bashrc
```

```bash
cd plink
plink -vcf ../dauren_gwas.vcf --double-id --pheno ../phenotypes.tsv --make-bed --out d1
awk -F'\t' '!seen[$1]++' ../rsid_conversion/InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > ../rsid_conversion/IIA-dictionary.txt
plink --bfile d1 --update-name ../rsid_conversion/IIA-dictionary.txt --make-bed --out d2
plink --bfile d2 --impute-sex ycount --make-bed --out d3
```

tr -s ' ' < d3.sexcheck > sexcheck.tsv

```python
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv('sexcheck.tsv',sep=" ")

plt.figure()
plt.scatter(df['F'], df['YCOUNT'])
plt.title("F vs. YCOUNT")
plt.xlabel("Y")
plt.ylabel("YCOUNT")
plt.savefig('f_vs_y.png')

plt.figure()
df['F'].hist(bins=20)
plt.title('F Distribution')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.savefig('f_hist.png') 

plt.figure()
df['YCOUNT'].hist(bins=20)
plt.title('Y count Distribution')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.savefig('ycount_hist.png')
```

I need to come up with thresholds to impute sex 

```bash
plink --bfile d2 --geno 0.05 --make-bed --out d3 --allow-no-sex
plink --bfile d3 --mind 0.02 --make-bed --out d4 --allow-no-sex
plink --bfile d4 --genome --min 0.2 --out pihat_min0.2 --allow-no-sex
plink --bfile d4 --missing --out missing_report --allow-no-sex
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt

echo "206667660001_R07C01
206752750007_R11C02
206767120002_R02C01
206767120002_R08C01
206767120002_R08C02
206767120003_R02C02
206767120003_R12C01" > related_remove.txt

cat related_remove.txt | awk '{print $1 "\t" $1}' > to_remove.txt

plink --bfile d4 --remove to_remove.txt --allow-no-sex --make-bed --out d5
plink --bfile d5 --maf 0.001 --make-bed --out d6
```

```bash
plink2 --bfile d6 --pca 10 --out pca
plink --bfile d6 --covar pca.eigenvec --allow-no-sex  --model --out model_mafs_filtered
cat model_mafs_filtered.model | awk '$10 != "NA" && $10 < 1e-3' | sort -gk 9,9
```
