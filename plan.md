~~1) genotype calling~~  
~~2) confirm gtc isn't silently dropping any samples~~  
~~4) Check REF/ALT consistency and strand orientation against the reference FASTA before merge (bcftools norm -check-ref).~~  
~~5) get stats for my vcf file (snps, indels, etc.)~~  
~~6) Compress merged output: bgzip dauren_gwas.vcf + index.~~  
~~7) splitting multiallelics~~  
~~8) annotation of variants - rsid and genes~~  
~~9) Filter/flag on WGT, SSR, and WTD - withdrawn, non-unique-mapping, and suspect variants~~  
~~10) remove full duplicates; identify duplicates by rsid and assign them different names if ref/alt/position are slightly different~~  
~~11) remove Ambiguous strand SNPs (A/T, C/G)~~  
12) Genotyping call rate per sample (--mind) — you have this, keep 0.02–0.05.  
13) Heterozygosity outliers — compute --het, exclude samples >±3-4 SD from mean heterozygosity (catches contamination/inbreeding artifacts). You're missing this entirely.  
14) Sex check — apply a documented threshold (e.g., F<0.2 = female, F>0.8 = male for X-het; or YCOUNT-based cutoff from your histogram — pick the valley between the two clusters) and flag/exclude discordant sex.  
15) Relatedness (--genome, KING, or plink2 --king-cutoff) — use a systematic pruning algorithm, not manual ID lists. Document the pi-hat threshold rationale (0.2 is reasonable for excluding 2nd-degree+).  
16) Population structure / ancestry outliers — project your samples onto 1000 Genomes or a reference panel via PCA to identify and exclude ancestry outliers before running the GWAS PCA. This is very likely missing and is a common reviewer objection.  
17) Call rate (--geno) — you have this.  
18) HWE filter — apply in controls only, threshold typically --hwe 1e-6 (stricter for cases can be inappropriate since disease-associated loci may deviate from HWE by design).  
19) MAF filter — 0.001 is very permissive; for a few hundred samples you have essentially no power to detect anything below MAF ~0.05 reliably. Consider MAF ≥0.01 or ≥0.05 for the primary GWAS and treat rarer variants separately with burden/SKAT-type tests if relevant; MAF 2×N alleles → minimum detectable non-zero MAF = 1/(2N)  
20) Differential missingness between cases/controls (--test-missing) — catches genotyping artifacts that masquerade as association signal.  
21) Recompute PCA on the final QC'd, ancestry-outlier-pruned dataset (LD-prune first: --indep-pairwise 50 5 0.2).  
22) Scree plot / eigenvalue inspection to decide how many PCs to retain as covariates (typically 5-10, but check where the elbow is).  
23) Logistic regression with covariates: sex, age (if available), top N PCs, and genotyping batch/plate if you have multiple runs.  
24) If any cryptic relatedness remains, switch to a mixed-model method (SAIGE, REGENIE, BOLT-LMM) instead of dropping related individuals, since relatedness pruning throws away sample size you probably can't spare.  
25) Consider Firth-corrected logistic regression for rare variants/small cell counts (plink2 --glm firth-fallback).  
26) λ_GC (genomic inflation factor) — should be close to 1.0; report it.  
27) QQ plot of observed vs expected p-values.  
28) Manhattan plot with genome-wide (5×10⁻⁸) and suggestive (1×10⁻⁵) significance lines.  
29) If imputing to genome-wide density, run imputation QC (INFO/R² score filtering, typically >0.3 or >0.8 for high-confidence variants).  
30) Regional association (LocusZoom-style) plots for any hits approaching significance.  
31) Functional annotation (VEP, ANNOVAR) and cross-reference with GWAS Catalog / immune-disease databases given this is an immune-focused array.  
32) eQTL/colocalization check if resources allow.  
33) Data/code availability — deposit summary statistics (GWAS Catalog), and ideally a reproducible pipeline (Snakemake/Nextflow) with versioned tools/containers, not a loose shell script history like this one.  
34) Clear justification of the array's scope — since this is a targeted immune-focused chip, be explicit in Methods about coverage limitations, and justify genome-wide claims only if you imputed to full coverage.  
35) Multiple testing transparency — state the significance threshold used and why (genome-wide 5×10⁻�8, or Bonferroni-adjusted for a targeted panel size if you keep it as array-only). 12) compute an effective number of independent tests (via LD-based correlation matrix, e.g. simpleM or PLINK --indep-pairwise derived eigenvalue counting) instead of blanket 5×10⁻⁸ applies directly to your model_mafs_filtered output  
36) Sample size/power justification — a priori or post hoc power calculation (e.g., via genpwr or CaTS).  
37) optional from me: use KGWAS to increase power: https://pmc.ncbi.nlm.nih.gov/articles/PMC11643201/  
38) sample/variant count log per QC step - halfway through  
39) Version-pin array-analysis-cli, bcftools, plink/plink2 and record in methods (reproducibility).  

