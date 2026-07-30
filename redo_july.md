1) use plink2 --king-cutoff instead of PI_HAT
2) use MAF 2×N alleles → minimum detectable non-zero MAF = 1/(2N)
3) remove Ambiguous strand SNPs (A/T, C/G)
4) duplicate-variant removal
5) excess-heterozygosity sample filter
6) sample/variant count log per QC step
7) compute an effective number of independent tests (via LD-based correlation matrix, e.g. simpleM or PLINK --indep-pairwise derived eigenvalue counting) instead of blanket 5×10⁻⁸ applies directly to your model_mafs_filtered output
8) PCA computed but not used
9) HWE filter
10) Sex imputation threshold undefined (Ad hoc relatedness removal — manually pasting 7 sample IDs into related_remove.txt with no documented rationale for why those 7 (vs. keeping the other member of each related pair) isn't defensible in a methods section.
11) --model is the wrong association test — for a Q1 paper you need logistic regression with covariates (sex, age, PCs, batch/plate if relevant), not the basic chi-square/trend/genotypic tests in --model. Use plink --logistic / plink2 --glm firth-fallback (or better, SAIGE/REGENIE if you have any relatedness or rare variants).
multiple testing correction framework
12) This is an Immuno BeadChip, not a genome-wide array - reframe as a targeted/candidate-region association study
13) No inflation/stratification diagnostics — no λ_GC, no QQ plot, no Manhattan plot anywhere in the pipeline.
14) Small apparent sample size — your top hit shows ~83 cases/174 controls-ish splits; that's underpowered for genome-wide significance for most effect sizes. Worth running a formal power calculation before you invest more in analysis.
15) Version-pin array-analysis-cli, bcftools, plink/plink2 and record in methods (reproducibility).
16) Check REF/ALT consistency and strand orientation against the reference FASTA before merge (bcftools norm -check-ref).
17) Compress merged output: bgzip dauren_gwas.vcf + index.
18) optional from me: use KGWAS to increase power: https://pmc.ncbi.nlm.nih.gov/articles/PMC11643201/


New workflow:
A. Genotype calling & VCF conversion (mostly fine, tighten these):

Confirm --output-folder ./gtc bug isn't silently dropping/renaming samples — diff sample counts between sample sheet and GTC output.
Version-pin array-analysis-cli, bcftools, plink/plink2 and record in methods (reproducibility).

B. VCF merge

Check REF/ALT consistency and strand orientation against the reference FASTA before merge (bcftools norm -check-ref).
Compress merged output: bgzip dauren_gwas.vcf + index.

C. Sample-level QC (order matters — do this before variant QC)

Genotyping call rate per sample (--mind) — you have this, keep 0.02–0.05.
Heterozygosity outliers — compute --het, exclude samples >±3-4 SD from mean heterozygosity (catches contamination/inbreeding artifacts). You're missing this entirely.
Sex check — apply a documented threshold (e.g., F<0.2 = female, F>0.8 = male for X-het; or YCOUNT-based cutoff from your histogram — pick the valley between the two clusters) and flag/exclude discordant sex.
Relatedness (--genome, KING, or plink2 --king-cutoff) — use a systematic pruning algorithm, not manual ID lists. Document the pi-hat threshold rationale (0.2 is reasonable for excluding 2nd-degree+).
Population structure / ancestry outliers — project your samples onto 1000 Genomes or a reference panel via PCA to identify and exclude ancestry outliers before running the GWAS PCA. This is very likely missing and is a common reviewer objection.

D. Variant-level QC

Call rate (--geno) — you have this.
HWE filter — apply in controls only, threshold typically --hwe 1e-6 (stricter for cases can be inappropriate since disease-associated loci may deviate from HWE by design).
MAF filter — 0.001 is very permissive; for a few hundred samples you have essentially no power to detect anything below MAF ~0.05 reliably. Consider MAF ≥0.01 or ≥0.05 for the primary GWAS and treat rarer variants separately with burden/SKAT-type tests if relevant.
Differential missingness between cases/controls (--test-missing) — catches genotyping artifacts that masquerade as association signal.

E. Population structure correction

Recompute PCA on the final QC'd, ancestry-outlier-pruned dataset (LD-prune first: --indep-pairwise 50 5 0.2).
Scree plot / eigenvalue inspection to decide how many PCs to retain as covariates (typically 5-10, but check where the elbow is).

F. Association testing

Logistic regression with covariates: sex, age (if available), top N PCs, and genotyping batch/plate if you have multiple runs.
If any cryptic relatedness remains, switch to a mixed-model method (SAIGE, REGENIE, BOLT-LMM) instead of dropping related individuals, since relatedness pruning throws away sample size you probably can't spare.
Consider Firth-corrected logistic regression for rare variants/small cell counts (plink2 --glm firth-fallback).

G. Post-GWAS diagnostics (required for any Q1 submission)

λ_GC (genomic inflation factor) — should be close to 1.0; report it.
QQ plot of observed vs expected p-values.
Manhattan plot with genome-wide (5×10⁻⁸) and suggestive (1×10⁻⁵) significance lines.
If imputing to genome-wide density, run imputation QC (INFO/R² score filtering, typically >0.3 or >0.8 for high-confidence variants).

H. Follow-up on top hits

Regional association (LocusZoom-style) plots for any hits approaching significance.
Functional annotation (VEP, ANNOVAR) and cross-reference with GWAS Catalog / immune-disease databases given this is an immune-focused array.
eQTL/colocalization check if resources allow.
Replication in an independent cohort — for a Q1 immunology/genetics journal this is close to mandatory for any novel locus claim, or at minimum explicit acknowledgment of the discovery-only nature with appropriate caveats.
3. What Q1 reviewers/editors will specifically expect
STREGA/STROBE reporting checklist compliance (genetic epidemiology reporting standard) — most gen/epi journals require this at submission.
Sample size/power justification — a priori or post hoc power calculation (e.g., via genpwr or CaTS).
Ethics statement — IRB approval and informed consent documentation.
Data/code availability — deposit summary statistics (GWAS Catalog), and ideally a reproducible pipeline (Snakemake/Nextflow) with versioned tools/containers, not a loose shell script history like this one.
Clear justification of the array's scope — since this is a targeted immune-focused chip, be explicit in Methods about coverage limitations, and justify genome-wide claims only if you imputed to full coverage.
Multiple testing transparency — state the significance threshold used and why (genome-wide 5×10⁻�8, or Bonferroni-adjusted for a targeted panel size if you keep it as array-only).
