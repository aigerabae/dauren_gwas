# Host Genetic Susceptibility to Tuberculosis in a Kazakh Cohort: Targeted Immune-Panel Association Pipeline

[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fbiomedicines14102387-blue)](https://doi.org/10.3390/biomedicines14102387)
![Status](https://img.shields.io/badge/status-published-brightgreen)

Analysis workflow behind a **peer-reviewed publication in *Biomedicines* (MDPI)**: a case-control genetic association study of tuberculosis (TB) susceptibility in Kazakh participants without conventional TB risk factors, using a targeted immune-gene genotyping array. The workflow runs from raw Illumina IDAT files through QC, population-structure control, Firth-corrected association testing, multiple-testing correction, gene-based testing, and polygenic risk scoring.

## Publication

> Serikzhan, A.; Abilova, Z.; Akhmetova, A.; Gabdulkayum, A.; Daniyarov, A.Z.; Rakhimova, S.E.; Kozhamkulov, U.; Kairov, U.; Akilzhanova, A.; Yerezhepov, D. **Genetic Susceptibility Within the Immune System Beyond Conventional Risk Factors: An Integrative Analysis in a Kazakh Cohort.** *Biomedicines* **2026**, *14*, 2387. https://doi.org/10.3390/biomedicines14102387

<details>
<summary>BibTeX</summary>

```bibtex
@article{serikzhan2026immune,
  author  = {Serikzhan, A. and Abilova, Z. and Akhmetova, A. and Gabdulkayum, A. and Daniyarov, A.Z. and Rakhimova, S.E. and Kozhamkulov, U. and Kairov, U. and Akilzhanova, A. and Yerezhepov, D.},
  title   = {Genetic Susceptibility Within the Immune System Beyond Conventional Risk Factors: An Integrative Analysis in a Kazakh Cohort},
  journal = {Biomedicines},
  year    = {2026},
  volume  = {14},
  pages   = {2387},
  doi     = {10.3390/biomedicines14102387}
}
```
</details>

## Study at a glance

| | |
|---|---|
| Design | Case-control, targeted immunogenetic association study |
| Array | Illumina Infinium Immune Array-24 v2 (GRCh38) |
| Samples | 192 genotyped, **164 after QC** (84 cases, 80 controls) |
| Variants | **181,595** after QC (180,388 autosomal) |
| Covariates | Genetically imputed sex + PC1–PC3 |
| Primary model | Firth-fallback logistic regression, additive (PLINK 2) |

## Key findings

- **No variant reached genome-wide significance** (5 × 10⁻⁸), nor the array-specific simpleM threshold (3.06 × 10⁻⁶; M_eff = 16,319.2).
- **Lead signal: rs1560036** (chr5:175,833,369): OR = 2.98, p = 2.42 × 10⁻⁵, λGC = 1.062. It was stable across sex-unadjusted, full-Firth and permutation analyses (empirical p = 1.1 × 10⁻⁵). The genome-wide max(T) permutation p was 0.656, so the signal is interpreted as hypothesis-generating.
- **Power:** at N = 164 the study is adequately powered (≥ 80%) at genome-wide significance only for large-effect common variants, so it is framed as exploratory.
- **Gene-based (MAGMA):** no gene passed the Bonferroni threshold (3.33 × 10⁻⁶ across 15,015 genes), using an ancestry-matched Kazakh WGS reference panel for LD.
- **PRS:** scores built from International TB Host Genetics Consortium summary statistics (East Asian and European bases) showed weak, non-significant discrimination (AUC 0.585 and 0.572; DeLong p = 0.74).

## Pipeline overview

| Step | What it does | Tools |
|---|---|---|
| 1 | Genotype calling (IDAT → GTC → VCF) | Illumina `array-analysis-cli` |
| 2 | Per-sample REF/ALT consistency checks | bcftools, samtools |
| 3–4 | Merge, realign, split multiallelics, rsID annotation (dbSNP), remove suspect/withdrawn, duplicate and strand-ambiguous variants | bcftools, SnpSift, awk |
| 5 | Sample QC: call rate (`--geno/--mind 0.05`), heterozygosity (±3 SD), sex check | PLINK 1.9 |
| 6 | Relatedness (KING 0.177), differential missingness (p < 1e-5), HWE in controls (p < 1e-6), MAF ≥ 0.01 | PLINK 2, PLINK 1.9 |
| 7 | LD-pruned PCA; iterative ancestry-outlier removal (±4 SD on PC1/PC2; N = 179 → 173 → 168 → 164); projection against worldwide reference populations | PLINK, Python |
| 8 | Association testing and covariate sensitivity (PC1, PC1–3, PC1–10, ± sex) | PLINK 2 (`--glm firth-fallback`) |
| 9 | Genomic inflation, QQ and Manhattan plots | Python (scipy, matplotlib) |
| 10, 13 | Multiple-testing thresholds: LD-pruned Bonferroni and simpleM effective tests | PLINK, R |
| 12 | Post-hoc power at genome-wide and suggestive thresholds | R (`genpwr`) |
| 14 | Sensitivity analyses: sex-unadjusted, full Firth, permutation, genotypic model | PLINK, PLINK 2 |
| 15 | Functional annotation of top hits, including GTEx eQTL look-up | Python |
| 16 | Polygenic risk scores (liftover to GRCh37, clumping, scoring, AUC, DeLong) | PLINK, R (`pROC`) |
| 17 | Gene-based testing with MAGMA (10/10 kb and 35/10 kb windows) | MAGMA 1.10 |

## Repository contents

- [`main_restructured.md`](main_restructured.md): the full analysis log, in order, with commands, QC thresholds, intermediate counts and output excerpts.

## Software versions

```
bcftools 1.24 · bgzip 1.24 · tabix 1.24
PLINK 1.9 and PLINK 2
Illumina array-analysis-cli 2.1.0 (linux x64)
SnpSift 5.4c
MAGMA 1.10
R (genpwr 1.0.4, pROC) · Python (pandas, scipy, matplotlib)
```

## Notes on reproducibility

- The workflow is documented as a chronological analysis log rather than a one-command pipeline. Paths (for example `~/biostar/dauren_gwas/...`) are from the original compute environment and need adapting.
- Individual-level genotype and phenotype data are **not** included, because they are human-subject data. See the paper's data availability statement for access.
- Public inputs: GRCh38 reference genome, dbSNP common variants (build 151), MAGMA GRCh38 gene locations, and ITHGC TB GWAS meta-analysis summary statistics (Schurz et al., 2024).

## Citation

If you use or adapt this workflow, please cite the paper above.

## Contact

Corresponding authors are listed in the publication. Questions about the code can be raised as a GitHub issue.
