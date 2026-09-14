```
conda create -n eva -c conda-forge -c bioconda eva-sub-cli
conda activate eva
eva-sub-cli.py --help
plink --bfile d13 --recode vcf --out tb_dataset

eva-sub-cli.py --metadata_xlsx EVA_TB.xlsx --submission_dir ./ --tasks VALIDATE
```

```
#!/usr/bin/env bash
# get_official_grch38.sh
# Downloads the official NCBI/INSDC GRCh38 primary assembly and rewrites
# the FASTA headers from GenBank accessions (CM000663.2 ...) to the plain
# chromosome names (1, 2, ... X, Y, MT) your VCF already uses.
#
# Run this on a machine with normal internet access (not sandboxed).
# Requires: wget or curl, awk, samtools

set -euo pipefail

ACCESSION="GCA_000001405.15"   # base GRCh38, primary chromosome sequences
                                 # are identical across .15/.17/.18/.20/.22/.25/.29
                                 # (patches only add extra scaffolds), so this
                                 # is safe even though your metadata says .18
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/${ACCESSION}_GRCh38"

echo "== downloading assembly report and genomic FASTA =="
wget -N "${BASE_URL}/${ACCESSION}_GRCh38_assembly_report.txt"
wget -N "${BASE_URL}/${ACCESSION}_GRCh38_genomic.fna.gz"
gunzip -f "${ACCESSION}_GRCh38_genomic.fna.gz"

echo "== building accession -> chromosome-name map from assembly report =="
# Assembly report columns (tab-separated, after comment lines):
# 1 Sequence-Name  2 Sequence-Role  3 Assigned-Molecule  4 ...  5 GenBank-Accn  ...
grep -v "^#" "${ACCESSION}_GRCh38_assembly_report.txt" \
    | awk -F'\t' '{print $5"\t"$1}' > accn_to_chr.txt

echo "== rewriting FASTA headers to plain chromosome names =="
awk -v mapfile=accn_to_chr.txt '
    BEGIN {
        while ((getline line < mapfile) > 0) {
            split(line, a, "\t");
            map[a[1]] = a[2];
        }
    }
    /^>/ {
        acc = substr($1, 2);
        if (acc in map) {
            print ">" map[acc];
        } else {
            print $0;   # leave unplaced/unlocalized scaffolds as-is
        }
        next;
    }
    { print }
' "${ACCESSION}_GRCh38_genomic.fna" > GRCh38_official.fa

samtools faidx GRCh38_official.fa

echo "== quick sanity check: primary chromosomes present =="
grep ">" GRCh38_official.fa | grep -E "^>(1|2|X|Y|MT)$" || true

echo "Done. Use GRCh38_official.fa in place of GRCh38_v2.fa from here on:"
echo
echo "  bcftools norm --check-ref e -f GRCh38_official.fa tb_dataset.vcf -Ou -o /dev/null"
echo "  bcftools +fixref tb_dataset.vcf -- -f GRCh38_official.fa"
echo
echo "If mismatches drop to ~0% (after the strand-flip step from before),"
echo "that confirms the Illumina FASTA was the actual problem. Also re-point"
echo "eva-sub-cli's Fasta File entry at GRCh38_official.fa and re-run validation."
```

```
# 1. Rename chr 23 -> X (skip if you already did this on your working copy)
cat > chr_rename_map.txt <<'EOF'
23	X
24	Y
25	X
26	MT
EOF
bcftools annotate --rename-chrs chr_rename_map.txt tb_dataset.vcf -Oz -o tb_dataset.chrfix.vcf.gz
tabix -p vcf tb_dataset.chrfix.vcf.gz

# 2. Actually apply the strand flip this time (-m flip), and drop unresolvable sites (-d)
bcftools +fixref tb_dataset.chrfix.vcf.gz -- -f GRCh38_official.fa -m flip -d \
    | bcftools sort -Oz -o tb_dataset.fixref.vcf.gz
tabix -p vcf tb_dataset.fixref.vcf.gz

# 3. Re-check the match rate on the corrected file
bcftools norm --check-ref e -f GRCh38_official.fa tb_dataset.fixref.vcf.gz -Ou -o /dev/null
```

let's refresh what we did here: 
1) i fixed Column "Taxonomy ID" is not populated by adding it manually into my metadata;  
2) i fixed sample names to be 207859430008_R02C01 and not 207859430008_R02C01_207859430008_R02C01 in the vcf (kept the same name tb_dataset.vcf);  
3) changed assemby accession in metadata from GCF_000001405.40 to GCA_000001405.15  
4) Assembly check - ref alt not matching - flipped those alleles  
5) Warning: Non-GCA reference found in metadata. Please provide the INSDC accession for your reference assembly Some sequences are not INSDC accessioned. For my dataset i used grch38.fa file supplied by illumina and i initially tried to verify my dataset against it but I decided to use legit GCA assembly instead and changed ref file in the metadata to be GRCh38_official.fa  

Validating again:
```
eva-sub-cli.py --metadata_xlsx EVA_TB.xlsx --submission_dir ./ --tasks VALIDATE
```
