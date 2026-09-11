```
conda create -n eva -c conda-forge -c bioconda eva-sub-cli
conda activate eva
eva-sub-cli.py --help
plink --bfile d13 --recode vcf --out tb_dataset

eva-sub-cli.py --metadata_xlsx EVA_TB.xlsx --submission_dir ./ --tasks VALIDATE

eva-sub-cli.py --metadata_json EVA_TB.json --submission_dir ./ --tasks VALIDATE
```

```
awk '
BEGIN{
  map["NC_000001.11"]="1"
  map["NC_000002.12"]="2"
  map["NC_000003.12"]="3"
  map["NC_000004.12"]="4"
  map["NC_000005.10"]="5"
  map["NC_000006.12"]="6"
  map["NC_000007.14"]="7"
  map["NC_000008.11"]="8"
  map["NC_000009.12"]="9"
  map["NC_000010.11"]="10"
  map["NC_000011.10"]="11"
  map["NC_000012.12"]="12"
  map["NC_000013.11"]="13"
  map["NC_000014.9"]="14"
  map["NC_000015.10"]="15"
  map["NC_000016.10"]="16"
  map["NC_000017.11"]="17"
  map["NC_000018.10"]="18"
  map["NC_000019.10"]="19"
  map["NC_000020.11"]="20"
  map["NC_000021.9"]="21"
  map["NC_000022.11"]="22"
  map["NC_000023.11"]="X"
  map["NC_000024.10"]="Y"
  map["NC_012920.1"]="MT"
}
# When we see a header line (starts with ">")
/^>/{
  # extract just the first token (e.g. >NC_000001.10)
  split($1, a, ">");
  old=a[2];
  if(old in map) {
    print ">"map[old];
  } else {
    print $0;  # keep header unchanged if not in map
  }
  next;
}
# For sequence lines, print as-is
{ print }
' GRCh38_genome.fa > GRCh38_v2.fa
```

bcftools +fixref tb_dataset.vcf -- -f GRCh38_v2.fa
bcftools norm --check-ref e -f GRCh38_v2.fa tb_dataset.vcf -Ou -o /dev/null
