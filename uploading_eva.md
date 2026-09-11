conda create -n eva -c conda-forge -c bioconda eva-sub-cli
conda activate eva
eva-sub-cli.py --help
plink --bfile d13 --recode vcf --out tb_dataset

eva-sub-cli.py --metadata_xlsx EVA_TB.xlsx --submission_dir ./ --tasks VALIDATE

eva-sub-cli.py --metadata_json EVA_TB.json --submission_dir ./ --tasks VALIDATE
