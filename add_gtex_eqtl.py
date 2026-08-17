#!/usr/bin/env python3
"""
add_gtex_eqtl.py

Adds GTEx cis-eQTL annotation columns to an existing annotated SNP table
(the kind of TSV you get from VEP + gnomAD + GWAS Catalog annotation).

Reads your table (must contain RSID, CHR, POS, REF, ALT columns), queries
the GTEx Portal API v2 for each variant, and writes a new TSV with the
original columns preserved plus GTEx columns appended at the end.

WORKFLOW PER VARIANT
---------------------
1. rsID -> GTEx variant ID  (GET /api/v2/dataset/variant?snpId=...)
   GTEx variant IDs look like: chr5_175833369_G_A_b38
2. GTEx variant ID -> precomputed significant single-tissue eQTLs
   (GET /api/v2/association/singleTissueEqtl?variantId=...&datasetId=gtex_v8)

BUILD NOTE
----------
GTEx v8 (datasetId="gtex_v8") is built on GRCh38. If your original CHR/POS
are GRCh38 (VEP_Assembly_x = GRCh38 in your header, so this looks correct
for your data), step 1 (rsID lookup) handles the coordinate translation for
you internally - you don't need to convert anything yourself. If a future
GTEx release changes the default dataset, update GTEX_DATASET_ID below.

USAGE
-----
    python add_gtex_eqtl.py --input annotated_snps.tsv --output annotated_snps_with_gtex.tsv

    Optional:
    --rsid-col RSID       (default: RSID)
    --test 3              (only process first N rows, for a dry run)
    --sleep 1.0           (seconds between API calls; don't lower this much)

OUTPUT COLUMNS ADDED
---------------------
    GTEx_VariantId              GTEx-format variant ID (or "NA" if not found)
    GTEx_eQTL_Significant        "Yes" / "No" / "NA (lookup failed)"
    GTEx_eQTL_Gene_Symbols        semicolon-separated genes this variant is an eQTL for
    GTEx_eQTL_Tissues             semicolon-separated tissues with a significant hit
    GTEx_eQTL_Min_PValue          smallest p-value across all tissue/gene hits
    GTEx_eQTL_Details             full semicolon list: GENE@TISSUE(p=...,nes=...)
    GTEx_eQTL_Hit_Count           number of significant tissue/gene eQTL hits

NOT EXECUTED AGAINST LIVE DATA: this environment has no network access to
gtexportal.org, so this script has been written carefully against the
documented API but not run end-to-end. Run with --test 3 first and inspect
the output before processing your full table.
"""

import argparse
import csv
import sys
import time
import requests

GTEX_BASE = "https://gtexportal.org/api/v2"
GTEX_DATASET_ID = "gtex_v8"


def get_gtex_variant_id(rsid, session, sleep_sec):
    """Step 1: rsID -> GTEx variant ID."""
    try:
        r = session.get(
            f"{GTEX_BASE}/dataset/variant",
            params={"snpId": rsid, "datasetId": GTEX_DATASET_ID},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"  [WARN] GTEx variant lookup failed for {rsid}: {e}", file=sys.stderr)
        time.sleep(sleep_sec)
        return None

    records = data.get("data", [])
    time.sleep(sleep_sec)
    if not records:
        return None
    return records[0].get("variantId")


def get_gtex_eqtls(gtex_variant_id, session, sleep_sec):
    """Step 2: GTEx variant ID -> significant single-tissue eQTLs."""
    try:
        r = session.get(
            f"{GTEX_BASE}/association/singleTissueEqtl",
            params={"variantId": gtex_variant_id, "datasetId": GTEX_DATASET_ID},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"  [WARN] GTEx eQTL lookup failed for {gtex_variant_id}: {e}", file=sys.stderr)
        time.sleep(sleep_sec)
        return None

    time.sleep(sleep_sec)
    return data.get("data", [])


def annotate_row(rsid, session, sleep_sec):
    """Returns a dict of the GTEx columns to add for one SNP."""
    empty = {
        "GTEx_VariantId": "NA",
        "GTEx_eQTL_Significant": "NA (no GTEx variant record)",
        "GTEx_eQTL_Gene_Symbols": "NA",
        "GTEx_eQTL_Tissues": "NA",
        "GTEx_eQTL_Min_PValue": "NA",
        "GTEx_eQTL_Details": "NA",
        "GTEx_eQTL_Hit_Count": 0,
    }

    if not rsid or rsid in (".", "NA", ""):
        empty["GTEx_eQTL_Significant"] = "NA (no rsID)"
        return empty

    gtex_variant_id = get_gtex_variant_id(rsid, session, sleep_sec)
    if gtex_variant_id is None:
        return empty

    eqtl_hits = get_gtex_eqtls(gtex_variant_id, session, sleep_sec)
    if eqtl_hits is None:
        return {**empty, "GTEx_VariantId": gtex_variant_id,
                "GTEx_eQTL_Significant": "NA (lookup failed)"}

    if len(eqtl_hits) == 0:
        return {**empty, "GTEx_VariantId": gtex_variant_id,
                "GTEx_eQTL_Significant": "No"}

    genes = sorted(set(h.get("geneSymbol", "?") for h in eqtl_hits))
    tissues = sorted(set(h.get("tissueSiteDetailId", "?") for h in eqtl_hits))
    pvals = [h.get("pValue") for h in eqtl_hits if h.get("pValue") is not None]
    min_p = min(pvals) if pvals else "NA"
    details = "; ".join(
        f"{h.get('geneSymbol','?')}@{h.get('tissueSiteDetailId','?')}"
        f"(p={h.get('pValue','NA')},nes={h.get('nes','NA')})"
        for h in eqtl_hits
    )

    return {
        "GTEx_VariantId": gtex_variant_id,
        "GTEx_eQTL_Significant": "Yes",
        "GTEx_eQTL_Gene_Symbols": ";".join(genes),
        "GTEx_eQTL_Tissues": ";".join(tissues),
        "GTEx_eQTL_Min_PValue": min_p,
        "GTEx_eQTL_Details": details,
        "GTEx_eQTL_Hit_Count": len(eqtl_hits),
    }


def main():
    parser = argparse.ArgumentParser(description="Add GTEx eQTL columns to an annotated SNP table.")
    parser.add_argument("--input", required=True, help="Path to your annotated TSV")
    parser.add_argument("--output", required=True, help="Path to write the augmented TSV")
    parser.add_argument("--rsid-col", default="RSID", help="Name of the rsID column (default: RSID)")
    parser.add_argument("--test", type=int, default=None, help="Only process the first N rows (dry run)")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds to sleep between API calls")
    args = parser.parse_args()

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        orig_fieldnames = reader.fieldnames

    if args.rsid_col not in orig_fieldnames:
        sys.exit(f"ERROR: column '{args.rsid_col}' not found in {args.input}. "
                  f"Available columns: {orig_fieldnames}")

    if args.test:
        rows = rows[:args.test]
        print(f"[TEST MODE] Only processing first {args.test} rows.")

    session = requests.Session()
    gtex_cols = [
        "GTEx_VariantId", "GTEx_eQTL_Significant", "GTEx_eQTL_Gene_Symbols",
        "GTEx_eQTL_Tissues", "GTEx_eQTL_Min_PValue", "GTEx_eQTL_Details",
        "GTEx_eQTL_Hit_Count",
    ]

    for i, row in enumerate(rows):
        rsid = row.get(args.rsid_col, "").strip()
        print(f"[{i+1}/{len(rows)}] {rsid} ...")
        gtex_data = annotate_row(rsid, session, args.sleep)
        row.update(gtex_data)

    out_fieldnames = orig_fieldnames + gtex_cols
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    n_sig = sum(1 for r in rows if r.get("GTEx_eQTL_Significant") == "Yes")
    print(f"\nDone. Wrote {len(rows)} rows to {args.output}")
    print(f"{n_sig} of {len(rows)} SNPs are significant GTEx eQTLs in at least one tissue.")


if __name__ == "__main__":
    main()
