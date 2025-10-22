#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge per-project metadata CSV/TSV files into:
  1) A single Excel workbook with:
        - One sheet per project (PRJNA*)
        - An 'ALL' sheet (all rows merged)
        - A 'SUMMARY' sheet (counts & coverage)
  2) Optionally, one merged CSV of all rows

Key features
------------
- Header-safe: keeps each file's actual columns; union across projects
- Handles CSV or TSV automatically (delimiter inference)
- Preserves text exactly (dtype=str), no type coercion
- Stable column order: starts from a preferred header order (if provided)
  then appends any unseen columns in discovery order

Usage
-----
python merge_project_metadata.py \
  --root /mnt/18T/chibao/gliomas/data/metadata_full \
  --pattern "*/*_metadata.tsv" \
  --out-xlsx /mnt/18T/chibao/gliomas/data/metadata_full/ALL_projects_metadata.xlsx \
  --out-csv  /mnt/18T/chibao/gliomas/data/metadata_full/ALL_projects_metadata.tsv

Notes
-----
- Excel sheet names are truncated to 31 chars (Excel limit) and deduplicated if needed.
- Very large datasets: Excel writing is slower; CSV is recommended for bulk pipelines.
"""

import argparse
import sys
import re
from pathlib import Path
from collections import OrderedDict, defaultdict

import pandas as pd

# --------------------------
# Preferred header ordering:
# (We will keep this order first, then add any unseen columns at the end)
PREFERRED_HEADERS = [
    "run_accession","study_accession","study_title","experiment_accession","experiment_title",
    "experiment_desc","organism_taxid","organism_name","library_name","library_strategy",
    "library_source","library_selection","library_layout","sample_accession","sample_title",
    "biosample","bioproject","instrument","instrument_model","instrument_model_desc",
    "total_spots","total_size","run_total_spots","run_total_bases","run_alias",
    "public_filename","public_size","public_date","public_md5","public_version",
    "public_semantic_name","public_supertype","public_sratoolkit","aws_url","aws_free_egress",
    "aws_access_type","public_url","ncbi_url","ncbi_free_egress","ncbi_access_type","gcp_url",
    "gcp_free_egress","gcp_access_type","experiment_alias","source_name","cell type","treatment",
    "time point","ena_fastq_http","ena_fastq_http_1","ena_fastq_http_2","ena_fastq_ftp",
    "ena_fastq_ftp_1","ena_fastq_ftp_2","disease","age","sex","tissue","cell_type"
]
# --------------------------


def infer_project_id(path: Path) -> str:
    """
    Extract project id (e.g., PRJNA578617) from the path.
    Assumes structure .../metadata_full/PRJNA*/... or falls back to parent name.
    """
    parts = path.parts
    for p in reversed(parts):
        if re.match(r"PRJ[ENSDGXA]\d+", p):
            return p
    # fallback: folder name one up
    return path.parent.name


def infer_sep(fname: Path):
    """
    Ask pandas to infer separator; favor TSV and CSV transparently.
    Using sep=None with engine='python' lets pandas sniff.
    """
    return None, "python"  # sep=None, engine='python'


def safe_sheet_name(name: str, used: set) -> str:
    """
    Excel sheet name <= 31 chars; deduplicate if needed.
    """
    base = name[:31]
    candidate = base
    i = 1
    while candidate in used:
        suffix = f"_{i}"
        candidate = (base[:31 - len(suffix)]) + suffix
        i += 1
    used.add(candidate)
    return candidate


def main():
    ap = argparse.ArgumentParser(description="Merge per-project metadata into Excel (sheets) and merged CSV.")
    ap.add_argument("--root", required=True, help="Root folder containing project subfolders (metadata_full)")
    ap.add_argument("--pattern", default="*/*_metadata.tsv",
                    help="Glob under root (e.g., '*/*_metadata.tsv' or '*/*_metadata.csv')")
    ap.add_argument("--out-xlsx", required=True, help="Path to output Excel workbook (*.xlsx)")
    ap.add_argument("--out-csv", default=None, help="(Optional) Path to write merged ALL rows (CSV/TSV inferred by extension)")
    ap.add_argument("--encoding", default="utf-8", help="Character encoding (default: utf-8)")
    ap.add_argument("--add-project-col", action="store_true",
                    help="Add a 'project_id' column to frames and ALL (default: off)")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    files = sorted(root.glob(args.pattern))
    if not files:
        print(f"[ERROR] No files matched {root}/{args.pattern}", file=sys.stderr)
        sys.exit(2)

    print(f"[INFO] Found {len(files)} file(s). Building column union...")

    # ----- First pass: collect column union and per-file col sets -----
    all_columns_order = OrderedDict((col, True) for col in PREFERRED_HEADERS)
    per_file_cols = {}
    project_to_files = defaultdict(list)

    for f in files:
        sep, eng = infer_sep(f)
        try:
            # read only header row cheaply
            df_head = pd.read_csv(f, sep=sep, engine=eng, nrows=0, dtype=str, encoding=args.encoding)
        except Exception as e:
            print(f"[WARN] Skipping {f}: {e}", file=sys.stderr)
            continue

        cols = list(df_head.columns)
        # register any columns not in PREFERRED_HEADERS
        for c in cols:
            if c not in all_columns_order:
                all_columns_order[c] = True

        per_file_cols[f] = cols
        project_to_files[infer_project_id(f)].append(f)

    union_cols = list(all_columns_order.keys())
    print(f"[INFO] Unioned columns: {len(union_cols)} total")

    # ----- Second pass: read full dataframes and align to union -----
    per_project_frames = {}
    project_counts = []
    for project, flist in sorted(project_to_files.items()):
        frames = []
        for f in flist:
            sep, eng = infer_sep(f)
            try:
                df = pd.read_csv(f, sep=sep, engine=eng, dtype=str, encoding=args.encoding)
            except Exception as e:
                print(f"[WARN] Could not read {f}: {e}", file=sys.stderr)
                continue

            # Reindex to union cols (missing columns → NaN)
            df = df.reindex(columns=union_cols)
            if args.add_project_col:
                df.insert(0, "project_id", project)
            frames.append(df)

        if not frames:
            print(f"[WARN] No readable frames for {project}, skipping sheet.", file=sys.stderr)
            continue

        proj_df = pd.concat(frames, ignore_index=True)
        per_project_frames[project] = proj_df
        project_counts.append({
            "project_id": project,
            "n_files": len(flist),
            "n_rows": int(proj_df.shape[0]),
            "n_cols": int(proj_df.shape[1])
        })
        print(f"[INFO] {project}: {len(flist)} file(s), {proj_df.shape[0]} rows, {proj_df.shape[1]} cols")

    if not per_project_frames:
        print("[ERROR] No project dataframes assembled. Nothing to write.", file=sys.stderr)
        sys.exit(3)

    # ----- Build ALL dataframe -----
    all_frames = list(per_project_frames.values())
    df_all = pd.concat(all_frames, ignore_index=True)
    if args.add_project_col and "project_id" not in df_all.columns:
        # (Shouldn't happen, but just in case)
        df_all.insert(0, "project_id", None)

    # ----- SUMMARY dataframe -----
    df_summary = pd.DataFrame(project_counts).sort_values(["n_rows", "project_id"], ascending=[False, True])
    df_summary["total_rows_all_projects"] = df_all.shape[0]
    df_summary["union_columns"] = df_all.shape[1]

    # Optional coverage: % non-null per key columns (example metrics)
    key_cols = ["run_accession", "sample_accession", "study_accession", "experiment_accession"]
    for kc in key_cols:
        if kc in df_all.columns:
            # coverage per project
            coverage = []
            for project, pdf in per_project_frames.items():
                denom = len(pdf)
                num = int(pdf[kc].notna().sum()) if kc in pdf.columns else 0
                pct = (100.0 * num / denom) if denom else 0.0
                coverage.append((project, pct))
            # put into summary as columns like coverage_run_accession_pct
            cdict = dict(coverage)
            df_summary[f"coverage_{kc}_pct"] = df_summary["project_id"].map(cdict)

    # ----- Write Excel workbook -----
    out_xlsx = Path(args.out_xlsx).expanduser().resolve()
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    used_names = set()

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xlw:
        # SUMMARY first (quick glance)
        df_summary.to_excel(xlw, index=False, sheet_name=safe_sheet_name("SUMMARY", used_names))

        # ALL rows
        # Excel row limits: ~1,048,576 rows. If you exceed, consider splitting or writing CSV only.
        df_all.to_excel(xlw, index=False, sheet_name=safe_sheet_name("ALL", used_names))

        # Per-project sheets
        for project, pdf in per_project_frames.items():
            sheet = safe_sheet_name(project, used_names)
            pdf.to_excel(xlw, index=False, sheet_name=sheet)

    print(f"[OK] Excel workbook written: {out_xlsx}")

    # ----- Optional merged CSV/TSV -----
    if args.out_csv:
        out_csv = Path(args.out_csv).expanduser().resolve()
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        # Infer delimiter from extension
        ext = out_csv.suffix.lower()
        if ext in [".tsv", ".txt"]:
            df_all.to_csv(out_csv, index=False, sep="\t")
        else:
            df_all.to_csv(out_csv, index=False)
        print(f"[OK] Merged ALL rows written: {out_csv}")


if __name__ == "__main__":
    # Make pandas not silently downcast
    pd.options.mode.copy_on_write = True
    main()
