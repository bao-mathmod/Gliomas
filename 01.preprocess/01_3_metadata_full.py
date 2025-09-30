#!/usr/bin/env python3
import argparse
import os
import sys
import time
import traceback
from datetime import datetime
import pandas as pd
from glob import glob

from pysradb import SRAweb

KEY_FIELDS = ["disease", "age", "sex", "tissue", "cell_type", "treatment"]

def ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S%z")

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def init_logger(log_path):
    ensure_dir(os.path.dirname(log_path))
    def log(msg):
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"[{ts()}] {msg}\n")
        print(msg, flush=True)
    return log

def load_project_list(csv_path):
    df = pd.read_csv(csv_path)
    if "PRJNA" not in df.columns:
        raise ValueError(f"CSV must contain a 'PRJNA' column; got columns: {list(df.columns)}")
    prjnas = (
        df["PRJNA"].astype(str).str.strip().str.upper().dropna().unique().tolist()
    )
    prjnas = [p for p in prjnas if p.startswith("PRJN")]
    return prjnas, df

def read_existing(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            return pd.read_csv(path, sep="\t")
        except Exception:
            return pd.read_csv(path)
    return None

def write_tsv(df, path):
    ensure_dir(os.path.dirname(path))
    df.to_csv(path, sep="\t", index=False)

def append_summary_row(summary_path, row_dict):
    ensure_dir(os.path.dirname(summary_path))
    row_df = pd.DataFrame([row_dict])
    if not os.path.exists(summary_path):
        row_df.to_csv(summary_path, sep="\t", index=False)
    else:
        existing = read_existing(summary_path)
        if existing is None or existing.empty:
            row_df.to_csv(summary_path, sep="\t", index=False)
        else:
            combined = pd.concat([existing, row_df], ignore_index=True)
            dedup_cols = ["project", "status"]
            combined = combined.drop_duplicates(subset=dedup_cols, keep="last")
            combined.to_csv(summary_path, sep="\t", index=False)

def mark_processed(ledger_path, project, status, n_rows, message):
    ensure_dir(os.path.dirname(ledger_path))
    row = {
        "timestamp": ts(),
        "project": project,
        "status": status,
        "rows": n_rows,
        "message": message,
    }
    row_df = pd.DataFrame([row])
    if not os.path.exists(ledger_path):
        row_df.to_csv(ledger_path, sep="\t", index=False)
    else:
        existing = read_existing(ledger_path)
        if existing is None or existing.empty:
            row_df.to_csv(ledger_path, sep="\t", index=False)
        else:
            combined = pd.concat([existing, row_df], ignore_index=True)
            combined.to_csv(ledger_path, sep="\t", index=False)

def add_if_missing_columns(df, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = pd.NA
    return df

def summarize_project_df(project, df):
    col_candidates = {
        "run": ["run_accession", "run"],
        "sample": ["sample_accession", "sample"],
        "experiment": ["experiment_accession", "experiment"],
    }
    def pick(colnames, default=0):
        for c in colnames:
            if c in df.columns:
                return df[c].nunique()
        return default
    n_runs = pick(col_candidates["run"])
    n_samples = pick(col_candidates["sample"])
    n_experiments = pick(col_candidates["experiment"])
    coverage = {}
    for k in KEY_FIELDS:
        candidates = [c for c in df.columns if c.lower() == k]
        if candidates:
            c = candidates[0]
            coverage[f"{k}_filled"] = int(df[c].notna().sum())
            coverage[f"{k}_unique"] = int(df[c].dropna().nunique())
        else:
            coverage[f"{k}_filled"] = 0
            coverage[f"{k}_unique"] = 0
    row = {
        "timestamp": ts(),
        "project": project,
        "status": "success",
        "n_runs": int(n_runs),
        "n_samples": int(n_samples),
        "n_experiments": int(n_experiments),
        **coverage,
    }
    return row

def expand_sample_attributes(db, project, log):
    try:
        df = db.sra_metadata(project, detailed=True, expand_sample_attributes=True)
        return df
    except TypeError:
        log(f"{project}: pysradb missing expand_sample_attributes=True; falling back to manual expand.")
        df = db.sra_metadata(project, detailed=True)
        if "sample_attribute" in df.columns:
            sa = df["sample_attribute"].fillna("")
            parsed_rows = []
            for s in sa:
                entry = {}
                for kv in str(s).split("||"):
                    kv = kv.strip().strip("|").strip()
                    if not kv:
                        continue
                    if ":" in kv:
                        k, v = kv.split(":", 1)
                        entry[k.strip()] = v.strip()
                parsed_rows.append(entry)
            attr_df = pd.DataFrame(parsed_rows)
            df = pd.concat([df.drop(columns=["sample_attribute"]), attr_df], axis=1)
        return df

def merge_dataframes_align(dfs):
    """Align columns across many dataframes and concatenate."""
    all_cols = set()
    for d in dfs:
        all_cols.update(d.columns)
    all_cols = list(all_cols)
    aligned = []
    for d in dfs:
        missing = [c for c in all_cols if c not in d.columns]
        for c in missing:
            d[c] = pd.NA
        aligned.append(d[all_cols])
    return pd.concat(aligned, ignore_index=True)

def merge_to_master(master_path, new_df, log):
    existing = read_existing(master_path)
    if existing is None or existing.empty:
        write_tsv(new_df, master_path)
        log(f"master: created {master_path} with {len(new_df)} rows")
        return len(new_df)
    else:
        combined = merge_dataframes_align([existing, new_df])
        key = "run_accession" if "run_accession" in combined.columns else None
        if key:
            before = len(combined)
            combined = combined.drop_duplicates(subset=[key], keep="last")
            deduped = before - len(combined)
            log(f"master: deduped {deduped} rows by {key}")
        else:
            combined = combined.drop_duplicates(keep="last")
        write_tsv(combined, master_path)
        added = len(combined) - len(existing)
        log(f"master: updated {master_path}; +{max(0, added)} new rows")
        return max(0, added)

def already_processed(ledger_path, project):
    if not os.path.exists(ledger_path):
        return False
    df = read_existing(ledger_path)
    if df is None or df.empty:
        return False
    last = df[df["project"] == project].tail(1)
    if last.empty:
        return False
    return last["status"].iloc[0] == "success"

def merge_all_per_project(outdir, master_path, log):
    """Scan outdir/**/<PRJNA>_metadata.tsv and merge them into master_path."""
    pattern = os.path.join(outdir, "*", "*_metadata.tsv")
    files = sorted(glob(pattern))
    if not files:
        log(f"merge: no per-project metadata found at pattern {pattern}")
        return 0

    log(f"merge: found {len(files)} per-project TSVs to merge")
    dfs = []
    for fp in files:
        try:
            df = read_existing(fp)
            if df is not None and not df.empty:
                dfs.append(df)
        except Exception as e:
            log(f"merge: skip {fp} ({e})")

    if not dfs:
        log("merge: nothing to merge (all files empty/unreadable)")
        return 0

    merged = merge_dataframes_align(dfs)
    key = "run_accession" if "run_accession" in merged.columns else None
    if key:
        before = len(merged)
        merged = merged.drop_duplicates(subset=[key], keep="last")
        log(f"merge: deduped {before - len(merged)} rows by {key}")
    else:
        merged = merged.drop_duplicates(keep="last")

    write_tsv(merged, master_path)
    log(f"merge: wrote {master_path} with {len(merged)} rows")
    return len(merged)

def main():
    ap = argparse.ArgumentParser(description="Harvest SRA/ENA metadata for multiple PRJNA projects via pysradb.")
    ap.add_argument("--csv", default="/mnt/12T/chibao/data/official_data/glioma_accession.csv",
                    help="CSV path containing at least a PRJNA column (default: %(default)s)")
    ap.add_argument("--outdir", default="./out", help="Output directory (default: %(default)s)")
    ap.add_argument("--force", action="store_true", help="Force re-fetch even if project marked success")
    ap.add_argument("--skip-master", action="store_true", help="Skip updating master while harvesting")
    ap.add_argument("--merge-only", action="store_true",
                    help="Do not harvest; only merge all per-project TSVs into the master file")
    ap.add_argument("--merge-at-end", action="store_true",
                    help="After harvesting, rescan outdir and rebuild master from per-project files")
    ap.add_argument("--master-name", default="all_metadata.tsv",
                    help="Filename for merged master (default: %(default)s)")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    log_path = os.path.join(outdir, "logs", "metadata_harvest.log")
    ledger_path = os.path.join(outdir, "processed_projects.tsv")
    master_path = os.path.join(outdir, args.master_name)
    summary_path = os.path.join(outdir, "summary_projects.tsv")

    log = init_logger(log_path)

    if args.merge_only:
        log(f"=== MERGE-ONLY mode === outdir={outdir} master={master_path}")
        merged_rows = merge_all_per_project(outdir, master_path, log)
        log(f"=== MERGE-ONLY done; rows={merged_rows} ===")
        return

    # Normal harvest mode
    log(f"=== START metadata harvest === csv={args.csv} outdir={outdir} force={args.force}")
    try:
        projects, _ = load_project_list(args.csv)
    except Exception as e:
        log(f"FATAL: cannot read CSV: {e}")
        sys.exit(2)

    if not projects:
        log("No PRJNA projects found in CSV.")
        if args.merge_at_end:
            log("merge-at-end requested; running merge over existing per-project files.")
            merge_all_per_project(outdir, master_path, log)
        sys.exit(0)

    db = SRAweb()
    total_added_master = 0

    for i, project in enumerate(projects, 1):
        proj_dir = os.path.join(outdir, project)
        per_project_path = os.path.join(proj_dir, f"{project}_metadata.tsv")
        t0 = time.time()

        if (not args.force) and already_processed(ledger_path, project) and os.path.exists(per_project_path):
            log(f"[{i}/{len(projects)}] {project}: already processed — skipping (use --force to re-run)")
            continue

        try:
            log(f"[{i}/{len(projects)}] {project}: querying pysradb...")
            df = expand_sample_attributes(db, project, log)

            if df is None or df.empty:
                msg = "pysradb returned no rows"
                log(f"{project}: {msg}")
                mark_processed(ledger_path, project, "empty", 0, msg)
                append_summary_row(summary_path, {
                    "timestamp": ts(),
                    "project": project,
                    "status": "empty",
                    "n_runs": 0, "n_samples": 0, "n_experiments": 0,
                    **{f"{k}_filled": 0 for k in KEY_FIELDS},
                    **{f"{k}_unique": 0 for k in KEY_FIELDS},
                })
                continue

            for must in ["run_accession", "sample_accession", "experiment_accession", "study_accession"]:
                if must not in df.columns:
                    df[must] = pd.NA

            df = add_if_missing_columns(df, KEY_FIELDS)

            write_tsv(df, per_project_path)
            log(f"{project}: wrote per-project metadata: {per_project_path} ({len(df)} rows)")

            added = 0
            if not args.skip_master:
                added = merge_to_master(master_path, df, log)
                total_added_master += max(0, added)

            row = summarize_project_df(project, df)
            append_summary_row(summary_path, row)

            mark_processed(ledger_path, project, "success", len(df), f"added_to_master={added}")
            dt = time.time() - t0
            log(f"{project}: OK in {dt:.1f}s")

        except Exception as e:
            err = f"{type(e).__name__}: {e}"
            log(f"{project}: ERROR {err}")
            log(traceback.format_exc())
            mark_processed(ledger_path, project, "error", 0, err)
            continue

    if args.merge_at_end:
        log("merge-at-end: rebuilding master from all per-project TSVs.")
        merge_all_per_project(outdir, master_path, log)

    log(f"=== DONE; master located at {master_path} ===")

if __name__ == "__main__":
    main()
