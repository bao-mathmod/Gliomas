#!/usr/bin/env bash
set -euo pipefail

# --------- Config ----------
BASE_DIR="/mnt/12T/chibao/data/cellranger_data"
OUT_SHEET="/mnt/12T/chibao/data/stuff_data/special/sample_metadata.tsv"
LOG="/mnt/12T/chibao/data/stuff_data/special/sample_metadata.log"
ERR="/mnt/12T/chibao/data/stuff_data/special/sample_metadata.err"
DEFAULT_REF="GRCh38"          # change if needed
CR_VERSION="9.0.1"            # fixed as per your runs
# ---------------------------

# Return the closest ancestor directory whose basename matches a regex
closest_ancestor_match() {
    local start="$1" regex="$2" max_up="${3:-10}"
    local cur="$start"
    for _ in $(seq 1 "$max_up"); do
        cur="$(dirname "$cur")"
        local base="$(basename "$cur")"
        if [[ "$base" =~ $regex ]]; then
            echo "$base"
            return 0
        fi
    done
    echo "NA"
    return 1
}

# Convenience wrappers for our four ID types
get_run_id()       { closest_ancestor_match "$1" '^SRR[0-9]+$' 12; }
get_experiment_id(){ closest_ancestor_match "$1" '^SRX[0-9]+$' 12; }
get_sample_id()    { closest_ancestor_match "$1" '^SAMN[0-9]+$' 12; }
get_project_id()   { closest_ancestor_match "$1" '^PRJNA[0-9]+$' 12; }

# Create header once
if [[ ! -f "$OUT_SHEET" ]]; then
    echo -e "Project_ID\tSample_ID\tExperiment_ID\tRun_ID\tCellRanger_Output_Dir\tCellRanger_Version\tReference_Genome\tEstimated_Number_of_Cells\tMean_Reads_per_Cell\tMedian_Genes_per_Cell\tTotal_Reads\tValid_Barcodes\tValid_UMIs\tSequencing_Saturation\tQ30_Bases_in_Barcode\tQ30_Bases_in_RNA_Read\tQ30_Bases_in_UMI\tReads_Mapped_to_Genome\tReads_Mapped_Confidently_to_Genome\tReads_Mapped_Confidently_to_Intergenic\tReads_Mapped_Confidently_to_Intronic\tReads_Mapped_Confidently_to_Exonic\tReads_Mapped_Confidently_to_Transcriptome\tReads_Mapped_Antisense_to_Gene\tFraction_Reads_in_Cells\tTotal_Genes_Detected\tMedian_UMI_Counts_per_Cell\tNotes" > "$OUT_SHEET"
fi

# CSV -> TSV metrics extractor (robust to quotes/commas/percent signs)
parse_metrics_csv() {
    local csv_path="$1"
    python3 - "$csv_path" <<'PY'
import csv, sys, re

csv_path = sys.argv[1]

# Canonical list (Cell Ranger labels). We will match keys case/space/punct-insensitively.
wanted = [
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median Genes per Cell",
    "Number of Reads",
    "Valid Barcodes",
    "Valid UMIs",
    "Sequencing Saturation",
    "Q30 Bases in Barcode",
    "Q30 Bases in RNA Read",
    "Q30 Bases in UMI",
    "Reads Mapped to Genome",
    "Reads Mapped Confidently to Genome",
    "Reads Mapped Confidently to Intergenic Regions",
    "Reads Mapped Confidently to Intronic Regions",
    "Reads Mapped Confidently to Exonic Regions",
    "Reads Mapped Confidently to Transcriptome",
    "Reads Mapped Antisense to Gene",
    "Fraction Reads in Cells",
    "Total Genes Detected",
    "Median UMI Counts per Cell",
]

def norm_key(s):
    # lower, remove non-alphanum, collapse spaces: makes matching robust
    return re.sub(r'[^a-z0-9]+', '_', s.lower()).strip('_')

wanted_norm = [norm_key(k) for k in wanted]

def clean_val(v):
    if v is None:
        return ""
    v = str(v).strip().strip('"').strip("'")
    if v == "" or v.lower() in {"na", "n/a", "nan"}:
        return ""
    # remove thousands separators inside numbers
    v = v.replace(",", "")
    # strip trailing % but keep the numeric value as 0-100 (easier to sort/filter)
    if v.endswith("%"):
        v = v[:-1]
    return v

try:
    with open(csv_path, newline='') as fh:
        reader = csv.reader(fh)
        rows = list(reader)
        if len(rows) < 2:
            print("", end="")  # empty -> caller can handle
            sys.exit(0)
        header = rows[0]
        values = rows[1]

        # build dict with normalized keys
        keymap = {norm_key(k): i for i, k in enumerate(header)}
        out_vals = []
        for k_norm in wanted_norm:
            if k_norm in keymap and keymap[k_norm] < len(values):
                out_vals.append(clean_val(values[keymap[k_norm]]))
            else:
                out_vals.append("")  # missing column -> blank cell

        # Output as a single TSV line to stdout
        print("\t".join(out_vals), end="")

except Exception as e:
    # On error, print nothing (bash side will log) and exit non-zero
    sys.stderr.write(f"[parse_error] {csv_path}: {e}\n")
    sys.exit(2)
PY
}

# Walk all metrics_summary.csv files
# Note: we handle spaces/newlines safely.
# ... (same header/config/functions as you posted)

while IFS= read -r -d '' METRICS; do
  OUTDIR="$(dirname "$METRICS")"

  RUN_ID="$(get_run_id "$OUTDIR")"
  EXP_ID="$(get_experiment_id "$OUTDIR")"
  SAMPLE_ID="$(get_sample_id "$OUTDIR")"
  PROJECT_ID="$(get_project_id "$OUTDIR")"

  if [[ "$EXP_ID" == "NA" ]]; then
    : # old layout .../count/<SRR>/outs → keep EXP_ID=NA
  fi
  if [[ "$SAMPLE_ID" == "NA" ]]; then
    : # add a custom rule here if you ever encode sample names differently
  fi

  [[ "$PROJECT_ID" == "NA" ]] && echo "$(date +'%F %T')  WARN missing PRJNA for: $OUTDIR" >> "$LOG"
  [[ "$RUN_ID" == "NA" ]] && echo "$(date +'%F %T')  WARN missing SRR for: $OUTDIR" >> "$LOG"

  if ! METRIC_LINE=$(parse_metrics_csv "$METRICS"); then
    echo "$(date +'%F %T')  ERROR parsing: $METRICS" >> "$ERR"
    continue
  fi
  if [[ -z "$METRIC_LINE" ]]; then
    echo "$(date +'%F %T')  WARNING empty metrics: $METRICS" >> "$LOG"
    continue
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t-\n" \
    "$PROJECT_ID" "$SAMPLE_ID" "$EXP_ID" "$RUN_ID" "$OUTDIR" "$CR_VERSION" "$DEFAULT_REF" "$METRIC_LINE" \
    >> "$OUT_SHEET"

done < <(find "$BASE_DIR" -type f -name "metrics_summary.csv" -print0)
