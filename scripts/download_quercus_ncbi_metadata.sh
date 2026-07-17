#!/usr/bin/env bash
set -euo pipefail

# download_quercus_ncbi_metadata.sh
#
# Output:
#   quercus_ncbi_metadata/
#     accessions.tsv
#     accessions.txt
#     quercus_sra_runinfo_raw.csv
#     quercus_sra_runinfo_with_species.csv
#     raw_parts/

OUTDIR="data/ncbi_metadata"
RAW_DIR="${OUTDIR}/raw_parts"
BATCH_SIZE=20

mkdir -p "$OUTDIR" "$RAW_DIR"

cat > "${OUTDIR}/accessions.tsv" <<'EOF'
Quercus_robur	ERR12033543
Quercus_robur	ERR12033542
Quercus_robur	ERR12033541
Quercus_robur	ERR12033540
Quercus_robur	ERR12033538
Quercus_petraea	ERR12400190
Quercus_petraea	ERR12400189
Quercus_petraea	ERR12400188
Quercus_petraea	ERR12405606
Quercus_petraea	ERR12400284
Quercus_acutissima	SRR16661031
Quercus_acutissima	SRR16248089
Quercus_acutissima	SRR16248111
Quercus_acutissima	SRR16248130
Quercus_acutissima	SRR16248136
Quercus_variabilis	SRR26589715
Quercus_variabilis	SRR26589716
Quercus_variabilis	SRR26589718
Quercus_variabilis	SRR26589720
Quercus_variabilis	SRR26589721
Quercus_lobata	SRR14632920
Quercus_lobata	SRR14632922
Quercus_lobata	SRR14632924
Quercus_lobata	SRR14632927
Quercus_lobata	SRR14632929
Quercus_fabri	SRR16661037
Quercus_fabri	SRR16661048
Quercus_rubra	SRR23696573
Quercus_rubra	SRR16661059
Quercus_rubra	SRR16661060
Quercus_rubra	SRR16661095
Quercus_rubra	ERR10015092
Quercus_rubra	SRR23696575
Quercus_rubra	SRR27110051
Quercus_rubra	SRR32267770
Quercus_rubra	SRR32267808
Quercus_rubra	SRR32267797
Quercus_guyavifolia	SRR16661119
Quercus_guyavifolia	SRR11398771
Quercus_guyavifolia	SRR11398778
Quercus_guyavifolia	SRR11398764
Quercus_guyavifolia	SRR11398769
Quercus_dentata	SRR16996306
Quercus_dentata	SRR16996295
Quercus_dentata	SRR16996303
Quercus_dentata	SRR16996319
Quercus_dentata	SRR16996281
Quercus_aquifolioides	SRR13569288
Quercus_aquifolioides	SRR16661008
Quercus_aquifolioides	SRR16661007
Quercus_aquifolioides	SRR11398804
Quercus_aquifolioides	SRR11398807
EOF

cut -f2 "${OUTDIR}/accessions.tsv" > "${OUTDIR}/accessions.txt"

RAW_COMBINED="${OUTDIR}/quercus_sra_runinfo_raw.csv"
FINAL_CSV="${OUTDIR}/quercus_sra_runinfo_with_species.csv"

rm -f "$RAW_COMBINED"

echo "Downloading NCBI SRA RunInfo metadata..."

split -l "$BATCH_SIZE" "${OUTDIR}/accessions.txt" "${RAW_DIR}/batch_"

part_id=0
for batch_file in "${RAW_DIR}"/batch_*; do
    part_id=$((part_id + 1))

    acc_csv=$(paste -sd, "$batch_file")
    part_csv="${RAW_DIR}/runinfo_part_${part_id}.csv"

    echo "Batch ${part_id}: ${acc_csv}"

    curl -L --retry 5 --retry-delay 5 --fail \
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${acc_csv}" \
        -o "$part_csv"

    if [[ ! -s "$part_csv" ]]; then
        echo "ERROR: empty metadata file for batch ${part_id}"
        exit 1
    fi

    if [[ ! -f "$RAW_COMBINED" ]]; then
        cat "$part_csv" > "$RAW_COMBINED"
    else
        tail -n +2 "$part_csv" >> "$RAW_COMBINED"
    fi

    # Be polite to NCBI
    sleep 1
done

echo "Adding requested species labels..."

python3 - <<'PY'
import csv
from pathlib import Path

outdir = Path("data/ncbi_metadata")
map_file = outdir / "accessions.tsv"
raw_csv = outdir / "quercus_sra_runinfo_raw.csv"
final_csv = outdir / "quercus_sra_runinfo_with_species.csv"

species_by_run = {}

with map_file.open(newline="") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        species, run = line.split("\t")
        species_by_run[run] = species

with raw_csv.open(newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    old_fields = reader.fieldnames or []

if not rows:
    raise SystemExit("ERROR: no metadata rows were downloaded.")

new_fields = ["requested_species"] + old_fields

with final_csv.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=new_fields)
    writer.writeheader()

    for row in rows:
        run = row.get("Run", "")
        row_out = {"requested_species": species_by_run.get(run, "UNKNOWN")}
        row_out.update(row)
        writer.writerow(row_out)

downloaded = {row.get("Run", "") for row in rows}
requested = set(species_by_run)

missing = sorted(requested - downloaded)

print(f"Downloaded metadata rows: {len(rows)}")
print(f"Requested accessions:      {len(requested)}")
print(f"Final file:                {final_csv}")

if missing:
    print("\nWARNING: These accessions were requested but not found in downloaded RunInfo:")
    for acc in missing:
        print(acc)
PY

echo
echo "Done."
echo "Raw metadata:"
echo "  ${RAW_COMBINED}"
echo "Metadata with species labels:"
echo "  ${FINAL_CSV}"
