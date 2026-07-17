#!/usr/bin/env bash
set -euo pipefail

# download_quercus_biosample_metadata.sh
#
# Input expected from previous script:
#   quercus_ncbi_metadata/quercus_sra_runinfo_with_species.csv
#
# Outputs:
#   quercus_ncbi_metadata/biosample/xml/*.xml
#   quercus_ncbi_metadata/biosample/quercus_biosample_attributes_long.csv
#   quercus_ncbi_metadata/biosample/quercus_biosample_attributes_wide.csv
#   quercus_ncbi_metadata/quercus_sra_runinfo_with_species_and_biosample.csv
#
# Optional, recommended by NCBI:
#   export NCBI_EMAIL="your.email@example.com"
#   export NCBI_API_KEY="your_ncbi_api_key"   # optional

OUTDIR="data/ncbi_metadata"
RUNINFO_CSV="${OUTDIR}/quercus_sra_runinfo_with_species.csv"

BIOSAMPLE_DIR="${OUTDIR}/biosample"
XML_DIR="${BIOSAMPLE_DIR}/xml"

mkdir -p "$BIOSAMPLE_DIR" "$XML_DIR"

if [[ ! -s "$RUNINFO_CSV" ]]; then
    echo "ERROR: Cannot find RunInfo metadata:"
    echo "  $RUNINFO_CSV"
    echo "Run the SRA RunInfo metadata script first."
    exit 1
fi

python3 - <<'PY'
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

OUTDIR = Path("data/ncbi_metadata")
RUNINFO_CSV = OUTDIR / "quercus_sra_runinfo_with_species.csv"

BIOSAMPLE_DIR = OUTDIR / "biosample"
XML_DIR = BIOSAMPLE_DIR / "xml"

LONG_CSV = BIOSAMPLE_DIR / "quercus_biosample_attributes_long.csv"
WIDE_CSV = BIOSAMPLE_DIR / "quercus_biosample_attributes_wide.csv"
MERGED_CSV = OUTDIR / "quercus_sra_runinfo_with_species_and_biosample.csv"
UID_MAP_TSV = BIOSAMPLE_DIR / "biosample_uid_map.tsv"
MISSING_TXT = BIOSAMPLE_DIR / "missing_biosamples.txt"

XML_DIR.mkdir(parents=True, exist_ok=True)

EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
TOOL = "quercus_biosample_downloader"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def make_url(endpoint, params):
    params = dict(params)
    params["tool"] = TOOL
    if EMAIL:
        params["email"] = EMAIL
    if API_KEY:
        params["api_key"] = API_KEY
    return f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)

def fetch_url(url, binary=True, retries=5, sleep_base=2):
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": f"{TOOL}/1.0"
                }
            )
            with urllib.request.urlopen(req, timeout=90) as response:
                data = response.read()
            # Be polite to NCBI. Without an API key, stay below 3 requests/sec.
            time.sleep(0.40)
            return data if binary else data.decode("utf-8", errors="replace")
        except Exception as e:
            last_error = e
            if attempt < retries:
                time.sleep(sleep_base * attempt)
            else:
                raise RuntimeError(f"Failed after {retries} attempts: {url}\n{e}") from e
    raise last_error

def clean_colname(x):
    x = (x or "").strip()
    x = re.sub(r"\s+", "_", x)
    x = re.sub(r"[^A-Za-z0-9_.-]", "_", x)
    x = re.sub(r"_+", "_", x)
    x = x.strip("_")
    return x or "unknown"

def get_text(elem, path):
    found = elem.find(path)
    if found is None:
        return ""
    return " ".join("".join(found.itertext()).split())

def find_biosample_elem(root):
    if root.tag == "BioSample":
        return root
    return root.find(".//BioSample")

def search_biosample_uid(accession):
    # First try an accession-specific query.
    queries = [
        f"{accession}[Accession]",
        accession
    ]

    for q in queries:
        url = make_url(
            "esearch.fcgi",
            {
                "db": "biosample",
                "term": q,
                "retmode": "json",
                "retmax": "5"
            }
        )
        txt = fetch_url(url, binary=False)
        data = json.loads(txt)
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            return ids[0]

    return ""

def download_biosample_xml(accession, uid):
    xml_path = XML_DIR / f"{accession}.xml"

    if xml_path.exists() and xml_path.stat().st_size > 0:
        return xml_path

    url = make_url(
        "efetch.fcgi",
        {
            "db": "biosample",
            "id": uid,
            "retmode": "xml"
        }
    )
    xml_bytes = fetch_url(url, binary=True)

    if not xml_bytes.strip():
        raise RuntimeError(f"Empty XML returned for {accession}")

    xml_path.write_bytes(xml_bytes)
    return xml_path

def parse_biosample_xml(accession, uid, xml_path):
    xml_text = xml_path.read_text(errors="replace")
    root = ET.fromstring(xml_text)
    bs = find_biosample_elem(root)

    if bs is None:
        raise RuntimeError(f"No BioSample element found in {xml_path}")

    row = {
        "BioSample": accession,
        "biosample_uid": uid,
        "xml_accession": bs.attrib.get("accession", ""),
        "xml_id": bs.attrib.get("id", ""),
        "access": bs.attrib.get("access", ""),
        "publication_date": bs.attrib.get("publication_date", ""),
        "submission_date": bs.attrib.get("submission_date", ""),
        "last_update": bs.attrib.get("last_update", ""),
        "sample_title": get_text(bs, ".//Description/Title"),
        "owner": get_text(bs, ".//Owner/Name"),
        "model": get_text(bs, ".//Models/Model"),
    }

    org = bs.find(".//Description/Organism")
    if org is not None:
        row["organism"] = org.attrib.get("taxonomy_name", "") or " ".join("".join(org.itertext()).split())
        row["tax_id"] = org.attrib.get("taxonomy_id", "")
    else:
        row["organism"] = ""
        row["tax_id"] = ""

    # IDs such as BioProject, SRA, Sample name, etc., when present.
    for id_elem in bs.findall(".//Ids/Id"):
        db = clean_colname(id_elem.attrib.get("db", "id")).lower()
        value = " ".join("".join(id_elem.itertext()).split())
        if value:
            key = f"id_{db}"
            if key in row and row[key]:
                row[key] += f";{value}"
            else:
                row[key] = value

    long_rows = []

    for attr in bs.findall(".//Attributes/Attribute"):
        raw_name = attr.attrib.get("attribute_name", "") or attr.attrib.get("harmonized_name", "") or attr.attrib.get("display_name", "")
        harmonized_name = attr.attrib.get("harmonized_name", "")
        display_name = attr.attrib.get("display_name", "")
        value = " ".join("".join(attr.itertext()).split())

        if not raw_name and not value:
            continue

        attr_col = "attr_" + clean_colname(harmonized_name or raw_name or display_name)

        if attr_col in row and row[attr_col]:
            row[attr_col] += f";{value}"
        else:
            row[attr_col] = value

        long_rows.append({
            "BioSample": accession,
            "biosample_uid": uid,
            "attribute_name": raw_name,
            "harmonized_name": harmonized_name,
            "display_name": display_name,
            "value": value
        })

    return row, long_rows

# ---------------------------------------------------------------------
# 1. Read RunInfo and extract BioSample accessions
# ---------------------------------------------------------------------

with RUNINFO_CSV.open(newline="") as f:
    reader = csv.DictReader(f)
    run_rows = list(reader)
    run_fields = reader.fieldnames or []

if "BioSample" not in run_fields:
    raise SystemExit("ERROR: RunInfo CSV does not contain a 'BioSample' column.")

biosamples = sorted({
    row.get("BioSample", "").strip()
    for row in run_rows
    if row.get("BioSample", "").strip()
})

if not biosamples:
    raise SystemExit("ERROR: No BioSample accessions found in RunInfo CSV.")

print(f"Unique BioSample accessions found: {len(biosamples)}")

# ---------------------------------------------------------------------
# 2. Download and parse BioSample XML
# ---------------------------------------------------------------------

wide_rows = []
long_rows_all = []
uid_records = []
missing = []

for i, accession in enumerate(biosamples, start=1):
    print(f"[{i}/{len(biosamples)}] BioSample {accession}")

    try:
        uid = search_biosample_uid(accession)
        if not uid:
            print(f"  WARNING: No BioSample UID found for {accession}", file=sys.stderr)
            missing.append(accession)
            continue

        xml_path = download_biosample_xml(accession, uid)
        wide_row, long_rows = parse_biosample_xml(accession, uid, xml_path)

        wide_rows.append(wide_row)
        long_rows_all.extend(long_rows)
        uid_records.append((accession, uid))

    except Exception as e:
        print(f"  ERROR: {accession}: {e}", file=sys.stderr)
        missing.append(accession)

# ---------------------------------------------------------------------
# 3. Write UID map
# ---------------------------------------------------------------------

with UID_MAP_TSV.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["BioSample", "biosample_uid"])
    writer.writerows(uid_records)

# ---------------------------------------------------------------------
# 4. Write long-format BioSample attributes
# ---------------------------------------------------------------------

long_fields = [
    "BioSample",
    "biosample_uid",
    "attribute_name",
    "harmonized_name",
    "display_name",
    "value"
]

with LONG_CSV.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=long_fields)
    writer.writeheader()
    for row in long_rows_all:
        writer.writerow(row)

# ---------------------------------------------------------------------
# 5. Write wide-format BioSample table
# ---------------------------------------------------------------------

base_fields = [
    "BioSample",
    "biosample_uid",
    "xml_accession",
    "xml_id",
    "access",
    "publication_date",
    "submission_date",
    "last_update",
    "sample_title",
    "owner",
    "model",
    "organism",
    "tax_id"
]

all_wide_fields = set()
for row in wide_rows:
    all_wide_fields.update(row.keys())

extra_fields = sorted(all_wide_fields - set(base_fields))
wide_fields = [x for x in base_fields if x in all_wide_fields] + extra_fields

with WIDE_CSV.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=wide_fields, extrasaction="ignore")
    writer.writeheader()
    for row in wide_rows:
        writer.writerow(row)

# ---------------------------------------------------------------------
# 6. Merge BioSample wide metadata back into RunInfo table
# ---------------------------------------------------------------------

biosample_by_acc = {
    row["BioSample"]: row
    for row in wide_rows
    if row.get("BioSample")
}

biosample_fields_to_add = [
    x for x in wide_fields
    if x != "BioSample"
]

merged_fields = run_fields + [f"bs_{x}" for x in biosample_fields_to_add]

with MERGED_CSV.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=merged_fields, extrasaction="ignore")
    writer.writeheader()

    for run_row in run_rows:
        out = dict(run_row)
        bs_acc = run_row.get("BioSample", "").strip()
        bs_row = biosample_by_acc.get(bs_acc, {})

        for field in biosample_fields_to_add:
            out[f"bs_{field}"] = bs_row.get(field, "")

        writer.writerow(out)

# ---------------------------------------------------------------------
# 7. Missing list
# ---------------------------------------------------------------------

with MISSING_TXT.open("w") as f:
    for accession in missing:
        f.write(accession + "\n")

print()
print("Done.")
print(f"BioSample XML files:       {XML_DIR}")
print(f"Long BioSample CSV:        {LONG_CSV}")
print(f"Wide BioSample CSV:        {WIDE_CSV}")
print(f"Merged RunInfo+BioSample:  {MERGED_CSV}")
print(f"Missing BioSamples:        {len(missing)}")
if missing:
    print(f"Missing list:              {MISSING_TXT}")
PY
