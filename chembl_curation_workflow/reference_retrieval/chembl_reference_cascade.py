#!/usr/bin/env python3
"""
Resolve a single primary reference identifier for each record from its
ChEMBL document ID, using a defined precedence hierarchy.

Background
----------
DOI, PubMed ID, and patent ID are all attributes of the same ChEMBL
`document` resource. They are therefore retrieved in a single pass, and
the precedence hierarchy is applied afterwards as a local selection step.
This reproduces the outcome of running three sequential retrieval scripts
while issuing roughly one third of the API requests, and it makes the
selection rule explicit rather than implicit in execution order.

Row retention
-------------
Each input row is an individual bioactivity measurement. A single source
document, and therefore a single DOI, PubMed ID, or patent ID, commonly
applies to many rows, since one publication may report measurements for
many compounds against several targets. Identifiers are retrieved once
per unique document ID and written back to every row referencing that
document.

Row count, row order, and the document ID column itself are unchanged
between input and output. All three are enforced as explicit checks
before the output is written, so the output can be placed alongside the
input column for column.

Precedence
----------
1. DOI          present for essentially all journal-published sources
2. PubMed ID    assigned where available; coverage is incomplete for
                non-MEDLINE-indexed and older literature
3. Patent ID    assigned only where neither of the above resolved, which
                corresponds to bioactivity data sourced from patents

The order is configurable via --order so the rule is visible and testable
rather than hard-coded.

Output columns
--------------
doi, pubmed_id, patent_id   raw values as returned by ChEMBL, all
                            populated wherever they exist
reference_id                the identifier selected under the hierarchy
reference_type              doi | pubmed_id | patent_id | unresolved

Keeping the raw fields alongside the selection allows a reader to
re-derive the second from the first and to audit the rule independently.

Usage
-----
Open a terminal in the folder containing the data and run:

    python chembl_reference_cascade.py

With no arguments the script reads `chembl_doc_id_input` from the
current folder (any of .xlsx, .xls, .csv, .tsv, .txt), detects the
document ID column automatically, and writes
`chembl_api_reference_cascade_output.xlsx` alongside it. File paths are
never required. To override any of these:

    python chembl_reference_cascade.py --in other_file.xlsx \
        --out results.xlsx --id-col "Document ChEMBL ID"
"""

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

# Default input file stem. The script looks for this name in the current
# working directory with any supported extension, so it can be run with
# no arguments at all from the folder holding the data.
DEFAULT_INPUT_STEM = "chembl_doc_id_input"

# Default output filename. Written to the current working directory.
DEFAULT_OUTPUT_NAME = "chembl_api_reference_cascade_output.xlsx"

# Extensions tried against the default stem, in order of preference.
SUPPORTED_EXTENSIONS = [".xlsx", ".xls", ".csv", ".tsv", ".txt"]

# Column names tried when --id-col is not given. Earlier retrieval
# scripts used different spellings for the same field, so all known
# variants are recognized.
ID_COLUMN_CANDIDATES = [
    "chembl_document_id",
    "document_chembl_id",
    "Document ChEMBL ID",
    "Document Chembl ID",
    "ChEMBL Document ID",
    "docid",
]

# Precedence order applied when selecting the primary reference.
DEFAULT_ORDER = ["doi", "pubmed_id", "patent_id"]

# Additional bibliographic fields retrieved for provenance and for
# manual verification of the resolved references.
CONTEXT_FIELDS = ["journal", "year", "title", "doc_type"]

MULTI_VALUE_SEP = "|"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resolve primary reference identifiers from ChEMBL "
                    "document IDs using a precedence hierarchy."
    )
    parser.add_argument("--in", dest="infile", default=None,
                        help="Input file in the current folder. Defaults to "
                             f"'{DEFAULT_INPUT_STEM}' with any supported "
                             "extension. A path is never required.")
    parser.add_argument("--out", dest="outfile", default=None,
                        help="Output file, written to the current folder. "
                             f"Defaults to '{DEFAULT_OUTPUT_NAME}'. If the "
                             "extension is .xlsx or .xls, output is written "
                             "as Excel; otherwise as CSV.")
    parser.add_argument("--id-col", default=None,
                        help="Column holding ChEMBL document IDs. Detected "
                             "automatically from known variants if omitted.")
    parser.add_argument("--order", nargs="+", default=DEFAULT_ORDER,
                        choices=["doi", "pubmed_id", "patent_id"],
                        help="Precedence order, highest priority first "
                             f"(default: {' '.join(DEFAULT_ORDER)}).")
    parser.add_argument("--no-context", action="store_true",
                        help="Omit journal, year, title, and doc_type "
                             "columns from the output.")
    parser.add_argument("--batch-size", type=int, default=50,
                        help="Document IDs per API request (default: 50).")
    parser.add_argument("--retries", type=int, default=3,
                        help="Retry attempts per batch (default: 3).")
    parser.add_argument("--retry-wait", type=float, default=5.0,
                        help="Seconds between retries, doubled each "
                             "attempt (default: 5).")
    return parser.parse_args()


def read_table(path):
    """Read CSV, TSV, or Excel based on file extension."""
    suffix = Path(path).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if suffix in (".csv", ".txt", ".tsv"):
        return pd.read_csv(path, sep="\t" if suffix == ".tsv" else ",")
    raise ValueError(f"Unsupported input format: {suffix}")


def resolve_input(explicit):
    """
    Locate the input file in the current working directory.

    Everything is resolved relative to the folder the script is run from,
    so a path is never needed. If --in is given without an extension, the
    supported extensions are tried in turn.
    """
    cwd = Path.cwd()

    if explicit:
        candidate = Path(explicit)
        # Ignore any directory component; work in the current folder.
        candidate = cwd / candidate.name
        if candidate.exists():
            return candidate
        if not candidate.suffix:
            for ext in SUPPORTED_EXTENSIONS:
                with_ext = candidate.with_suffix(ext)
                if with_ext.exists():
                    return with_ext
        sys.exit(f"Error: '{candidate.name}' not found in {cwd}.\n"
                 f"{describe_available_files(cwd)}")

    for ext in SUPPORTED_EXTENSIONS:
        candidate = cwd / f"{DEFAULT_INPUT_STEM}{ext}"
        if candidate.exists():
            return candidate

    sys.exit(f"Error: no input file given and no '{DEFAULT_INPUT_STEM}' "
             f"file found in {cwd}.\n"
             f"{describe_available_files(cwd)}\n"
             f"Run with --in <filename> to choose one.")


def describe_available_files(folder):
    """List data files in the folder to help the user pick one."""
    found = sorted(
        p.name for p in folder.iterdir()
        if p.is_file() and p.suffix.lower() in SUPPORTED_EXTENSIONS
    )
    if not found:
        return "No CSV or Excel files were found in this folder."
    listing = "\n  ".join(found[:20])
    more = f"\n  ... and {len(found) - 20} more" if len(found) > 20 else ""
    return f"Files available here:\n  {listing}{more}"


def resolve_output(explicit, infile):
    """Derive the output name from the input unless one is given.

    Default output is the fixed name in DEFAULT_OUTPUT_NAME. If --out is
    provided, it is used as given (in the current folder).
    """
    if explicit:
        return Path.cwd() / Path(explicit).name
    return Path.cwd() / DEFAULT_OUTPUT_NAME


def resolve_id_column(df, explicit):
    """
    Identify the ChEMBL document ID column.

    Matching is case-insensitive and ignores spaces and underscores, so
    the various spellings used across earlier scripts all resolve.
    """
    if explicit:
        if explicit in df.columns:
            return explicit
        sys.exit(f"Error: column '{explicit}' not found.\n"
                 f"Columns present: {', '.join(map(str, df.columns))}")

    def normalize(name):
        return str(name).lower().replace(" ", "").replace("_", "")

    normalized = {normalize(c): c for c in df.columns}
    for candidate in ID_COLUMN_CANDIDATES:
        match = normalized.get(normalize(candidate))
        if match:
            return match

    sys.exit("Error: could not identify the ChEMBL document ID column.\n"
             f"Columns present: {', '.join(map(str, df.columns))}\n"
             "Specify it with --id-col \"<column name>\".")


def collect_unique_ids(series):
    """
    Flatten a column that may contain comma-separated ID lists into an
    ordered set of unique, stripped ChEMBL document IDs.
    """
    unique, seen = [], set()
    for cell in series:
        if pd.isna(cell):
            continue
        for token in str(cell).split(","):
            token = token.strip()
            if token and token not in seen:
                seen.add(token)
                unique.append(token)
    return unique


def fetch_batch(doc_ids, fields, retries, retry_wait):
    """
    Retrieve one batch of documents, keyed by document ID.

    Raises on persistent failure. A transient API error must never be
    recorded as an absent identifier, since that is indistinguishable
    from a genuine absence and would silently corrupt the cascade.
    """
    requested = list(dict.fromkeys(fields + ["document_chembl_id"]))
    wait, last_error = retry_wait, None

    for attempt in range(1, retries + 1):
        try:
            records = new_client.document.filter(
                document_chembl_id__in=doc_ids
            ).only(requested)
            return {
                rec["document_chembl_id"]: rec
                for rec in records
                if rec.get("document_chembl_id")
            }
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                tqdm.write(f"  batch failed (attempt {attempt}/{retries}): "
                           f"{exc}; retrying in {wait:.0f}s")
                time.sleep(wait)
                wait *= 2

    raise RuntimeError(
        f"Batch of {len(doc_ids)} document IDs failed after {retries} "
        f"attempts. Last error: {last_error}"
    ) from last_error


def clean_value(value):
    """Normalize a field value to a stripped string; None becomes empty."""
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in ("none", "nan") else text


def select_reference(values, order):
    """
    Apply the precedence hierarchy to one record's retrieved values.

    Returns (reference_id, reference_type). A record with no identifier
    of any kind is reported as 'unresolved' rather than left blank, so
    that these rows can be counted and inspected.
    """
    for field in order:
        value = values.get(field, "")
        if value:
            return value, field
    return "", "unresolved"


def assert_row_count(expected, actual, stage):
    """
    Enforce the row-count invariant.

    A ChEMBL document is a source, not a measurement. A single paper or
    patent commonly reports many bioactivity values, so one document ID
    legitimately recurs across many rows. Deduplication is applied only
    to the set of API requests; it must never reach the output. Any
    change in row count means identifiers have been collapsed onto the
    wrong key and the result is invalid.
    """
    if expected != actual:
        sys.exit(
            f"Error: row count changed at stage '{stage}': "
            f"{expected} in, {actual} out. Retrieval must not add or "
            f"remove rows. Records sharing a document identifier are "
            f"distinct measurements and are retained individually."
        )


def get_api_status():
    """Record the ChEMBL release so retrieval can be reproduced."""
    try:
        status = new_client.status.get()
        if isinstance(status, list):
            status = status[0]
        return {
            "chembl_db_version": status.get("chembl_db_version"),
            "chembl_release_date": status.get("chembl_release_date"),
            "api_version": status.get("api_version"),
        }
    except Exception as exc:
        return {"error": f"Could not retrieve ChEMBL status: {exc}"}


def main():
    args = parse_args()

    fields = list(args.order)
    if not args.no_context:
        fields += CONTEXT_FIELDS

    infile = resolve_input(args.infile)
    outfile = resolve_output(args.outfile, infile)

    print(f"Working folder: {Path.cwd()}")
    print(f"Input:  {infile.name}")
    print(f"Output: {outfile.name}")

    df = read_table(infile)
    id_col = resolve_id_column(df, args.id_col)
    print(f"Document ID column: {id_col}")

    unique_ids = collect_unique_ids(df[id_col])
    if not unique_ids:
        sys.exit(f"Error: no ChEMBL document IDs found in '{id_col}'.")

    print(f"Resolving {len(unique_ids)} unique document IDs from "
          f"{len(df)} input rows.")
    print(f"Precedence: {' > '.join(args.order)}")

    # Single retrieval pass. Each unique ID is fetched exactly once,
    # however many rows reference it.
    lookup = {}
    batches = [unique_ids[i:i + args.batch_size]
               for i in range(0, len(unique_ids), args.batch_size)]
    for batch in tqdm(batches, desc="Querying ChEMBL"):
        lookup.update(fetch_batch(batch, fields, args.retries,
                                  args.retry_wait))

    unresolved_ids = [i for i in unique_ids if i not in lookup]
    if unresolved_ids:
        print(f"Warning: {len(unresolved_ids)} document ID(s) returned no "
              f"record and are reported as empty, not as missing data. "
              f"First few: {', '.join(unresolved_ids[:5])}")

    # Map retrieved values back onto the original rows, then apply the
    # precedence hierarchy per row.
    out = df.copy()
    per_field = {field: [] for field in fields}
    reference_ids, reference_types = [], []

    for cell in df[id_col]:
        tokens = ([t.strip() for t in str(cell).split(",") if t.strip()]
                  if not pd.isna(cell) else [])

        # Collect values for each field across all IDs in the cell.
        row_values = {}
        for field in fields:
            collected = []
            for token in tokens:
                record = lookup.get(token)
                value = clean_value(record.get(field)) if record else ""
                if value:
                    collected.append(value)
            joined = MULTI_VALUE_SEP.join(dict.fromkeys(collected))
            row_values[field] = joined
            per_field[field].append(joined)

        ref_id, ref_type = select_reference(row_values, args.order)
        reference_ids.append(ref_id)
        reference_types.append(ref_type)

    for field in fields:
        out[field] = per_field[field]
    out["reference_id"] = reference_ids
    out["reference_type"] = reference_types

    # Invariant: retrieval annotates rows, it never adds, removes, or
    # reorders them. Each input row is an individual bioactivity
    # measurement, whereas a document, PMID, DOI, or patent identifier may
    # legitimately apply to many rows. Deduplication is therefore valid
    # only when issuing API requests, never in the output. Values are
    # assigned positionally rather than merged on the document ID, since a
    # merge on a duplicated key can both reorder and multiply rows.
    if len(out) != len(df):
        sys.exit(f"Row count changed: {len(df)} in, {len(out)} out. "
                 "Retrieval must annotate rows, never add or remove them.")

    if not out.index.equals(df.index):
        sys.exit("Row order or index changed during retrieval. "
                 "Retrieval must preserve the original row order.")

    original = df[id_col].astype(str).tolist()
    final = out[id_col].astype(str).tolist()
    if original != final:
        first = next(i for i, (a, b) in enumerate(zip(original, final))
                     if a != b)
        sys.exit(f"Document ID column changed during retrieval. "
                 f"First difference at row {first}: "
                 f"'{original[first]}' became '{final[first]}'.")

    # Write the output in the format matching the output filename's
    # extension. XLSX is written via openpyxl; CSV via pandas.
    suffix = outfile.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        out.to_excel(outfile, index=False)
    else:
        out.to_csv(outfile, index=False)

    # Stage-by-stage counts mirror what the sequential scripts produced
    # and give a quick sanity check on the cascade.
    counts = out["reference_type"].value_counts().to_dict()
    rows_per_doc = (len(df) / len(unique_ids)) if unique_ids else 0
    print(f"\nWrote {len(out)} rows to {outfile.name} "
          f"(input rows: {len(df)}, unchanged)")
    print(f"{len(unique_ids)} unique document IDs across {len(df)} rows "
          f"({rows_per_doc:.1f} rows per document on average); "
          f"all rows retained.")
    print("Reference type breakdown:")
    for field in args.order + ["unresolved"]:
        n = counts.get(field, 0)
        pct = (100 * n / len(out)) if len(out) else 0
        print(f"  {field:<12} {n:>7}  ({pct:5.1f}%)")

    provenance = {
        "input_file": infile.name,
        "output_file": outfile.name,
        "id_column": id_col,
        "precedence_order": args.order,
        "fields_retrieved": fields,
        "input_rows": len(df),
        "output_rows": len(out),
        "rows_retained_unchanged": len(out) == len(df),
        "row_order_preserved": True,
        "unique_document_ids": len(unique_ids),
        "mean_rows_per_document": round(rows_per_doc, 3),
        "document_ids_with_no_record": len(unresolved_ids),
        "reference_type_counts": counts,
        "retrieval_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                 time.gmtime()),
        "chembl_status": get_api_status(),
    }
    provenance_path = outfile.with_suffix(outfile.suffix + ".provenance.json")
    with open(provenance_path, "w") as handle:
        json.dump(provenance, handle, indent=2)
    print(f"Wrote retrieval provenance to {provenance_path.name}")


if __name__ == "__main__":
    main()
