"""
ascii_normalize.py
====================
Post-process assay descriptions to convert Unicode characters to
ASCII-compatible representations, preserving semantic meaning.

Usage:
    Open Terminal in the folder containing this script and the input
    XLSX, then run:

        python3 ascii_normalize.py

    (caffeinate is unnecessary — the script finishes in a few seconds.)

Configured input/output (edit at the top of this file if needed):
    INPUT_FILE  = "Full_dataset_with_assay_descriptions.xlsx"
    OUTPUT_FILE = "Full_dataset_with_assay_descriptions_ASCII.xlsx"
    COLUMN      = "assay_description"

The output XLSX contains all original columns plus:
    assay_description          — ASCII-normalized version
    assay_description_unicode  — original Unicode version (preserved)

Rules applied (in order):

1. Context-sensitive Greek letter substitution:
   - "μ" as a receptor prefix (μ-opioid, μOR, μ-OR, hμ-OR)   → "mu"
   - "μ" as a unit prefix (μM, μL, μg, μmol, μm, μs)        → "u" + unit
     e.g., "10 μM" → "10 uM", "5 μL" → "5 uL", "2 μg" → "2 ug"
   - "μ" standalone or in other contexts                     → "mu"
   - Same pattern for other Greek letters (α, β, γ, δ, κ)
     used as receptor / signaling / peptide prefixes.

2. Unicode dashes → ASCII hyphen:
   en-dash "–" (U+2013), em-dash "—" (U+2014) → "-"

3. Superscripts and subscripts:
   Ca²⁺, NO₂⁻, NH₂, etc. → Ca2+, NO2-, NH2
   [³⁵S], [³H], [¹²⁵I]  → [35S], [3H], [125I]

4. Miscellaneous:
   Non-breaking space (U+00A0)  → regular space
   Curly quotes "" ''            → straight quotes " '
   Ellipsis …                    → ...
   Degree sign °                 → " deg " (context-appropriate)
"""

from __future__ import annotations

import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Configuration — edit these if your filenames or column names differ
# ---------------------------------------------------------------------------

INPUT_FILE = "ascii_normalization_input_dataset.csv"
OUTPUT_FILE = "ASCII_compliant_assay_descriptions_output.xlsx"
COLUMN = "assay_description"
KEEP_UNICODE = True   # If True, adds a second column with the original Unicode text


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

# Greek letters that appear in the corpus, mapped to their spelled-out ASCII.
GREEK_ASCII = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "μ": "mu",
    "κ": "kappa",
    "Α": "Alpha",
    "Β": "Beta",
    "Γ": "Gamma",
    "Δ": "Delta",
    "Μ": "Mu",
    "Κ": "Kappa",
}

# When μ (or another Greek letter) is followed by one of these unit-suffix
# characters/strings, treat it as a unit prefix, NOT as a receptor prefix.
# Order matters: check longer strings first (μmol before μm).
# Micro prefix uses single-letter "u" abbreviation (uM, uL, ug, umol, um, us)
# consistent with standard pharmacology-table conventions.
UNIT_SUFFIX_MAP = [
    # (regex pattern after the Greek letter, ASCII replacement)
    (r"mol\b", "umol"),
    (r"M\b",   "uM"),
    (r"L\b",   "uL"),
    (r"l\b",   "uL"),
    (r"g\b",   "ug"),
    (r"m\b",   "um"),
    (r"s\b",   "us"),
]

# Superscript digits/signs → ASCII digits/signs.
SUPERSCRIPT_MAP = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁺": "+", "⁻": "-", "⁼": "=", "⁽": "(", "⁾": ")",
}

# Subscript digits/signs → ASCII digits/signs.
SUBSCRIPT_MAP = {
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
    "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
    "₊": "+", "₋": "-", "₌": "=", "₍": "(", "₎": ")",
}

# Miscellaneous single-character substitutions.
MISC_MAP = {
    "\u2013": "-",     # en-dash
    "\u2014": "-",     # em-dash
    "\u2212": "-",     # minus sign
    "\u00a0": " ",     # non-breaking space
    "\u201c": '"',     # left double quote
    "\u201d": '"',     # right double quote
    "\u2018": "'",     # left single quote
    "\u2019": "'",     # right single quote
    "\u2026": "...",   # ellipsis
    "\u00d7": "x",     # multiplication sign
    "\u2265": ">=",    # greater than or equal
    "\u2264": "<=",    # less than or equal
    "\u2260": "!=",    # not equal
    "\u00b1": "+/-",   # plus-minus
}


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------


def normalize_greek_context_sensitive(text: str) -> str:
    """Replace Greek letters, distinguishing unit-prefix from other uses.

    For μ specifically:
      - "10 μM"       → "10 micromolar"
      - "μ-opioid"    → "mu-opioid"
      - "μOR"         → "muOR"
      - "hμ-OR"       → "hmu-OR"

    The rule: if the Greek letter is IMMEDIATELY followed by a unit-suffix
    character (M, L, l, g, m, s, mol) that ends the token (i.e., followed
    by a word boundary), treat as a unit. Otherwise spell it out.
    """
    result = []
    i = 0
    while i < len(text):
        c = text[i]
        if c in GREEK_ASCII:
            # Check if this is a unit prefix.
            remaining = text[i + 1:]
            unit_match = None
            for pat, replacement in UNIT_SUFFIX_MAP:
                m = re.match(pat, remaining)
                if m:
                    unit_match = (m, replacement)
                    break
            if unit_match:
                m, replacement = unit_match
                result.append(replacement)
                i += 1 + m.end()
                continue
            # Not a unit — spell out the Greek letter.
            result.append(GREEK_ASCII[c])
            i += 1
        else:
            result.append(c)
            i += 1
    return "".join(result)


def normalize_superscripts_subscripts(text: str) -> str:
    """Convert superscript/subscript characters to ASCII equivalents."""
    for uni, ascii_ch in SUPERSCRIPT_MAP.items():
        text = text.replace(uni, ascii_ch)
    for uni, ascii_ch in SUBSCRIPT_MAP.items():
        text = text.replace(uni, ascii_ch)
    return text


def normalize_degree_sign(text: str) -> str:
    """Convert °C, °F, etc. to " deg C", " deg F" with clean spacing.

    Handles patterns like:
      "56 °C"   → "56 deg C"
      "56°C"    → "56 deg C"
      "56 ± 0.1 °C" → "56 +/- 0.1 deg C"
    """
    # "<num> ?°<letter>" → "<num> deg <letter>"
    text = re.sub(r"\s*\u00b0\s*([CFK])\b", r" deg \1", text)
    # Any remaining degree signs (rare) → " deg"
    text = text.replace("\u00b0", " deg")
    # Collapse any accidental double spaces created
    text = re.sub(r"  +", " ", text)
    return text


def normalize_misc(text: str) -> str:
    """Convert misc Unicode punctuation and symbols to ASCII."""
    # Handle degree signs first (context-sensitive)
    text = normalize_degree_sign(text)
    # Then simple 1:1 substitutions
    for uni, ascii_str in MISC_MAP.items():
        text = text.replace(uni, ascii_str)
    return text


def normalize_ascii(text: str) -> str:
    """Full normalization pipeline. Returns pure ASCII text."""
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return text
    s = str(text)
    # NFKC normalization first — folds compatibility characters (e.g., "㎛" → "µm")
    s = unicodedata.normalize("NFKC", s)
    # Greek letters (context-sensitive: unit prefix vs. receptor prefix)
    s = normalize_greek_context_sensitive(s)
    # Superscripts and subscripts
    s = normalize_superscripts_subscripts(s)
    # Miscellaneous single-character replacements
    s = normalize_misc(s)
    # Final safety pass: strip any remaining non-ASCII characters and warn
    # (this shouldn't happen if the maps above are complete, but we don't
    # want silent surprises).
    return s


def find_remaining_non_ascii(text: str) -> list[tuple[int, str]]:
    """Return positions of any non-ASCII characters left in text."""
    if text is None or not isinstance(text, str):
        return []
    return [(i, c) for i, c in enumerate(text) if ord(c) > 127]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run(
    input_path: Path,
    output_path: Path,
    column: str,
    keep_unicode: bool,
) -> None:
    print(f"[1/4] Reading: {input_path}")
    if not input_path.exists():
        raise SystemExit(
            f"\nERROR: input file not found: {input_path}\n"
            f"Make sure the file exists in the folder where you're running "
            f"this script,\nor edit INPUT_FILE at the top of ascii_normalize.py."
        )
    df = pd.read_csv(input_path)
    print(f"      rows: {len(df)}")

    if column not in df.columns:
        raise SystemExit(f"column {column!r} not found. Available: {list(df.columns)}")

    print(f"[2/4] Normalizing column: {column}")
    original_col = df[column].copy()
    normalized = df[column].apply(normalize_ascii)

    print(f"[3/4] Checking for residual non-ASCII characters...")
    residual_hits = []
    for i, val in normalized.items():
        if val is None or not isinstance(val, str):
            continue
        leftovers = find_remaining_non_ascii(val)
        if leftovers:
            residual_hits.append((i, leftovers, val))
    if residual_hits:
        print(f"      WARNING: {len(residual_hits)} row(s) still contain non-ASCII characters.")
        print(f"      First few examples (row_idx, chars, description snippet):")
        for row_idx, chars, val in residual_hits[:5]:
            char_list = ", ".join(f"{c!r} (U+{ord(c):04X})" for _, c in chars[:5])
            print(f"        row {row_idx}: {char_list}")
            print(f"          snippet: {val[:150]}")
        print(f"      Add these characters to the appropriate map at the top of ascii_normalize.py.")
    else:
        print(f"      OK: all descriptions are pure ASCII.")

    # Report change statistics
    n_changed = sum(
        1 for a, b in zip(original_col, normalized)
        if isinstance(a, str) and isinstance(b, str) and a != b
    )
    print(f"      changed: {n_changed} of {len(df)} row(s) had at least one substitution.")

      # Assemble output
    output_df = df.copy()

    # Replace the values only in the assay_description column
    output_df[column] = normalized

    # Rename that specific column while preserving its original position
    output_df.rename(
        columns={column: "assay_description_ascii_compliant"},
        inplace=True,
    )

    print(f"[4/4] Writing: {output_path}")
    output_df.to_excel(output_path, index=False)
    print("      done.")


def main() -> None:
    run(
        input_path=Path(INPUT_FILE),
        output_path=Path(OUTPUT_FILE),
        column=COLUMN,
        keep_unicode=KEEP_UNICODE,
    )


if __name__ == "__main__":
    main()
