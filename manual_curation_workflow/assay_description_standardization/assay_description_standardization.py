"""
run_curation.py
================

Generate assay descriptions for a dataset of experimental protocols
using Anthropic's Claude API.

Design (matches the locked rules in curation_rules.md):
  - Dedup key:  (experimental_protocol, protein_name)
  - Output:     one description per input row, in original row order
  - Guarantees: identical (protocol, protein) pairs receive byte-identical
                descriptions; row count and order preserved

Quick start (with default demo file names):

    Ensure Assay_standardization_input_dataset.csv and system_prompt.txt
    are in the current directory, then run:

        python run_curation.py

    This reads Assay_standardization_input_dataset.csv and writes
    Assay_standardization_output_dataset.xlsx in production mode.

Custom file names:

    python run_curation.py \\
        --input  my_data.csv \\
        --output my_output.xlsx \\
        --system-prompt my_prompt.txt \\
        --mode production

Validation pass (returns JSON metadata for spot-checking each unique
(protocol, protein) tuple):

    python run_curation.py --mode validation

Required environment variable:
    ANTHROPIC_API_KEY

Required input columns (case-sensitive, configurable via --col-* flags):
    experimental_protocol, assay_name, protein_name

Output columns:
    All input columns + assay_description
    (validation mode also adds: assay_class, template_used, platform,
     stimulating_agent, labeled_ligand, readout, incubation_retained)

Dependencies:
    pip install anthropic pandas openpyxl tenacity tqdm
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
from anthropic import Anthropic, APIError, APIStatusError
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "claude-opus-4-7"
DEFAULT_MAX_TOKENS = 1024
DEFAULT_TEMPERATURE = None  # Opus 4.7 deprecates this; pass None to omit it.

# Banned tokens that must not appear in any generated description.
# These belong in the standard_endpoint column, not the description.
BANNED_TOKENS = [
    r"\bIC50\b",
    r"\bEC50\b",
    r"\bKi\b",
    r"\bKd\b",
    r"\bEmax\b",
    r"\bpIC50\b",
    r"\bpEC50\b",
    r"\bpKi\b",
    r"\bMPE\b",
    r"\bAD50\b",
    r"\bED50\b",
    r"according to",
    r"published procedure",
]
BANNED_PATTERN = re.compile("|".join(BANNED_TOKENS), re.IGNORECASE)

# Heuristic markers used to detect in vivo (Template C) descriptions for
# the concentration-response / dose-response consistency check.
IN_VIVO_MARKERS = re.compile(
    r"\b(mice|mouse|rat|rats|swiss webster|c57|sprague|wistar|"
    r"subcutaneous|s\.c\.|intraperitoneal|i\.p\.|intravenous|i\.v\.|"
    r"oral(ly)?|p\.o\.|intrathecal|i\.c\.v\.|"
    r"tail[- ]withdrawal|tail[- ]flick|tail[- ]immersion|hot[- ]plate|"
    r"formalin test|antinociception|antinociceptive|locomotor)\b",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DedupKey:
    """The locked dedup key: (experimental_protocol, protein_name)."""

    experimental_protocol: str
    protein_name: str

    @classmethod
    def from_row(
        cls, row: pd.Series, protocol_col: str, protein_col: str
    ) -> "DedupKey":
        proto = "" if pd.isna(row[protocol_col]) else str(row[protocol_col])
        prot = "" if pd.isna(row[protein_col]) else str(row[protein_col])
        # Light normalization: strip + collapse internal whitespace.
        # The original values stay in the dataframe; this only affects the key.
        proto_norm = re.sub(r"\s+", " ", proto).strip()
        prot_norm = re.sub(r"\s+", " ", prot).strip()
        return cls(experimental_protocol=proto_norm, protein_name=prot_norm)


# ---------------------------------------------------------------------------
# Anthropic client wrapper
# ---------------------------------------------------------------------------


class CurationClient:
    def __init__(
        self,
        system_prompt: str,
        model: str = DEFAULT_MODEL,
        max_tokens: int = DEFAULT_MAX_TOKENS,
        temperature: float = DEFAULT_TEMPERATURE,
    ):
        self.client = Anthropic()  # picks up ANTHROPIC_API_KEY from env
        self.system_prompt = system_prompt
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

    def _user_message(
        self,
        mode: str,
        protocol: str,
        protein_name: str,
        assay_name: str,
    ) -> str:
        mode_str = "JSON" if mode == "validation" else "TEXT"
        return (
            f"Output mode: {mode_str}\n\n"
            f"experimental_protocol:\n{protocol}\n\n"
            f"protein_name: {protein_name}\n"
            f"assay_name: {assay_name}"
        )

    @retry(
        retry=retry_if_exception_type((APIError, APIStatusError)),
        wait=wait_exponential(multiplier=2, min=2, max=60),
        stop=stop_after_attempt(5),
        reraise=True,
    )
    def call(
        self,
        mode: str,
        protocol: str,
        protein_name: str,
        assay_name: str,
    ) -> str:
        # Build kwargs conditionally — newer models (e.g., Opus 4.7) deprecate
        # the `temperature` parameter and reject requests that include it.
        kwargs = {
            "model": self.model,
            "max_tokens": self.max_tokens,
            "system": self.system_prompt,
            "messages": [
                {
                    "role": "user",
                    "content": self._user_message(
                        mode, protocol, protein_name, assay_name
                    ),
                }
            ],
        }
        if self.temperature is not None:
            kwargs["temperature"] = self.temperature
        msg = self.client.messages.create(**kwargs)
        # Concatenate all text blocks.
        text = "".join(
            block.text for block in msg.content if getattr(block, "type", "") == "text"
        )
        return text.strip()


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_description(desc: str) -> Optional[str]:
    """Return an error message if the description violates rules, else None."""
    if not desc or not desc.strip():
        return "empty description"
    if BANNED_PATTERN.search(desc):
        match = BANNED_PATTERN.search(desc)
        return f"banned token present: {match.group(0)!r}"
    if "```" in desc or desc.strip().startswith("{"):
        return "looks like JSON or Markdown fence in production output"

    # Concentration-response vs. dose-response consistency check.
    is_in_vivo = bool(IN_VIVO_MARKERS.search(desc))
    has_dose = "dose-response" in desc.lower() or "dose\u2013response" in desc.lower()
    has_conc = (
        "concentration-response" in desc.lower()
        or "concentration\u2013response" in desc.lower()
    )

    if is_in_vivo and has_conc and not has_dose:
        return ("in vivo description uses 'concentration-response' "
                "instead of 'dose-response'")
    if not is_in_vivo and has_dose and not has_conc:
        return ("in vitro description uses 'dose-response' "
                "instead of 'concentration-response'")

    return None


def assert_invariants(input_df: pd.DataFrame, output_df: pd.DataFrame) -> None:
    """Hard checks before writing output."""
    assert len(output_df) == len(input_df), (
        f"row count mismatch: input={len(input_df)} output={len(output_df)}"
    )
    assert (output_df.index == input_df.index).all(), "row order changed"
    assert output_df["assay_description"].notna().all(), "null descriptions present"
    assert (output_df["assay_description"].astype(str).str.len() > 0).all(), (
        "empty descriptions present"
    )


def assert_duplicate_consistency(
    df: pd.DataFrame, protocol_col: str, protein_col: str
) -> None:
    """All rows sharing (protocol, protein) must have identical descriptions."""
    grouped = df.groupby([protocol_col, protein_col])["assay_description"].nunique()
    bad = grouped[grouped > 1]
    if not bad.empty:
        sample = bad.head(3).to_dict()
        raise AssertionError(
            f"duplicate inconsistency: {len(bad)} (protocol, protein) groups "
            f"have multiple distinct descriptions. Sample: {sample}"
        )


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def parse_validation_json(raw: str) -> dict:
    """Strip code fences/preamble and parse JSON from validation-mode output.

    Tolerates: markdown code fences, preamble before the JSON, trailing text
    after the JSON, single-quoted JSON, and missing outer braces if the model
    returns just the inner schema.
    """
    cleaned = raw.strip()
    # Remove markdown code fences if present
    cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
    cleaned = re.sub(r"\s*```$", "", cleaned)
    cleaned = cleaned.strip()

    # First try a clean parse
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        pass

    # Fallback: find the first balanced JSON object in the text and parse that.
    start = cleaned.find("{")
    if start == -1:
        raise json.JSONDecodeError("no JSON object found in output", cleaned, 0)
    depth = 0
    end = -1
    in_string = False
    escape = False
    for i in range(start, len(cleaned)):
        c = cleaned[i]
        if in_string:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                in_string = False
            continue
        if c == '"':
            in_string = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                end = i
                break
    if end == -1:
        raise json.JSONDecodeError(
            "unbalanced JSON braces", cleaned, start
        )
    return json.loads(cleaned[start:end + 1])


def run(
    input_path: Path,
    output_path: Path,
    system_prompt_path: Path,
    mode: str,
    protocol_col: str,
    assay_name_col: str,
    protein_col: str,
    model: str,
    max_tokens: int,
    temperature: float,
    sleep_between_calls: float = 0.0,
) -> None:
    print(f"[1/6] Reading input: {input_path}")
    input_df = pd.read_csv(input_path)

    for col in (protocol_col, assay_name_col, protein_col):
        if col not in input_df.columns:
            raise SystemExit(f"missing required column: {col!r}")

    print(f"      input rows: {len(input_df)}")

    print("[2/6] Loading system prompt")
    system_prompt = system_prompt_path.read_text(encoding="utf-8")

    print("[3/6] Building dedup key cache")
    keys = [
        DedupKey.from_row(row, protocol_col, protein_col)
        for _, row in input_df.iterrows()
    ]
    unique_keys = list(dict.fromkeys(keys))  # preserves first-seen order
    print(f"      unique (protocol, protein) tuples: {len(unique_keys)}")
    print(f"      duplication factor: {len(keys) / max(len(unique_keys), 1):.1f}x")

    # Map representative assay_name per unique key (first occurrence).
    rep_assay_name: dict[DedupKey, str] = {}
    for k, (_, row) in zip(keys, input_df.iterrows()):
        if k not in rep_assay_name:
            an = row[assay_name_col]
            rep_assay_name[k] = "" if pd.isna(an) else str(an)

    print(f"[4/6] Calling Anthropic API ({mode} mode, model={model})")
    client = CurationClient(
        system_prompt=system_prompt,
        model=model,
        max_tokens=max_tokens,
        temperature=temperature,
    )

    cache: dict[DedupKey, dict] = {}
    errors: list[tuple[DedupKey, str]] = []
    first_error_printed = False

    for key in tqdm(unique_keys, desc="curating"):
        try:
            raw = client.call(
                mode=mode,
                protocol=key.experimental_protocol,
                protein_name=key.protein_name,
                assay_name=rep_assay_name[key],
            )
            if mode == "validation":
                try:
                    obj = parse_validation_json(raw)
                except Exception as parse_err:
                    if not first_error_printed:
                        print(f"\n\n[DEBUG] First JSON parse failure on "
                              f"protein={key.protein_name!r}")
                        print(f"[DEBUG] Parse error: {parse_err!r}")
                        print(f"[DEBUG] Raw output (first 500 chars):")
                        print(f"---\n{raw[:500]}\n---\n")
                        first_error_printed = True
                    raise
                desc = obj.get("assay_description", "").strip()
                cache[key] = {"_raw": raw, "_meta": obj, "assay_description": desc}
            else:
                cache[key] = {"_raw": raw, "_meta": None, "assay_description": raw}

            err = validate_description(cache[key]["assay_description"])
            if err:
                errors.append((key, err))

        except Exception as e:
            if not first_error_printed:
                import traceback
                print(f"\n\n[DEBUG] First API failure on "
                      f"protein={key.protein_name!r}")
                print(f"[DEBUG] Exception type: {type(e).__name__}")
                print(f"[DEBUG] Exception repr: {e!r}")
                print(f"[DEBUG] Traceback:")
                traceback.print_exc()
                print()
                first_error_printed = True
            errors.append((key, f"API call failed: {type(e).__name__}: {e}"))
            cache[key] = {"_raw": "", "_meta": None, "assay_description": ""}

        if sleep_between_calls > 0:
            time.sleep(sleep_between_calls)

    # AUTO-FIX: response-phrase mismatch.
    # If a description was generated with the wrong response-phrase
    # (e.g., "concentration-response" in an in vivo description), apply
    # the correct phrase mechanically rather than aborting the whole run.
    # This is safe because the substitution is unambiguous in either
    # direction and only fires when the validator detects a mismatch.
    auto_fixed_keys: list[DedupKey] = []
    for key in unique_keys:
        desc = cache[key]["assay_description"]
        if not isinstance(desc, str) or not desc:
            continue
        is_in_vivo = bool(IN_VIVO_MARKERS.search(desc))
        if is_in_vivo:
            new_desc = re.sub(
                r"concentration([-\u2013])response",
                r"dose\1response",
                desc,
                flags=re.IGNORECASE,
            )
        else:
            new_desc = re.sub(
                r"dose([-\u2013])response",
                r"concentration\1response",
                desc,
                flags=re.IGNORECASE,
            )
        if new_desc != desc:
            cache[key]["assay_description"] = new_desc
            if cache[key].get("_meta"):
                cache[key]["_meta"]["assay_description"] = new_desc
            auto_fixed_keys.append(key)
    if auto_fixed_keys:
        print(f"\n[AUTO-FIX] Corrected response-phrase wording in "
              f"{len(auto_fixed_keys)} description(s).")

    # Re-validate after auto-fixes
    errors = []
    for key in unique_keys:
        desc = cache[key]["assay_description"]
        if not desc:
            errors.append((key, "empty description (API call failed)"))
            continue
        err = validate_description(desc)
        if err:
            errors.append((key, err))

    if errors:
        print(f"\n[!!] {len(errors)} validation errors remaining out of "
              f"{len(unique_keys)} unique tuples:")
        for k, err in errors[:10]:
            print(f"     - protein={k.protein_name!r}: {err}")
        if len(errors) > 10:
            print(f"     ... and {len(errors) - 10} more")

    print("[5/6] Mapping descriptions back to all input rows")
    output_df = input_df.copy()
    output_df["assay_description"] = [
        cache[k]["assay_description"] for k in keys
    ]

    if mode == "validation":
        # Add per-tuple metadata columns to support spot-checking.
        meta_cols = [
            "assay_class",
            "template_used",
            "target",
            "target_resolution",
            "biological_system",
            "platform",
            "stimulating_agent",
            "labeled_ligand",
            "readout",
            "incubation_retained",
            "route_of_administration",
            "stimulus_parameters",
            "mode_handling",
            "reference_agonist",
            "response_phrase",
        ]
        for col in meta_cols:
            output_df[col] = [
                (cache[k]["_meta"] or {}).get(col) for k in keys
            ]

    # Always include the raw model output so failed rows can be inspected.
    output_df["_raw_model_output"] = [cache[k]["_raw"] for k in keys]

    # Write a debug file BEFORE running invariants/abort. If anything goes
    # wrong downstream, this file always preserves the API responses so
    # nothing is lost.
    debug_path = output_path.with_name(output_path.stem + "_debug" + output_path.suffix)
    print(f"      writing debug file: {debug_path}")
    output_df.to_excel(debug_path, index=False)

    # If any remaining (post-auto-fix) errors exist in production mode,
    # abort here — but only AFTER the debug file is safely written.
    if errors and mode == "production":
        raise SystemExit(
            f"\n[ABORT] {len(errors)} validation error(s) remain after "
            f"auto-fix. Debug file saved to: {debug_path}\n"
            "Inspect the _raw_model_output column to see what the model "
            "produced. Either fix the prompt and re-run, or manually "
            "correct the offending descriptions in the debug file and "
            "use it as your final output."
        )

    print("[6/6] Asserting invariants")
    try:
        assert_invariants(input_df, output_df)
        assert_duplicate_consistency(output_df, protocol_col, protein_col)
    except AssertionError as e:
        print(f"\n[!!] Invariant check failed: {e}")
        print(f"      Debug file with raw outputs was saved to: {debug_path}")
        print(f"      Inspect the _raw_model_output column for empty/failed rows.")
        # Show which rows have empty descriptions
        empty_mask = output_df["assay_description"].astype(str).str.len() == 0
        if empty_mask.any():
            empty_rows = output_df[empty_mask].head(5)
            print(f"\n      First few empty rows (showing protein_name and "
                  f"raw output preview):")
            for idx, row in empty_rows.iterrows():
                raw_preview = str(row["_raw_model_output"])[:200]
                print(f"        row {idx}: protein={row[protein_col]!r}")
                print(f"          raw[:200]={raw_preview!r}")
        raise SystemExit(1)

    # Drop the debug column for the clean final output
    final_df = output_df.drop(columns=["_raw_model_output"])
    print(f"      writing: {output_path}")
    final_df.to_excel(output_path, index=False)
    print(f"      done. {len(final_df)} rows written, "
          f"{len(unique_keys)} unique tuples curated.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path,
                   default=Path("Assay_standardization_input_dataset.csv"),
                   help="Input CSV file (default: Assay_standardization_input_dataset.csv)")
    p.add_argument("--output", type=Path,
                   default=Path("Assay_standardization_output_dataset.xlsx"),
                   help="Output XLSX file (default: Assay_standardization_output_dataset.xlsx)")
    p.add_argument("--system-prompt", type=Path,
                   default=Path("system_prompt.txt"),
                   help="Plain-text file containing the system prompt body "
                        "(default: system_prompt.txt).")
    p.add_argument("--mode", choices=["production", "validation"],
                   default="production")
    p.add_argument("--col-protocol", default="experimental_protocol")
    p.add_argument("--col-assay-name", default="assay_name")
    p.add_argument("--col-protein", default="protein_name")
    p.add_argument("--model", default=DEFAULT_MODEL)
    p.add_argument("--max-tokens", type=int, default=DEFAULT_MAX_TOKENS)
    p.add_argument("--temperature", type=float, default=DEFAULT_TEMPERATURE,
                   help="Sampling temperature. Default is unset (omitted from "
                        "the API call) because Opus 4.7 deprecates this "
                        "parameter. Older models default to 0 if you pass "
                        "--temperature 0 explicitly.")
    p.add_argument("--sleep", type=float, default=0.0,
                   help="Seconds to sleep between API calls (rate limiting).")
    args = p.parse_args()

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("ERROR: ANTHROPIC_API_KEY environment variable is not set.",
              file=sys.stderr)
        sys.exit(2)

    run(
        input_path=args.input,
        output_path=args.output,
        system_prompt_path=args.system_prompt,
        mode=args.mode,
        protocol_col=args.col_protocol,
        assay_name_col=args.col_assay_name,
        protein_col=args.col_protein,
        model=args.model,
        max_tokens=args.max_tokens,
        temperature=args.temperature,
        sleep_between_calls=args.sleep,
    )


if __name__ == "__main__":
    main()
