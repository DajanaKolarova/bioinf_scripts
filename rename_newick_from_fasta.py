#!/usr/bin/env python3
"""Rename terminal Newick leaves from organism names stored in FASTA headers.

Expected FASTA header form, for example:
    >WP_012345678.1 ATP-dependent ligase [Bacillus subtilis]

The first whitespace-delimited token is treated as the sequence/accession ID.
The last text in square brackets is treated as the organism name.

python rename_newick_from_fasta.py \
    vstupni_sekvence.fasta \
    puvodni_strom.nwk \
    strom_organismy.nwk
"""

from __future__ import annotations

import argparse
import collections
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


RANGE_SUFFIX_RE = re.compile(r"/\d+-\d+$")
VERSION_SUFFIX_RE = re.compile(r"\.\d+$")


def iter_fasta_headers(path: Path) -> Iterable[str]:
    """Yield FASTA header text without the leading '>'."""
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}.")
                yield header


def parse_fasta_mapping(path: Path) -> List[Tuple[str, str]]:
    """Return ordered (accession, organism) pairs parsed from FASTA headers."""
    records: List[Tuple[str, str]] = []
    seen_accessions = set()

    for header in iter_fasta_headers(path):
        accession = header.split(maxsplit=1)[0]
        bracket_values = re.findall(r"\[([^\[\]]+)\]", header)
        if not bracket_values:
            raise ValueError(
                f"No organism in square brackets for FASTA ID {accession!r}:\n"
                f">{header}"
            )
        organism = bracket_values[-1].strip()
        if not organism:
            raise ValueError(f"Empty organism name for FASTA ID {accession!r}.")
        if accession in seen_accessions:
            raise ValueError(f"Duplicate FASTA ID/accession: {accession!r}.")
        seen_accessions.add(accession)
        records.append((accession, organism))

    if not records:
        raise ValueError("The input file contains no FASTA headers.")
    return records


def safe_label(text: str) -> str:
    """Create a broadly compatible unquoted Newick/Jalview label."""
    text = text.strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[\(\)\[\]\{\},:;=\"']+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "Unknown_organism"


def build_label_mapping(
    records: List[Tuple[str, str]],
    duplicate_policy: str,
) -> Dict[str, str]:
    """Build accession -> display label mapping."""
    normalized_organisms = [safe_label(org) for _, org in records]
    counts = collections.Counter(normalized_organisms)
    occurrence = collections.Counter()
    mapping: Dict[str, str] = {}

    for (accession, _organism), organism_label in zip(records, normalized_organisms):
        occurrence[organism_label] += 1
        if counts[organism_label] == 1:
            new_label = organism_label
        elif duplicate_policy == "error":
            raise ValueError(
                f"Organism {organism_label!r} occurs {counts[organism_label]} times. "
                "Use --duplicates accession or --duplicates number."
            )
        elif duplicate_policy == "accession":
            new_label = f"{organism_label}__{accession}"
        elif duplicate_policy == "number":
            new_label = f"{organism_label}__{occurrence[organism_label]}"
        else:  # protected by argparse
            raise AssertionError(duplicate_policy)

        mapping[accession] = new_label

    return mapping


def aliases_for_accession(accession: str) -> Iterable[str]:
    """Yield conservative aliases that a tree program may have used."""
    yield accession

    no_range = RANGE_SUFFIX_RE.sub("", accession)
    if no_range != accession:
        yield no_range

    # Common UniProt-style FASTA IDs: sp|P12345|NAME or tr|A0A...|NAME
    parts = accession.split("|")
    if len(parts) >= 2 and parts[1]:
        yield parts[1]

    # Only used as an alias when unique; handled below.
    no_version = VERSION_SUFFIX_RE.sub("", accession)
    if no_version != accession:
        yield no_version


def build_alias_mapping(mapping: Dict[str, str]) -> Dict[str, str]:
    """Add only unambiguous accession aliases."""
    candidates: Dict[str, set[str]] = collections.defaultdict(set)
    for accession, label in mapping.items():
        for alias in aliases_for_accession(accession):
            candidates[alias].add(label)

    return {
        alias: next(iter(labels))
        for alias, labels in candidates.items()
        if len(labels) == 1
    }


def read_quoted_label(text: str, start: int) -> Tuple[str, int]:
    """Read a single-quoted Newick label; doubled quotes escape a quote."""
    assert text[start] == "'"
    i = start + 1
    chars: List[str] = []
    while i < len(text):
        if text[i] == "'":
            if i + 1 < len(text) and text[i + 1] == "'":
                chars.append("'")
                i += 2
                continue
            return "".join(chars), i + 1
        chars.append(text[i])
        i += 1
    raise ValueError("Unterminated single-quoted label in Newick file.")


def quote_newick_label(label: str) -> str:
    """Quote a label only if required by Newick syntax."""
    if label and not re.search(r"[\s\(\)\[\],:;']", label):
        return label
    return "'" + label.replace("'", "''") + "'"


def split_range_suffix(label: str) -> Tuple[str, str]:
    match = RANGE_SUFFIX_RE.search(label)
    if match:
        return label[: match.start()], match.group(0)
    return label, ""


def relabel_newick_leaves(
    text: str,
    alias_mapping: Dict[str, str],
) -> Tuple[str, List[str], int]:
    """Rename terminal leaf labels without touching internal node labels."""
    out: List[str] = []
    unmatched: List[str] = []
    replaced = 0
    i = 0
    expecting_subtree = True  # start of tree, or immediately after '(' / ','

    while i < len(text):
        char = text[i]

        if char.isspace():
            out.append(char)
            i += 1
            continue

        if char == "[":  # Newick comment; preserve, including nested brackets.
            depth = 1
            j = i + 1
            while j < len(text) and depth:
                if text[j] == "[":
                    depth += 1
                elif text[j] == "]":
                    depth -= 1
                j += 1
            if depth:
                raise ValueError("Unterminated comment in Newick file.")
            out.append(text[i:j])
            i = j
            continue

        if char == "(":
            out.append(char)
            expecting_subtree = True
            i += 1
            continue

        if char == ",":
            out.append(char)
            expecting_subtree = True
            i += 1
            continue

        if char == ")":
            out.append(char)
            expecting_subtree = False
            i += 1
            continue

        if char in ":;":
            out.append(char)
            expecting_subtree = False
            i += 1
            continue

        # Read either a quoted or an unquoted token.
        if char == "'":
            token, end = read_quoted_label(text, i)
            original_rendering = text[i:end]
        else:
            end = i
            while end < len(text):
                current = text[end]
                if current.isspace() or current in "(),:;[]":
                    break
                end += 1
            token = text[i:end]
            original_rendering = token

        if not token:
            out.append(original_rendering)
            i = end
            continue

        if expecting_subtree:
            base_token, range_suffix = split_range_suffix(token)
            replacement = alias_mapping.get(token) or alias_mapping.get(base_token)
            if replacement is not None:
                out.append(quote_newick_label(replacement + range_suffix))
                replaced += 1
            else:
                out.append(original_rendering)
                unmatched.append(token)
            expecting_subtree = False
        else:
            # Internal node label, support value, branch length, etc.
            out.append(original_rendering)

        i = end

    return "".join(out), unmatched, replaced


def write_mapping(path: Path, records: List[Tuple[str, str]], mapping: Dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("accession\torganism\tnewick_label\n")
        for accession, organism in records:
            handle.write(f"{accession}\t{organism}\t{mapping[accession]}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Replace terminal accession labels in a Newick tree with organism "
            "names parsed from the original multifasta headers."
        )
    )
    parser.add_argument("fasta", type=Path, help="Original multifasta file")
    parser.add_argument("newick", type=Path, help="Input Newick tree")
    parser.add_argument("output", type=Path, help="Output Newick tree")
    parser.add_argument(
        "--duplicates",
        choices=("accession", "number", "error"),
        default="accession",
        help=(
            "How to distinguish repeated organism names: append accession "
            "(default), append a number, or stop with an error"
        ),
    )
    parser.add_argument(
        "--mapping-tsv",
        type=Path,
        help="Optional TSV report with accession, organism and final label",
    )
    parser.add_argument(
        "--allow-unmatched",
        action="store_true",
        help="Keep unmatched Newick leaf labels instead of exiting with an error",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        records = parse_fasta_mapping(args.fasta)
        mapping = build_label_mapping(records, args.duplicates)
        alias_mapping = build_alias_mapping(mapping)
        tree_text = args.newick.read_text(encoding="utf-8-sig")
        renamed_tree, unmatched, replaced = relabel_newick_leaves(tree_text, alias_mapping)

        unique_unmatched = list(dict.fromkeys(unmatched))
        if unique_unmatched and not args.allow_unmatched:
            preview = ", ".join(repr(x) for x in unique_unmatched[:10])
            more = "" if len(unique_unmatched) <= 10 else f" ... and {len(unique_unmatched) - 10} more"
            raise ValueError(
                "Some terminal Newick labels were not found in the FASTA IDs: "
                f"{preview}{more}. Use --allow-unmatched to leave them unchanged."
            )

        args.output.write_text(renamed_tree, encoding="utf-8")
        if args.mapping_tsv:
            write_mapping(args.mapping_tsv, records, mapping)

        print(f"Renamed terminal leaves: {replaced}")
        print(f"Output tree: {args.output}")
        if unique_unmatched:
            print(
                f"Warning: {len(unique_unmatched)} unmatched leaf label(s) were kept unchanged.",
                file=sys.stderr,
            )
        return 0
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
