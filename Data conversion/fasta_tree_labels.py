#!/usr/bin/env python3
"""
Prepare readable, unique FASTA identifiers for MAFFT/Jalview and optionally
relabel an existing aligned FASTA file or Newick tree using the generated map.

Expected input header example:
    >WP_012345678.1 description [Bacillus subtilis]

Default output identifier:
    Bacillus_subtilis__WP_012345678.1

    python fasta_tree_labels.py prepare \
    input_sequences.fasta \
    mafft_input_named.fasta \
    labels.tsv

"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import textwrap
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, TextIO, Tuple


ORGANISM_RE = re.compile(r"\[([^\[\]]+)\]")
JALVIEW_RANGE_RE = re.compile(r"(/\d+-\d+)$")
UNSAFE_LABEL_RE = re.compile(r"[^A-Za-z0-9_.-]+")
MULTIPLE_UNDERSCORES_RE = re.compile(r"_+")
SAFE_NEWICK_LABEL_RE = re.compile(r"^[A-Za-z0-9_.\-/]+$")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


class UserInputError(ValueError):
    """Raised for invalid FASTA, mapping, or Newick input."""


def open_text(path: str, mode: str) -> TextIO:
    """Open a UTF-8 text file, or stdin/stdout when path is '-'."""
    if path == "-":
        if "r" in mode:
            return nullcontext(sys.stdin)  # type: ignore[return-value]
        return nullcontext(sys.stdout)  # type: ignore[return-value]
    return open(Path(path), mode, encoding="utf-8", newline="")


def read_fasta(handle: TextIO) -> Iterator[FastaRecord]:
    header: str | None = None
    sequence_parts: list[str] = []
    line_number = 0

    for raw_line in handle:
        line_number += 1
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                sequence = "".join(sequence_parts)
                if not sequence:
                    raise UserInputError(
                        f"FASTA record '{header}' has no sequence."
                    )
                yield FastaRecord(header=header, sequence=sequence)

            header = line[1:].strip()
            if not header:
                raise UserInputError(
                    f"Empty FASTA header at line {line_number}."
                )
            sequence_parts = []
        else:
            if header is None:
                raise UserInputError(
                    f"Sequence data found before the first FASTA header "
                    f"at line {line_number}."
                )
            # Remove all whitespace from sequence lines while preserving symbols.
            sequence_parts.append("".join(line.split()))

    if header is not None:
        sequence = "".join(sequence_parts)
        if not sequence:
            raise UserInputError(f"FASTA record '{header}' has no sequence.")
        yield FastaRecord(header=header, sequence=sequence)


def write_fasta_record(
    handle: TextIO,
    identifier: str,
    sequence: str,
    description: str = "",
    width: int = 80,
) -> None:
    header = identifier if not description else f"{identifier} {description}"
    handle.write(f">{header}\n")
    if width <= 0:
        handle.write(f"{sequence}\n")
    else:
        for chunk in textwrap.wrap(sequence, width=width):
            handle.write(f"{chunk}\n")


def first_token(header: str) -> str:
    token = header.split(maxsplit=1)[0]
    if not token:
        raise UserInputError(f"Cannot obtain an identifier from header: {header}")
    return token


def description_after_first_token(header: str) -> str:
    parts = header.split(maxsplit=1)
    return parts[1] if len(parts) == 2 else ""


def extract_organism(header: str) -> str | None:
    """Return the last bracketed value in a FASTA header."""
    matches = ORGANISM_RE.findall(header)
    if not matches:
        return None
    organism = matches[-1].strip()
    return organism or None


def sanitize_label(value: str) -> str:
    """Convert a value into a conservative FASTA/Newick-safe token."""
    cleaned = UNSAFE_LABEL_RE.sub("_", value.strip())
    cleaned = MULTIPLE_UNDERSCORES_RE.sub("_", cleaned).strip("_.-")
    return cleaned or "Unknown"


def build_base_label(accession: str, organism: str, style: str) -> str:
    accession_label = sanitize_label(accession)
    organism_label = sanitize_label(organism)

    if style == "organism-accession":
        return f"{organism_label}__{accession_label}"
    if style == "accession-organism":
        return f"{accession_label}__{organism_label}"
    if style == "organism":
        return organism_label
    raise AssertionError(f"Unknown label style: {style}")


def unique_label(base: str, used: set[str]) -> str:
    if base not in used:
        used.add(base)
        return base

    index = 2
    while f"{base}__{index}" in used:
        index += 1
    label = f"{base}__{index}"
    used.add(label)
    return label


def prepare_fasta(args: argparse.Namespace) -> int:
    mapping_rows: list[dict[str, str]] = []
    used_new_ids: set[str] = set()
    seen_old_ids: set[str] = set()

    with open_text(args.input, "r") as input_handle, open_text(
        args.output, "w"
    ) as output_handle:
        count = 0
        for record in read_fasta(input_handle):
            accession = first_token(record.header)
            if accession in seen_old_ids:
                raise UserInputError(
                    f"Duplicate FASTA identifier '{accession}'. Tree leaves cannot "
                    "be mapped unambiguously; make the original identifiers unique."
                )
            seen_old_ids.add(accession)

            organism = extract_organism(record.header)
            if organism is None:
                if args.on_missing_organism == "error":
                    raise UserInputError(
                        f"No bracketed organism name found in header: {record.header}"
                    )
                if args.on_missing_organism == "accession":
                    organism = accession
                else:
                    organism = "Unknown_organism"

            base_label = build_base_label(accession, organism, args.style)
            new_id = unique_label(base_label, used_new_ids)

            description = record.header if args.keep_description else ""
            write_fasta_record(
                output_handle,
                identifier=new_id,
                sequence=record.sequence,
                description=description,
                width=args.width,
            )

            mapping_rows.append(
                {
                    "old_id": accession,
                    "new_id": new_id,
                    "organism": organism,
                    "original_header": record.header,
                }
            )
            count += 1

    with open_text(args.mapping, "w") as mapping_handle:
        writer = csv.DictWriter(
            mapping_handle,
            fieldnames=["old_id", "new_id", "organism", "original_header"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    print(
        f"Prepared {count} sequences. FASTA: {args.output}; mapping: {args.mapping}",
        file=sys.stderr,
    )
    return 0


def load_mapping(path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    new_ids: set[str] = set()

    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise UserInputError(f"Mapping file '{path}' is empty.")

        required = {"old_id", "new_id"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise UserInputError(
                f"Mapping file '{path}' is missing columns: {', '.join(sorted(missing))}"
            )

        for row_number, row in enumerate(reader, start=2):
            old_id = (row.get("old_id") or "").strip()
            new_id = (row.get("new_id") or "").strip()
            if not old_id or not new_id:
                raise UserInputError(
                    f"Empty old_id or new_id in mapping row {row_number}."
                )
            if old_id in mapping:
                raise UserInputError(
                    f"Duplicate old_id '{old_id}' in mapping file."
                )
            if new_id in new_ids:
                raise UserInputError(
                    f"Duplicate new_id '{new_id}' in mapping file."
                )
            mapping[old_id] = new_id
            new_ids.add(new_id)

    if not mapping:
        raise UserInputError(f"Mapping file '{path}' contains no records.")
    return mapping


def map_identifier(identifier: str, mapping: Dict[str, str]) -> str | None:
    """Map an exact ID and preserve an optional Jalview /start-end suffix."""
    if identifier in mapping:
        return mapping[identifier]

    match = JALVIEW_RANGE_RE.search(identifier)
    if match:
        suffix = match.group(1)
        base = identifier[: match.start()]
        if base in mapping:
            return f"{mapping[base]}{suffix}"

    return None


def apply_mapping_to_fasta(args: argparse.Namespace) -> int:
    mapping = load_mapping(args.mapping)
    unmatched: set[str] = set()
    count = 0

    with open_text(args.input, "r") as input_handle, open_text(
        args.output, "w"
    ) as output_handle:
        for record in read_fasta(input_handle):
            old_id = first_token(record.header)
            new_id = map_identifier(old_id, mapping)
            if new_id is None:
                unmatched.add(old_id)
                new_id = old_id

            description = (
                description_after_first_token(record.header)
                if args.keep_description
                else ""
            )
            write_fasta_record(
                output_handle,
                identifier=new_id,
                sequence=record.sequence,
                description=description,
                width=args.width,
            )
            count += 1

    if unmatched and not args.allow_unmapped:
        preview = ", ".join(sorted(unmatched)[:10])
        extra = "" if len(unmatched) <= 10 else f" ... (+{len(unmatched) - 10} more)"
        raise UserInputError(
            f"No mapping found for {len(unmatched)} FASTA identifier(s): "
            f"{preview}{extra}. Output was written, but treat it as incomplete. "
            "Use --allow-unmapped only when this is intentional."
        )

    print(
        f"Relabelled {count - len(unmatched)} of {count} FASTA records. "
        f"Output: {args.output}",
        file=sys.stderr,
    )
    return 0


def read_comment(text: str, start: int) -> Tuple[str, int]:
    """Read a possibly nested Newick/NEXUS comment beginning with '['."""
    depth = 0
    index = start
    while index < len(text):
        char = text[index]
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                return text[start : index + 1], index + 1
        index += 1
    raise UserInputError("Unterminated Newick comment.")


def read_quoted_newick_label(text: str, start: int) -> Tuple[str, int]:
    """Decode a single-quoted Newick label, including doubled apostrophes."""
    assert text[start] == "'"
    chars: list[str] = []
    index = start + 1

    while index < len(text):
        char = text[index]
        if char == "'":
            if index + 1 < len(text) and text[index + 1] == "'":
                chars.append("'")
                index += 2
                continue
            return "".join(chars), index + 1
        chars.append(char)
        index += 1

    raise UserInputError("Unterminated quoted Newick label.")


def quote_newick_label(label: str) -> str:
    if SAFE_NEWICK_LABEL_RE.fullmatch(label):
        return label
    return "'" + label.replace("'", "''") + "'"


def relabel_newick_text(
    text: str,
    mapping: Dict[str, str],
) -> Tuple[str, int, set[str]]:
    """
    Relabel only leaf labels, not internal node/bootstrap labels.

    A leaf is expected at the beginning of a tree, after '(' or after ','.
    Comments and whitespace are preserved.
    """
    output: list[str] = []
    unmatched: set[str] = set()
    replaced = 0
    index = 0
    expect_leaf = True

    while index < len(text):
        char = text[index]

        if char == "[":
            comment, index = read_comment(text, index)
            output.append(comment)
            continue

        if char == "(":
            output.append(char)
            index += 1
            expect_leaf = True
            continue

        if char == ",":
            output.append(char)
            index += 1
            expect_leaf = True
            continue

        if char == ")":
            output.append(char)
            index += 1
            expect_leaf = False
            continue

        if char == ";":
            output.append(char)
            index += 1
            expect_leaf = True
            continue

        if expect_leaf and char.isspace():
            output.append(char)
            index += 1
            continue

        if expect_leaf and char == "'":
            label, next_index = read_quoted_newick_label(text, index)
            mapped = map_identifier(label, mapping)
            if mapped is None:
                unmatched.add(label)
                output.append(text[index:next_index])
            else:
                output.append(quote_newick_label(mapped))
                replaced += 1
            index = next_index
            expect_leaf = False
            continue

        if expect_leaf and char not in ":,();":
            next_index = index
            while next_index < len(text):
                next_char = text[next_index]
                if next_char.isspace() or next_char in "[]():,;":
                    break
                next_index += 1

            label = text[index:next_index]
            if not label:
                output.append(char)
                index += 1
                continue

            mapped = map_identifier(label, mapping)
            if mapped is None:
                unmatched.add(label)
                output.append(label)
            else:
                output.append(quote_newick_label(mapped))
                replaced += 1
            index = next_index
            expect_leaf = False
            continue

        output.append(char)
        index += 1

    return "".join(output), replaced, unmatched


def apply_mapping_to_tree(args: argparse.Namespace) -> int:
    mapping = load_mapping(args.mapping)

    with open_text(args.input, "r") as input_handle:
        tree_text = input_handle.read()
    if not tree_text.strip():
        raise UserInputError(f"Newick file '{args.input}' is empty.")

    output_text, replaced, unmatched = relabel_newick_text(tree_text, mapping)

    with open_text(args.output, "w") as output_handle:
        output_handle.write(output_text)
        if output_text and not output_text.endswith("\n"):
            output_handle.write("\n")

    if unmatched and not args.allow_unmapped:
        preview = ", ".join(sorted(unmatched)[:10])
        extra = "" if len(unmatched) <= 10 else f" ... (+{len(unmatched) - 10} more)"
        raise UserInputError(
            f"No mapping found for {len(unmatched)} Newick leaf label(s): "
            f"{preview}{extra}. Output was written, but treat it as incomplete. "
            "Use --allow-unmapped only when this is intentional."
        )

    print(
        f"Relabelled {replaced} Newick leaves. Output: {args.output}",
        file=sys.stderr,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create readable organism-based FASTA IDs for MAFFT/Jalview and "
            "relabel aligned FASTA or Newick files consistently."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        help="Rename an original FASTA and create a TSV mapping file.",
    )
    prepare.add_argument("input", help="Input multi-FASTA file, or '-' for stdin.")
    prepare.add_argument("output", help="Output renamed FASTA file, or '-' for stdout.")
    prepare.add_argument("mapping", help="Output TSV mapping file.")
    prepare.add_argument(
        "--style",
        choices=["organism-accession", "accession-organism", "organism"],
        default="organism-accession",
        help=(
            "Identifier style. Default: organism-accession, e.g. "
            "Bacillus_subtilis__WP_012345.1."
        ),
    )
    prepare.add_argument(
        "--on-missing-organism",
        choices=["error", "unknown", "accession"],
        default="error",
        help="What to do when no [organism] is present. Default: error.",
    )
    prepare.add_argument(
        "--keep-description",
        action="store_true",
        help="Keep the complete original header after the new first-token ID.",
    )
    prepare.add_argument(
        "--width",
        type=int,
        default=80,
        help="Sequence line width; use 0 for one line. Default: 80.",
    )
    prepare.set_defaults(func=prepare_fasta)

    fasta = subparsers.add_parser(
        "fasta",
        help="Apply an existing mapping to a FASTA/aligned FASTA file.",
    )
    fasta.add_argument("input", help="Input FASTA file.")
    fasta.add_argument("mapping", help="TSV mapping made by the prepare command.")
    fasta.add_argument("output", help="Output relabelled FASTA file.")
    fasta.add_argument(
        "--keep-description",
        action="store_true",
        help="Preserve text after the first identifier token.",
    )
    fasta.add_argument(
        "--allow-unmapped",
        action="store_true",
        help="Leave unknown identifiers unchanged instead of returning an error.",
    )
    fasta.add_argument(
        "--width",
        type=int,
        default=80,
        help="Sequence line width; use 0 for one line. Default: 80.",
    )
    fasta.set_defaults(func=apply_mapping_to_fasta)

    tree = subparsers.add_parser(
        "tree",
        help="Apply an existing mapping to leaf labels in a Newick tree.",
    )
    tree.add_argument("input", help="Input Newick file.")
    tree.add_argument("mapping", help="TSV mapping made by the prepare command.")
    tree.add_argument("output", help="Output relabelled Newick file.")
    tree.add_argument(
        "--allow-unmapped",
        action="store_true",
        help="Leave unknown leaf labels unchanged instead of returning an error.",
    )
    tree.set_defaults(func=apply_mapping_to_tree)

    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if hasattr(args, "width") and args.width < 0:
        parser.error("--width must be 0 or a positive integer.")

    try:
        return int(args.func(args))
    except (OSError, UserInputError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
