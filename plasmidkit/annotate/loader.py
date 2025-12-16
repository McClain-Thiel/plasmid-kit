from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Iterator, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class InvalidSequence(ValueError):
    """Raised when an input sequence cannot be parsed."""


def normalise_sequence(sequence: str) -> str:
    cleaned = sequence.upper().replace("\n", "").replace("\r", "")
    allowed = set("ACGTN")
    filtered = [base for base in cleaned if base in allowed]
    if not filtered:
        raise InvalidSequence("Sequence contains no valid DNA bases")
    return "".join(filtered)


def load_record(source: str | Path | SeqRecord, is_sequence: Optional[bool] = None) -> SeqRecord:
    if isinstance(source, SeqRecord):
        return source
    
    # Path objects are always treated as files
    if isinstance(source, Path):
        if not source.exists():
            raise FileNotFoundError(f"File not found: {source}")
        return next(SeqIO.parse(str(source), infer_format(source)))

    if isinstance(source, str):
        # If caller specifies interpretation, honor it
        if is_sequence is True:
            return SeqRecord(Seq(normalise_sequence(source)), id="sequence")
        
        if is_sequence is False:
            path = Path(source)
            if not path.exists():
                raise FileNotFoundError(f"File not found: {source}")
            return next(SeqIO.parse(str(path), infer_format(path)))

        # Heuristic when not specified: 
        # If it looks like a file path that exists, treat as file.
        path = Path(source)
        if path.exists():
            return next(SeqIO.parse(str(path), infer_format(path)))
        
        # If it's not a file, do not guess. Require explicit is_sequence=True.
        # This prevents ambiguity and potential issues with long strings being treated as paths.
        raise ValueError(
            "Input string is not an existing file. "
            "If you are passing a raw DNA sequence, you must specify 'is_sequence=True'."
        )

    raise TypeError(f"Unsupported source type: {type(source)!r}")


def infer_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna"}:
        return "fasta"
    if suffix in {".gb", ".gbk", ".genbank"}:
        return "genbank"
    raise InvalidSequence(f"Unsupported file format: {path}")


def iter_records(source: str | Path | Iterable[str | Path | SeqRecord], is_sequence: Optional[bool] = None) -> Iterator[SeqRecord]:
    if isinstance(source, (str, Path, SeqRecord)):
        yield load_record(source, is_sequence=is_sequence)
        return
    for item in source:
        yield load_record(item)
