from __future__ import annotations

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
        if source.exists():
            return next(SeqIO.parse(str(source), infer_format(source)))
        # If an explicit Path does not exist, fall back to sequence interpretation
        return SeqRecord(Seq(normalise_sequence(str(source))), id="sequence")
    if isinstance(source, str):
        # If caller specifies interpretation, honor it
        if is_sequence is True:
            return SeqRecord(Seq(normalise_sequence(source)), id="sequence")
        if is_sequence is False:
            path = Path(source)
            if path.exists():
                return next(SeqIO.parse(str(path), infer_format(path)))
            # Backward-compatible fallback: interpret as sequence if file not found
            return SeqRecord(Seq(normalise_sequence(source)), id="sequence")
        # Heuristic when not specified: long strings are likely raw sequences
        if len(source) >= 1000:
            return SeqRecord(Seq(normalise_sequence(source)), id="sequence")
        # Otherwise prefer file if present, else interpret as sequence
        path = Path(source)
        if path.exists():
            return next(SeqIO.parse(str(path), infer_format(path)))
        return SeqRecord(Seq(normalise_sequence(source)), id="sequence")
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
