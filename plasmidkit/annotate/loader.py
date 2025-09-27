from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator

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


def load_record(source: str | Path | SeqRecord) -> SeqRecord:
    if isinstance(source, SeqRecord):
        return source
    if isinstance(source, (str, Path)):
        path = Path(source)
        if path.exists():
            return next(SeqIO.parse(str(path), infer_format(path)))
        return SeqRecord(Seq(normalise_sequence(str(source))), id="sequence")
    raise TypeError(f"Unsupported source type: {type(source)!r}")


def infer_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna"}:
        return "fasta"
    if suffix in {".gb", ".gbk", ".genbank"}:
        return "genbank"
    raise InvalidSequence(f"Unsupported file format: {path}")


def iter_records(source: str | Path | Iterable[str | Path | SeqRecord]) -> Iterator[SeqRecord]:
    if isinstance(source, (str, Path, SeqRecord)):
        yield load_record(source)
        return
    for item in source:
        yield load_record(item)
