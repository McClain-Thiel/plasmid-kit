from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Mapping, Optional

from Bio.SeqRecord import SeqRecord

from .detectors import run_detectors
from .loader import load_record
from .types import Feature


__all__ = ["annotate_record", "load_record", "Feature"]


def annotate_record(
    record: SeqRecord | str | Path,
    db: Mapping[str, object],
    detectors: Iterable[str] | None = None,
    is_sequence: Optional[bool] = None,
) -> List[Feature]:
    normalized_record = record if isinstance(record, SeqRecord) else load_record(record, is_sequence=is_sequence)
    sequence = str(normalized_record.seq)
    return run_detectors(sequence, db, detectors)
