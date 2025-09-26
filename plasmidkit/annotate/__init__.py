from __future__ import annotations

from typing import Iterable, List, Mapping

from Bio.SeqRecord import SeqRecord

from .detectors import run_detectors
from .loader import load_record
from .types import Feature


__all__ = ["annotate_record", "load_record", "Feature"]


def annotate_record(record: SeqRecord, db: Mapping[str, object], detectors: Iterable[str] | None = None) -> List[Feature]:
    sequence = str(record.seq)
    return run_detectors(sequence, db, detectors)
