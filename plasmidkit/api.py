from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Mapping, Sequence, Optional

from Bio.SeqRecord import SeqRecord

from .annotate import annotate_record, load_record
from .annotate.types import Feature
from .cache import manager
from .exporters import export_gff3, export_json, export_minimal_genbank
from .scoring.calculator import compute_score

__all__ = [
    "load_record",
    "annotate",
    "score",
    "annotate_and_score",
    "export_json",
    "export_gff3",
    "export_minimal_genbank",
    "set_cache_dir",
    "set_offline",
    "add_registry",
]


def annotate(
    record: SeqRecord | str | Path,
    db: str = "engineered-core@1.0.0",
    detectors: Iterable[str] | None = None,
    is_sequence: Optional[bool] = None,
) -> List[Feature]:
    artifacts = manager.get_artifacts(db)
    return annotate_record(record, artifacts, detectors, is_sequence=is_sequence)


def score(
    record: SeqRecord | str | Path,
    annotations: Sequence[Feature] | None = None,
    db: str = "engineered-core@1.0.0",
    is_sequence: Optional[bool] = None,
) -> Mapping[str, object]:
    artifacts = manager.get_artifacts(db)
    normalized_record = record if isinstance(record, SeqRecord) else load_record(record, is_sequence=is_sequence)
    features = list(annotations) if annotations is not None else annotate_record(normalized_record, artifacts, None)
    return compute_score(normalized_record, features, artifacts)


def annotate_and_score(
    record: SeqRecord | str | Path,
    db: str = "engineered-core@1.0.0",
    detectors: Iterable[str] | None = None,
    is_sequence: Optional[bool] = None,
) -> Mapping[str, object]:
    normalized_record = record if isinstance(record, SeqRecord) else load_record(record, is_sequence=is_sequence)
    annotations = annotate(normalized_record, db=db, detectors=detectors)
    score_report = score(normalized_record, annotations=annotations, db=db)
    return {
        "sequence_id": normalized_record.id,
        "length": len(normalized_record.seq),
        "annotations": [feature.to_dict() for feature in annotations],
        "score": score_report,
        "db": db,
    }


set_cache_dir = manager.set_cache_dir
set_offline = manager.set_offline
add_registry = manager.add_registry
