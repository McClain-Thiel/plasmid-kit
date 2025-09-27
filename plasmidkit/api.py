from __future__ import annotations

from typing import Iterable, List, Mapping, Sequence

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


def annotate(record: SeqRecord, db: str = "engineered-core@1.0.0", detectors: Iterable[str] | None = None) -> List[Feature]:
    artifacts = manager.get_artifacts(db)
    return annotate_record(record, artifacts, detectors)


def score(
    record: SeqRecord,
    annotations: Sequence[Feature] | None = None,
    db: str = "engineered-core@1.0.0",
) -> Mapping[str, object]:
    artifacts = manager.get_artifacts(db)
    features = list(annotations) if annotations is not None else annotate_record(record, artifacts, None)
    return compute_score(record, features, artifacts)


def annotate_and_score(
    record: SeqRecord,
    db: str = "engineered-core@1.0.0",
    detectors: Iterable[str] | None = None,
) -> Mapping[str, object]:
    annotations = annotate(record, db=db, detectors=detectors)
    score_report = score(record, annotations=annotations, db=db)
    return {
        "sequence_id": record.id,
        "length": len(record.seq),
        "annotations": [feature.to_dict() for feature in annotations],
        "score": score_report,
        "db": db,
    }


set_cache_dir = manager.set_cache_dir
set_offline = manager.set_offline
add_registry = manager.add_registry
