from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Mapping, Optional, Any

from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from .annotate import annotate_record, load_record
from .annotate.types import Feature
from .cache import manager
from .exporters import export_gff3, export_json, export_minimal_genbank

__all__ = [
    "load_record",
    "annotate",
    "analyze",
    "export_json",
    "export_gff3",
    "export_minimal_genbank",
    "set_cache_dir",
    "set_offline",
    "bootstrap_data",
]


def annotate(
    record: SeqRecord | str | Path,
    db: str = "engineered-core@1.0.0",
    detectors: Iterable[str] | None = None,
    is_sequence: Optional[bool] = None,
    skip_prodigal: bool = False,
) -> List[Feature]:
    """
    Annotates a plasmid sequence with features from the database and ORF prediction.
    Automatically resolves overlaps to prioritize specific markers over generic ORFs.
    """
    artifacts = manager.get_artifacts(db)
    raw_features = annotate_record(record, artifacts, detectors, is_sequence=is_sequence, skip_prodigal=skip_prodigal)
    return _resolve_overlaps(raw_features)


def analyze(
    record: SeqRecord | str | Path,
    db: str = "engineered-core@1.0.0",
    detectors: Iterable[str] | None = None,
    is_sequence: Optional[bool] = None,
    skip_prodigal: bool = False,
) -> Mapping[str, Any]:
    """
    Analyzes a plasmid sequence, producing a comprehensive report.
    
    Returns a dictionary containing:
    - sequence_id: ID of the sequence
    - length: Length of the sequence
    - gc_content: GC content percentage (0-100)
    - annotations: List of detected features (overlap resolved)
    - feature_counts: Summary counts of feature types
    - db: Database version used
    """
    normalized_record = record if isinstance(record, SeqRecord) else load_record(record, is_sequence=is_sequence)
    
    # Annotate (includes overlap resolution)
    annotations = annotate(
        normalized_record, db=db, detectors=detectors, skip_prodigal=skip_prodigal
    )
    
    # Calculate Statistics
    gc = round(gc_fraction(normalized_record.seq) * 100, 2)
    
    return {
        "sequence_id": normalized_record.id,
        "length": len(normalized_record.seq),
        "gc_content": gc,
        "annotations": [feature.to_dict() for feature in annotations],
        "feature_counts": _count_features(annotations),
        "db": db,
    }


def _count_features(features: List[Feature]) -> Mapping[str, int]:
    counts = {}
    for f in features:
        counts[f.type] = counts.get(f.type, 0) + 1
    return counts


def _resolve_overlaps(features: List[Feature]) -> List[Feature]:
    """
    Resolves overlaps between high-confidence specific features (markers, ORIs)
    and generic predicted ORFs (Prodigal CDS).
    
    If a generic CDS significantly overlaps (>50%) with a specific feature,
    the generic CDS is removed in favor of the specific annotation.
    """
    specific_types = {"marker", "rep_origin", "promoter", "terminator", "resistance_marker"}
    
    high_value_features = []
    generic_cds = []
    others = []

    for f in features:
        if f.type == "CDS" and f.method == "pyrodigal":
            generic_cds.append(f)
        elif f.type in specific_types or f.method != "pyrodigal":
            high_value_features.append(f)
        else:
            others.append(f)

    kept_cds = []
    
    for cds in generic_cds:
        is_clobbered = False
        cds_len = cds.end - cds.start
        if cds_len <= 0:
            continue
            
        for hv in high_value_features:
            # Check overlap
            overlap_len = max(0, min(cds.end, hv.end) - max(cds.start, hv.start))
            
            # Simple linear overlap fraction of the CDS
            overlap_fraction = overlap_len / cds_len
            
            # Threshold: if >50% of the CDS is explained by a specific marker, drop the CDS
            if overlap_fraction > 0.5:
                is_clobbered = True
                break
        
        if not is_clobbered:
            kept_cds.append(cds)
            
    # Return sorted by start position for consistency
    all_features = high_value_features + kept_cds + others
    all_features.sort(key=lambda f: f.start)
    
    return all_features


set_cache_dir = manager.set_cache_dir
set_offline = manager.set_offline

def bootstrap_data(cache_dir: Optional[str] = None, offline: Optional[bool] = None) -> Mapping[str, object]:
    """Prepare local caches and warm up the default database."""
    if cache_dir is not None:
        manager.set_cache_dir(cache_dir)
    if offline is not None:
        manager.set_offline(bool(offline))
    manager.ensure_cache_ready()
    return manager.get_artifacts("engineered-core@1.0.0")