from __future__ import annotations

from importlib import import_module
from typing import Dict, Iterable, List, Mapping

from ..types import Feature


_DEFAULT_ORDER = ["ori", "marker", "promoter", "terminator", "mcs", "orf_prodigal"]


def get_detector(name: str):
    module = import_module(f"plasmidkit.annotate.detectors.{name}")
    if not hasattr(module, "detect"):
        raise ValueError(f"Detector {name!r} does not define a detect() function")
    return module.detect


def run_detectors(sequence: str, db: Mapping[str, object], detectors: Iterable[str] | None = None) -> List[Feature]:
    order = list(detectors) if detectors else _DEFAULT_ORDER
    features: List[Feature] = []
    for name in order:
        detector_fn = get_detector(name)
        features.extend(detector_fn(sequence, db))
    return features
