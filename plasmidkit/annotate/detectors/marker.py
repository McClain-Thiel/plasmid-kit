from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    marker_entries = db.get("markers", [])
    for entry in marker_entries:
        motifs = entry.get("motifs", [])
        for motif in motifs:
            for pos in find_motifs(sequence, [motif]):
                features.append(
                    Feature(
                        type="cds",
                        id=entry["id"],
                        start=pos,
                        end=pos + len(motif),
                        strand="+",
                        method="motif",
                        confidence=0.85,
                        evidence={"motif": motif, "position": pos},
                    )
                )
    return features
