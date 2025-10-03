from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    terminator_entries = db.get("terminators", [])
    for entry in terminator_entries:
        motifs = entry.get("motifs", [])
        for motif in motifs:
            for pos in find_motifs(sequence, [motif], circular=True):
                features.append(
                    Feature(
                        type="terminator",
                        id=entry["id"],
                        start=pos,
                        end=pos + len(motif),
                        strand="+",
                        method="motif",
                        confidence=0.75,
                        evidence={"motif": motif, "position": pos},
                    )
                )
    return features
