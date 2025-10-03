from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    promoter_entries = db.get("promoters", [])
    for entry in promoter_entries:
        motifs = entry.get("motifs", [])
        for pos, motif in find_motifs_tagged(sequence, motifs, circular=True):
            features.append(
                Feature(
                    type="promoter",
                    id=entry["id"],
                    start=pos,
                    end=pos + len(motif),
                    strand="+",
                    method="motif",
                    confidence=0.8,
                    evidence={"motif": motif, "position": pos},
                )
            )
    return features
