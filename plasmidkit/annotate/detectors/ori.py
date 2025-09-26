from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    ori_entries = db.get("ori", [])
    for entry in ori_entries:
        motifs = entry.get("motifs", [])
        for motif in motifs:
            motif_hits = find_motifs(sequence, [motif])
            for pos in motif_hits:
                start = pos
                end = pos + len(motif)
                features.append(
                    Feature(
                        type="rep_origin",
                        id=entry["id"],
                        start=start,
                        end=end,
                        method="motif",
                        confidence=0.9,
                        evidence={"motif": motif, "position": pos},
                    )
                )
    return features
