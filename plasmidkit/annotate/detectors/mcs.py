from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    sites = db.get("mcs_sites", [])
    for site in sites:
        motif = site.get("sequence")
        if not motif:
            continue
        positions = find_motifs(sequence, [motif], circular=True)
        for pos in positions:
            features.append(
                Feature(
                    type="restriction_site",
                    id=site["id"],
                    start=pos,
                    end=pos + len(motif),
                    strand="+",
                    method="motif",
                    confidence=0.7,
                    evidence={"motif": motif, "copies": len(positions)},
                )
            )
    return features
