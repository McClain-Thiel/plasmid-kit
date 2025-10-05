from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    marker_entries = db.get("markers", [])
    for entry in marker_entries:
        motifs = entry.get("motifs", [])
        hits = find_motifs_fuzzy_tagged(
            sequence,
            motifs,
            max_mismatches=1,
            circular=True,
            include_rc=True,
        )
        hits = sorted(hits, key=lambda t: (t[3], -len(t[1]), t[0]))
        for pos, motif, strand, mismatches in hits:
            features.append(
                Feature(
                    type="cds",
                    id=entry["id"],
                    start=pos,
                    end=pos + len(motif),
                    strand=strand,
                    method="motif_fuzzy",
                    confidence=0.85,
                    evidence={
                        "motif": motif,
                        "position": pos,
                        "mismatches": mismatches,
                    },
                )
            )
    return features
