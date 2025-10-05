from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    promoter_entries = db.get("promoters", [])
    MAX_MISMATCHES = 1
    for entry in promoter_entries:
        motifs = entry.get("motifs", [])
        hits = find_motifs_fuzzy_tagged(
            sequence,
            motifs,
            max_mismatches=MAX_MISMATCHES,
            circular=True,
            include_rc=True,
        )
        # Prefer fewer mismatches, longer motifs
        hits = sorted(hits, key=lambda t: (t[3], -len(t[1]), t[0]))
        for pos, motif, strand, mismatches in hits:
            features.append(
                Feature(
                    type="promoter",
                    id=entry["id"],
                    start=pos,
                    end=pos + len(motif),
                    strand=strand,
                    method="motif_fuzzy",
                    confidence=0.8,
                    evidence={
                        "motif": motif,
                        "position": pos,
                        "mismatches": mismatches,
                        "max_mismatches": MAX_MISMATCHES,
                    },
                )
            )
    return features
