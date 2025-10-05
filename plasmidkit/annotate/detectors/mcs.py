from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    sites = db.get("mcs_sites", [])
    for site in sites:
        motif = site.get("sequence")
        if not motif:
            continue
        hits = find_motifs_fuzzy_tagged(
            sequence,
            [motif],
            max_mismatches=0,  # exact by default for restriction sites
            circular=True,
            include_rc=True,
        )
        for pos, _motif, strand, mismatches in hits:
            features.append(
                Feature(
                    type="restriction_site",
                    id=site["id"],
                    start=pos,
                    end=pos + len(motif),
                    strand=strand,
                    method="motif_fuzzy" if mismatches else "motif",
                    confidence=0.7,
                    evidence={"motif": motif, "mismatches": mismatches},
                )
            )
    return features
