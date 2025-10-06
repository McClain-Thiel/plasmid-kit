from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    marker_entries = db.get("markers", [])
    # Collect hits across markers, then select non-overlapping globally
    all_hits: list[tuple[int, str, str, int, str]] = []  # pos, motif, strand, mm, id
    for entry in marker_entries:
        motifs = entry.get("motifs", [])
        entry_id = entry.get("id", "marker")
        hits = find_motifs_fuzzy_tagged(
            sequence,
            motifs,
            max_mismatches=1,
            circular=True,
            include_rc=True,
        )
        for pos, motif, strand, mismatches in hits:
            all_hits.append((pos, motif, strand, mismatches, entry_id))

    all_hits.sort(key=lambda t: (t[3], -len(t[1]), t[0]))
    occupied: list[tuple[int, int]] = []
    seen_spans: set[tuple[int, int]] = set()
    for pos, motif, strand, mismatches, entry_id in all_hits:
        start = pos
        end = pos + len(motif)
        span = (start, end)
        if span in seen_spans:
            continue
        if any(not (end <= s or start >= e) for s, e in occupied):
            continue
        features.append(
            Feature(
                type="cds",
                id=entry_id,
                start=start,
                end=end,
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
        occupied.append(span)
        seen_spans.add(span)
    return features
