from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged, calculate_motif_confidence


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    terminator_entries = db.get("terminators", [])
    all_hits: list[tuple[int, str, str, int, str]] = []
    for entry in terminator_entries:
        motifs = entry.get("motifs", [])
        entry_id = entry.get("id", "terminator")
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
        confidence = calculate_motif_confidence(len(motif), mismatches, max_mismatches=1)
        features.append(
            Feature(
                type="terminator",
                id=entry_id,
                start=start,
                end=end,
                strand=strand,
                method="motif_fuzzy",
                confidence=confidence,
                evidence={
                    "motif": motif,
                    "position": pos,
                    "mismatches": mismatches,
                    "pct_identity": round((len(motif) - mismatches) / len(motif) * 100, 2),
                },
            )
        )
        occupied.append(span)
        seen_spans.add(span)
    return features
