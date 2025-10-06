from __future__ import annotations

from typing import Dict, List

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    promoter_entries = db.get("promoters", [])
    MAX_MISMATCHES = 1
    # Collect hits across all entries first
    all_hits: list[tuple[int, str, str, int, str]] = []  # pos, motif, strand, mm, id
    for entry in promoter_entries:
        motifs = entry.get("motifs", [])
        entry_id = entry.get("id", "promoter")
        hits = find_motifs_fuzzy_tagged(
            sequence,
            motifs,
            max_mismatches=MAX_MISMATCHES,
            circular=True,
            include_rc=True,
        )
        for pos, motif, strand, mismatches in hits:
            all_hits.append((pos, motif, strand, mismatches, entry_id))

    # Prefer fewer mismatches, longer motifs, then position
    all_hits.sort(key=lambda t: (t[3], -len(t[1]), t[0]))

    # Greedy non-overlapping selection across entries; collapse strand duplicates
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
                type="promoter",
                id=entry_id,
                start=start,
                end=end,
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
        occupied.append(span)
        seen_spans.add(span)
    return features
