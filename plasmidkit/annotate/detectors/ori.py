from __future__ import annotations

from typing import Dict, List, Tuple

from ..types import Feature
from .utils import find_motifs_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    ori_entries = db.get("ori", [])

    # Biological heuristics:
    # - Filter out trivially short motifs (e.g., single bases) that create many false positives
    # - Prefer the longest, non-overlapping matches per ori entry
    # - Optionally respect provided length_range to downweight/skip absurd hits

    MIN_MOTIF_LEN = 12  # typical ori sub-motifs (RNAI, AT-rich region) are longer than ~10bp

    for entry in ori_entries:
        raw_motifs: List[str] = entry.get("motifs", [])
        length_range: List[int] = entry.get("length_range", [])

        # Filter motifs by length
        motifs = [m for m in raw_motifs if isinstance(m, str) and len(m) >= MIN_MOTIF_LEN]
        if not motifs:
            continue

        # Find tagged hits (position, motif)
        hits: List[Tuple[int, str]] = find_motifs_tagged(sequence, motifs, circular=True)
        if not hits:
            continue

        # Sort hits by motif length desc, then stable by position to prefer longer motifs
        hits_sorted = sorted(hits, key=lambda t: (-len(t[1]), t[0]))

        # Greedily select non-overlapping longest matches
        selected: List[Tuple[int, str]] = []
        occupied: List[Tuple[int, int]] = []  # list of [start, end)
        for pos, motif in hits_sorted:
            start = pos
            end = pos + len(motif)
            overlaps = any(not (end <= s or start >= e) for s, e in occupied)
            if overlaps:
                continue
            # If length_range is provided, ensure candidate span size is plausible
            if length_range and len(motif) < min(length_range) * 0.05:  # motific chunk should be >=5% of ori span
                continue
            selected.append((pos, motif))
            occupied.append((start, end))

        # Emit features for selected hits only
        for pos, motif in selected:
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
                    evidence={"motif": motif, "position": pos, "filtered_min_len": MIN_MOTIF_LEN},
                )
            )

    return features
