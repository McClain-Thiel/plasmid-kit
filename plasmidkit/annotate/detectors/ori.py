from __future__ import annotations

from typing import Dict, List, Tuple

from ..types import Feature
from .utils import find_motifs_fuzzy_tagged


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    features: List[Feature] = []
    ori_entries = db.get("ori", [])

    # Biological heuristics:
    # - Filter out trivially short motifs (e.g., single bases) that create many false positives
    # - Prefer the longest, non-overlapping matches per ori entry
    # - Optionally respect provided length_range to downweight/skip absurd hits

    MIN_MOTIF_LEN = 12  # typical ori sub-motifs (RNAI, AT-rich region) are longer than ~10bp
    MAX_MISMATCHES = 1  # keep ori tolerance tight (0-1 mismatches)

    for entry in ori_entries:
        raw_motifs: List[str] = entry.get("motifs", [])
        length_range: List[int] = entry.get("length_range", [])

        # Filter motifs by length
        motifs = [m for m in raw_motifs if isinstance(m, str) and len(m) >= MIN_MOTIF_LEN]
        if not motifs:
            continue

        # Find fuzzy tagged hits (position, motif, strand, mismatches)
        raw_hits = find_motifs_fuzzy_tagged(
            sequence,
            motifs,
            max_mismatches=MAX_MISMATCHES,
            circular=True,
            include_rc=True,
        )
        hits: List[Tuple[int, str, str, int]] = raw_hits
        if not hits:
            continue

        # Sort hits by fewest mismatches, then motif length desc, then position
        hits_sorted = sorted(hits, key=lambda t: (t[3], -len(t[1]), t[0]))

        # Greedily select non-overlapping longest matches
        selected: List[Tuple[int, str, str, int]] = []
        occupied: List[Tuple[int, int]] = []  # list of [start, end)
        for pos, motif, strand, mismatches in hits_sorted:
            start = pos
            end = pos + len(motif)
            overlaps = any(not (end <= s or start >= e) for s, e in occupied)
            if overlaps:
                continue
            # If length_range is provided, ensure candidate span size is plausible
            if length_range and len(motif) < min(length_range) * 0.05:  # motific chunk should be >=5% of ori span
                continue
            selected.append((pos, motif, strand, mismatches))
            occupied.append((start, end))

        # Emit features for selected hits only
        for pos, motif, strand, mismatches in selected:
            start = pos
            end = pos + len(motif)
            features.append(
                Feature(
                    type="rep_origin",
                    id=entry["id"],
                    start=start,
                    end=end,
                    strand=strand,
                    method="motif_fuzzy",
                    confidence=0.9,
                    evidence={
                        "motif": motif,
                        "position": pos,
                        "filtered_min_len": MIN_MOTIF_LEN,
                        "mismatches": mismatches,
                        "max_mismatches": MAX_MISMATCHES,
                    },
                )
            )

    return features
