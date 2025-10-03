from __future__ import annotations

from typing import Dict, List, Tuple

from ..types import Feature


START_CODONS = ("ATG",)
STOP_CODONS = ("TAA", "TAG", "TGA")


def _find_orfs_linear(seq: str, min_aa: int = 60) -> List[Tuple[int, int, str]]:
    hits: List[Tuple[int, int, str]] = []
    n = len(seq)
    for frame in range(3):
        i = frame
        while i + 3 <= n:
            codon = seq[i : i + 3]
            if codon in START_CODONS:
                j = i + 3
                while j + 3 <= n:
                    stop = seq[j : j + 3]
                    if stop in STOP_CODONS:
                        length_aa = (j + 3 - i) // 3 - 1
                        if length_aa >= min_aa:
                            hits.append((i, j + 3, "+"))
                        i = j + 3
                        break
                    j += 3
                else:
                    # no stop found; advance to end
                    i = n
                    continue
            else:
                i += 3
    return hits


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    # Simple, fast ORF finder on plus strand; circular wrap handled by detectors orchestrator if needed
    seq = sequence.upper()
    features: List[Feature] = []
    min_aa = int(db.get("orf_min_aa", 60)) if isinstance(db, dict) else 60
    for start, end, strand in _find_orfs_linear(seq, min_aa=min_aa):
        features.append(
            Feature(
                type="CDS",
                id=f"orf_{start}_{end}",
                start=start,
                end=end,
                strand=strand,
                method="orf",
                confidence=0.5,
                evidence={"min_aa": min_aa},
            )
        )
    return features


