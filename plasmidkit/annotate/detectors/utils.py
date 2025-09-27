from __future__ import annotations

from typing import List, Sequence


def find_motifs(sequence: str, motifs: Sequence[str]) -> List[int]:
    sequence_upper = sequence.upper()
    hits: List[int] = []
    for motif in motifs:
        motif_upper = motif.upper()
        start = 0
        while True:
            idx = sequence_upper.find(motif_upper, start)
            if idx == -1:
                break
            hits.append(idx)
            start = idx + 1
    return hits


def gc_content(sequence: str) -> float:
    sequence_upper = sequence.upper()
    if not sequence_upper:
        return 0.0
    gc = sum(1 for base in sequence_upper if base in {"G", "C"})
    return gc / len(sequence_upper)


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.upper().translate(complement)[::-1]
