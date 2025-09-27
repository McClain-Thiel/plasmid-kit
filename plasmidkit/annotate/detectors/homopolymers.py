from __future__ import annotations

from typing import Dict


def analyse(sequence: str, min_run: int = 8) -> Dict[str, float]:
    seq = sequence.upper()
    longest = 0
    total = 0
    if not seq:
        return {"longest": 0.0, "count": 0.0}
    current = seq[0]
    run = 1
    for char in seq[1:]:
        if char == current:
            run += 1
        else:
            if run >= min_run:
                longest = max(longest, run)
                total += 1
            current = char
            run = 1
    if run >= min_run:
        longest = max(longest, run)
        total += 1
    return {"longest": float(longest), "count": float(total)}
