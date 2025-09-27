from __future__ import annotations

from collections import defaultdict
from typing import Dict

from .utils import reverse_complement


def analyse(sequence: str, k: int = 12) -> Dict[str, float]:
    seq = sequence.upper()
    counts: Dict[str, int] = defaultdict(int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        counts[kmer] += 1
    total_repeat_bases = sum((count - 1) * k for count in counts.values() if count > 1)
    longest_palindrome = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        rc = reverse_complement(kmer)
        if rc == kmer:
            longest_palindrome = max(longest_palindrome, k)
    return {
        "repeat_bases": float(total_repeat_bases),
        "longest_palindrome": float(longest_palindrome),
    }
