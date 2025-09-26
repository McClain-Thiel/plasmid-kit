from __future__ import annotations

from typing import Dict

from .utils import gc_content


def analyse(sequence: str) -> Dict[str, float]:
    return {
        "length": float(len(sequence)),
        "gc": gc_content(sequence),
    }
