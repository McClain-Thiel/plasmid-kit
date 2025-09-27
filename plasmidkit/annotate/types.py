from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class Feature:
    """Represents an annotated feature on a plasmid sequence."""

    type: str
    id: str
    start: int
    end: int
    strand: str = "+"
    method: str = "heuristic"
    confidence: float = 0.5
    evidence: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        data = {
            "type": self.type,
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "method": self.method,
            "confidence": self.confidence,
        }
        if self.evidence:
            data["evidence"] = self.evidence
        return data


def clamp_position(value: int, length: int) -> int:
    if value < 0:
        return 0
    if value > length:
        return length
    return value
