from __future__ import annotations

from typing import Dict, Iterable, Mapping, Sequence

from Bio.SeqRecord import SeqRecord

from ..annotate.types import Feature
from . import rules


def compute_score(record: SeqRecord, annotations: Sequence[Feature], db: Mapping[str, object]) -> Dict[str, object]:
    synth = rules.synthesise_components(record, annotations, db)
    assembly = rules.assembly_components(record, annotations)
    components: Dict[str, float] = {}
    components.update(synth)
    components.update(assembly)
    total = rules.combine_components(components)
    return {
        "total": round(total, 2),
        "components": {key: round(value, 2) for key, value in components.items()},
    }
