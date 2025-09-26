from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Mapping

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .annotate.types import Feature


def export_json(data: Mapping[str, object], path: str | Path) -> None:
    with Path(path).open("w", encoding="utf8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def export_gff3(record: SeqRecord, annotations: Iterable[Feature], path: str | Path) -> None:
    lines = ["##gff-version 3"]
    seq_id = record.id or "sequence"
    for feature in annotations:
        attributes = f"ID={feature.id};Note={feature.type}"
        lines.append(
            "\t".join(
                [
                    seq_id,
                    feature.method,
                    feature.type,
                    str(feature.start + 1),
                    str(feature.end),
                    f"{feature.confidence:.2f}",
                    feature.strand,
                    ".",
                    attributes,
                ]
            )
        )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf8")


def export_minimal_genbank(record: SeqRecord, annotations: Iterable[Feature], path: str | Path) -> None:
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    seq = record[:]  # copy
    seq.features = []
    for feature in annotations:
        location = FeatureLocation(feature.start, feature.end, strand=1 if feature.strand == "+" else -1)
        seq.features.append(SeqFeature(location=location, type=feature.type, qualifiers={"label": feature.id}))
    with Path(path).open("w", encoding="utf8") as handle:
        SeqIO.write(seq, handle, "genbank")
