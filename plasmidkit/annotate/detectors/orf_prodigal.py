from __future__ import annotations

from typing import Dict, List

from ..types import Feature


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    try:
        import pyrodigal
    except Exception:  # pragma: no cover - optional dependency not installed
        return []

    # Pyrodigal recommends passing bytes for performance
    seq_bytes = sequence.upper().encode("ascii", "ignore")
    min_len_aa = int(db.get("orf_min_aa", 60)) if isinstance(db, dict) else 60

    # Use Prodigal in metagenomic mode for plasmids lacking training signal; prefer bacterial code 11.
    # Support multiple pyrodigal versions by trying both parameter names.
    try:
        orf_finder = pyrodigal.GeneFinder(meta=True, closed=False, genetic_code=11)
    except TypeError:
        try:
            orf_finder = pyrodigal.GeneFinder(meta=True, closed=False, translation_table=11)
        except TypeError:
            orf_finder = pyrodigal.GeneFinder(meta=True, closed=False)
    genes = orf_finder.find_genes(seq_bytes)

    features: List[Feature] = []
    for g in genes:
        # Derive coordinates robustly across pyrodigal versions
        start = getattr(g, "begin", getattr(g, "start", None))
        end = getattr(g, "end", getattr(g, "stop", None))
        if start is None or end is None:
            continue
        nt_length = int(end) - int(start)
        length_aa = nt_length // 3
        if length_aa < min_len_aa:
            continue
        strand_val = getattr(g, "strand", 1)
        strand = "+" if strand_val == 1 else "-"
        features.append(
            Feature(
                type="CDS",
                id=f"prodigal_orf_{start}_{end}",
                start=int(start),
                end=int(end),
                strand=strand,
                method="pyrodigal",
                confidence=0.7,
                evidence={"min_aa": min_len_aa},
            )
        )
    return features


