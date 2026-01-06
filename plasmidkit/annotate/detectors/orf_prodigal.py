from __future__ import annotations

from typing import Dict, List, Optional
import pyrodigal

from ..types import Feature


def detect(sequence: str, db: Dict[str, object]) -> List[Feature]:
    seq_bytes = sequence.upper().encode("ascii", "ignore")
    
    min_len_aa = int(db.get("orf_min_aa", 50)) if isinstance(db, dict) else 50
    min_len_nt = int(db.get("orf_min_nt", 150)) if isinstance(db, dict) else 150

    # Attempt to initialize GeneFinder with compatible arguments across versions
    try:
        orf_finder = pyrodigal.GeneFinder(
            meta=True,
            closed=False,
            genetic_code=11,
            min_gene=min_len_nt,
        )
    except TypeError:
        try:
            orf_finder = pyrodigal.GeneFinder(
                meta=True,
                closed=False,
                translation_table=11,
                min_gene=min_len_nt,
            )
        except TypeError:
            try:
                orf_finder = pyrodigal.GeneFinder(meta=True, closed=False, min_gene=min_len_nt)
            except TypeError:
                orf_finder = pyrodigal.GeneFinder(meta=True, closed=False)
                
    genes = orf_finder.find_genes(seq_bytes)

    features: List[Feature] = []
    for g in genes:
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
        
        # Extract Prodigal raw score if available
        raw_score = getattr(g, "score", None)
        
        evidence = {
            "min_aa": min_len_aa, 
            "min_nt": min_len_nt,
            "prodigal_score": raw_score
        }

        features.append(
            Feature(
                type="CDS",
                id=f"prodigal_orf_{start}_{end}",
                start=int(start),
                end=int(end),
                strand=strand,
                method="pyrodigal",
                confidence=None,
                evidence=evidence,
            )
        )
    return features