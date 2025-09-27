from __future__ import annotations

from collections import Counter
from typing import Dict, Mapping, Sequence

from Bio.SeqRecord import SeqRecord

from ..annotate.detectors import gc_length, homopolymers, repeats
from ..annotate.types import Feature


def _length_score(length: float) -> float:
    if 2000 <= length <= 6000:
        return 15.0
    if 1000 <= length < 2000:
        return 8.0 + (length - 1000) / 1000 * 7.0
    if 6000 < length <= 10000:
        return 15.0 - (length - 6000) / 4000 * 7.0
    if length < 1000:
        return max(0.0, length / 1000 * 8.0)
    if length > 10000:
        over = min(length - 10000, 5000)
        return max(0.0, 8.0 - over / 5000 * 8.0)
    return 0.0


def _gc_score(gc_fraction: float) -> float:
    if 0.45 <= gc_fraction <= 0.55:
        return 10.0
    if gc_fraction <= 0.3 or gc_fraction >= 0.7:
        return 0.0
    diff = min(abs(gc_fraction - 0.5), 0.2)
    penalty = (diff / 0.2) ** 2 * 10.0
    return max(0.0, 10.0 - penalty)


def _repeat_penalty(repeat_bases: float, length: float) -> float:
    if length == 0:
        return 0.0
    ratio = repeat_bases / length
    return min(12.0, ratio * 12.0 * 5.0)


def _palindrome_penalty(longest: float) -> float:
    if longest < 12:
        return 0.0
    excess = longest - 12
    return min(8.0, 2.0 + excess / 6.0)


def _homopolymer_penalty(count: float) -> float:
    return min(6.0, count * 2.0)


def _forbidden_penalty(sequence: str, motifs: Sequence[Mapping[str, str]]) -> Dict[str, float]:
    penalties: Dict[str, float] = {}
    seq_upper = sequence.upper()
    total = 0.0
    for entry in motifs:
        motif = entry.get("sequence", "").upper()
        if not motif:
            continue
        occurrences = seq_upper.count(motif)
        if occurrences:
            penalty = min(9.0 - total, float(occurrences))
            total += penalty
            penalties[entry.get("id", motif)] = -penalty
        if total >= 9.0:
            break
    if total < 9.0:
        penalties["forbidden_motifs"] = -total
    return penalties


def synthesise_components(record: SeqRecord, annotations: Sequence[Feature], db: Mapping[str, object]) -> Dict[str, float]:
    sequence = str(record.seq)
    stats = gc_length.analyse(sequence)
    repeat_stats = repeats.analyse(sequence)
    homopolymer_stats = homopolymers.analyse(sequence)

    components: Dict[str, float] = {}
    components["length"] = round(_length_score(stats["length"]), 2)
    components["gc"] = round(_gc_score(stats["gc"]), 2)
    components["repeats"] = round(-_repeat_penalty(repeat_stats["repeat_bases"], stats["length"]), 2)
    components["palindromes"] = round(-_palindrome_penalty(repeat_stats["longest_palindrome"]), 2)
    components["homopolymers"] = round(-_homopolymer_penalty(homopolymer_stats["count"]), 2)

    forbidden = _forbidden_penalty(sequence, db.get("forbidden_motifs", []))
    for key, value in forbidden.items():
        components[key] = round(value, 2)
    return components


def _score_range(total: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, total))


def _ori_marker_scores(annotations: Sequence[Feature]) -> Dict[str, float]:
    components: Dict[str, float] = {}
    ori_present = any(feature.type == "rep_origin" for feature in annotations)
    marker_present = any(feature.type == "cds" for feature in annotations)
    components["ori_recognition"] = 8.0 if ori_present else -10.0
    components["marker_recognition"] = 6.0 if marker_present else -8.0
    promoter_present = any(feature.type == "promoter" for feature in annotations)
    terminator_present = any(feature.type == "terminator" for feature in annotations)
    promoter_score = 4.0 if promoter_present else 0.0
    terminator_score = 3.0 if terminator_present else 0.0
    components["promoter_terminator"] = promoter_score + terminator_score
    return components


def _mcs_score(annotations: Sequence[Feature]) -> float:
    restriction_sites = [feature.id for feature in annotations if feature.type == "restriction_site"]
    counts = Counter(restriction_sites)
    unique_sites = sum(1 for count in counts.values() if count == 1)
    capped = min(unique_sites, 3)
    return capped * 2.0


def _burden_penalty(length: float, annotations: Sequence[Feature]) -> float:
    ori_ids = {feature.id for feature in annotations if feature.type == "rep_origin"}
    promoter_ids = {feature.id for feature in annotations if feature.type == "promoter"}
    payload = length - sum(feature.end - feature.start for feature in annotations if feature.type == "cds")
    penalty = 0.0
    high_copy = bool(ori_ids & {"ColE1", "pMB1", "p15A"})
    strong_promoter = bool(promoter_ids & {"T7", "CMV"})
    if high_copy and strong_promoter and length > 5000:
        penalty -= 3.0
    if payload > 6000:
        penalty -= 2.0
    return penalty


def assembly_components(record: SeqRecord, annotations: Sequence[Feature]) -> Dict[str, float]:
    sequence_length = float(len(record.seq))
    components = _ori_marker_scores(annotations)
    components["mcs_uniqueness"] = _mcs_score(annotations)
    components["burden"] = _burden_penalty(sequence_length, annotations)
    return components


def combine_components(components: Mapping[str, float]) -> float:
    total = sum(components.values())
    return _score_range(total, 0.0, 100.0)
