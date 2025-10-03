from __future__ import annotations

from typing import List, Sequence, Tuple

# Optional dependency: pyahocorasick
try:
    import ahocorasick as _ahocorasick  # type: ignore

    _HAS_AHOCORASICK = True
except Exception:
    _ahocorasick = None  # type: ignore
    _HAS_AHOCORASICK = False


def find_motifs(sequence: str, motifs: Sequence[str], circular: bool = True) -> List[int]:
    # Use naive fallback for small motif sets; pyahocorasick for large sets
    sequence_upper = sequence.upper()
    hits: List[int] = []
    if not sequence_upper:
        return hits
    search_space = sequence_upper + (sequence_upper if circular else "")
    seq_len = len(sequence_upper)

    if len(motifs) < 8:
        for motif in motifs:
            motif_upper = motif.upper()
            start = 0
            while True:
                idx = search_space.find(motif_upper, start)
                if idx == -1:
                    break
                if idx < seq_len:
                    hits.append(idx)
                start = idx + 1
        return sorted(set(hits))

    if not _HAS_AHOCORASICK:
        # Fallback to naive if automaton not available
        for motif in motifs:
            motif_upper = motif.upper()
            start = 0
            while True:
                idx = search_space.find(motif_upper, start)
                if idx == -1:
                    break
                if idx < seq_len:
                    hits.append(idx)
                start = idx + 1
        return sorted(set(hits))

    automaton = _ahocorasick.Automaton()
    for motif in motifs:
        motif_upper = motif.upper()
        automaton.add_word(motif_upper, motif_upper)
    automaton.make_automaton()
    for end_idx, _ in automaton.iter(search_space):
        mlen = 0  # compute matched length by checking all motif candidates; store max
        # pyahocorasick doesn't directly expose match length with default value; recheck set
        for motif in motifs:
            mu = motif.upper()
            if search_space.endswith(mu, 0, end_idx + 1):
                mlen = max(mlen, len(mu))
        start_idx = end_idx + 1 - mlen
        if start_idx < seq_len:
            hits.append(start_idx)
    return sorted(set(hits))


def find_motifs_tagged(sequence: str, motifs: Sequence[str], circular: bool = True) -> List[Tuple[int, str]]:
    sequence_upper = sequence.upper()
    if not sequence_upper or not motifs:
        return []
    search_space = sequence_upper + (sequence_upper if circular else "")
    seq_len = len(sequence_upper)

    # Naive path for small motif sets
    if len(motifs) < 8:
        hits: List[Tuple[int, str]] = []
        for motif in motifs:
            mu = motif.upper()
            start = 0
            while True:
                idx = search_space.find(mu, start)
                if idx == -1:
                    break
                if idx < seq_len:
                    hits.append((idx, motif))
                start = idx + 1
        # Deduplicate and sort
        return sorted(set(hits), key=lambda t: (t[0], t[1]))

    # Aho-Corasick for larger motif sets if available
    if not _HAS_AHOCORASICK:
        # Fallback to naive if automaton not available
        hits: List[Tuple[int, str]] = []
        for motif in motifs:
            mu = motif.upper()
            start = 0
            while True:
                idx = search_space.find(mu, start)
                if idx == -1:
                    break
                if idx < seq_len:
                    hits.append((idx, motif))
                start = idx + 1
        return sorted(set(hits), key=lambda t: (t[0], t[1]))

    automaton = _ahocorasick.Automaton()
    for motif in motifs:
        mu = motif.upper()
        automaton.add_word(mu, (mu, len(mu), motif))  # store upper, length, original
    automaton.make_automaton()

    tagged_hits: List[Tuple[int, str]] = []
    for end_idx, value in automaton.iter(search_space):
        mu, mlen, original = value
        start_idx = end_idx + 1 - mlen
        if start_idx < seq_len:
            tagged_hits.append((start_idx, original))
    return sorted(set(tagged_hits), key=lambda t: (t[0], t[1]))


def gc_content(sequence: str) -> float:
    sequence_upper = sequence.upper()
    if not sequence_upper:
        return 0.0
    gc = sum(1 for base in sequence_upper if base in {"G", "C"})
    return gc / len(sequence_upper)


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.upper().translate(complement)[::-1]
