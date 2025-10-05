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


def find_motifs_fuzzy_tagged(
    sequence: str,
    motifs: Sequence[str],
    max_mismatches: int = 2,
    circular: bool = True,
    include_rc: bool = True,
) -> List[Tuple[int, str, str, int]]:
    """Fuzzy motif search (Hamming distance) with optional reverse-complement scanning.

    Optimization: seed-and-verify using the pigeonhole principle.
    Split motif into (d+1) exact seeds (d = max_mismatches). Search each seed
    exactly (forward and reverse-complement) to get candidate starts, then
    verify full Hamming distance at those starts. This avoids full sliding.

    Returns a list of tuples: (start_position, original_motif, strand, mismatches)
    where strand is "+" for forward and "-" for reverse-complement matches.
    """
    sequence_upper = sequence.upper()
    if not sequence_upper or not motifs:
        return []
    search_space = sequence_upper + (sequence_upper if circular else "")
    seq_len = len(sequence_upper)

    def reverse_complement_local(seq: str) -> str:
        tbl = str.maketrans("ACGT", "TGCA")
        return seq.translate(tbl)[::-1]

    def hamming_at(space: str, start: int, pattern: str, limit: int) -> int | None:
        mismatches = 0
        plen = len(pattern)
        # Bounds check first
        if start < 0 or start + plen > len(space):
            return None
        for j in range(plen):
            if space[start + j] != pattern[j]:
                mismatches += 1
                if mismatches > limit:
                    return None
        return mismatches

    def iter_seed_hits(space: str, seed: str) -> List[int]:
        pos_list: List[int] = []
        start = 0
        while True:
            idx = space.find(seed, start)
            if idx == -1:
                break
            pos_list.append(idx)
            start = idx + 1
        return pos_list

    best: dict[Tuple[int, str, str], int] = {}

    # Cap mismatches to small integer for performance guarantees
    d = 0 if max_mismatches <= 0 else 1 if max_mismatches == 1 else 2

    for motif in motifs:
        if not motif:
            continue
        mu = motif.upper()
        mlen = len(mu)
        if mlen == 0 or mlen > len(search_space):
            continue

        # Choose number of seeds = d+1, with minimum seed length
        num_seeds = d + 1
        # Avoid too-short seeds which create many candidates
        min_seed_len = 8
        # Compute seed spans (roughly equal partitions)
        base = mlen // num_seeds
        extra = mlen % num_seeds
        seeds: List[Tuple[int, int]] = []  # list of (offset, length)
        offset = 0
        for i in range(num_seeds):
            seg_len = base + (1 if i < extra else 0)
            if seg_len <= 0:
                continue
            seeds.append((offset, seg_len))
            offset += seg_len
        # If any seed is too short, merge adjacent segments to reach min_seed_len where possible
        merged: List[Tuple[int, int]] = []
        i = 0
        while i < len(seeds):
            off, ln = seeds[i]
            while ln < min_seed_len and i + 1 < len(seeds):
                # merge with next
                n_off, n_ln = seeds[i + 1]
                ln = (n_off + n_ln) - off
                i += 1
            merged.append((off, ln))
            i += 1
        if merged:
            seeds = merged

        patterns: List[Tuple[str, str]] = [(mu, "+")]
        rc = reverse_complement_local(mu) if include_rc else None
        if rc is not None:
            patterns.append((rc, "-"))

        for pat, strand in patterns:
            for seed_off, seed_len in seeds:
                seed = pat[seed_off : seed_off + seed_len]
                if not seed:
                    continue
                for hit in iter_seed_hits(search_space, seed):
                    start_idx = hit - seed_off
                    if start_idx >= seq_len:
                        # Only record candidates starting in first copy for circular handling
                        continue
                    mm = hamming_at(search_space, start_idx, pat, d)
                    if mm is None:
                        continue
                    key = (start_idx, motif, strand)
                    prev = best.get(key)
                    if prev is None or mm < prev:
                        best[key] = mm

    hits: List[Tuple[int, str, str, int]] = [
        (pos, motif, strand, mm) for (pos, motif, strand), mm in best.items()
    ]
    hits.sort(key=lambda t: (t[3], -len(t[1]), t[0], t[2]))
    return hits


def gc_content(sequence: str) -> float:
    sequence_upper = sequence.upper()
    if not sequence_upper:
        return 0.0
    gc = sum(1 for base in sequence_upper if base in {"G", "C"})
    return gc / len(sequence_upper)


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.upper().translate(complement)[::-1]
