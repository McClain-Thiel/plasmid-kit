from __future__ import annotations

from pathlib import Path

import pytest


def test_analyze(tmp_path: Path) -> None:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plasmidkit as pk

    record = pk.load_record(Path("tests/data/pUC19.fasta"))
    
    # Test analyze API replaces score
    report = pk.analyze(record)
    
    assert report["sequence_id"]
    assert report["length"] > 0
    assert 0 <= report["gc_content"] <= 100
    assert len(report["annotations"]) > 0
    assert report["feature_counts"]
    
    exported = tmp_path / "report.json"
    pk.export_json(report, exported)
    assert exported.exists()


def _list_fasta_files() -> list[Path]:
    data_dir = Path(__file__).parent / "data"
    return sorted(data_dir.glob("*.fasta"))


@pytest.mark.parametrize("fasta_name", [p.name for p in _list_fasta_files()])
def test_benchmark_all_fastas(fasta_name: str, tmp_path: Path) -> None:
    import sys
    import time

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plasmidkit as pk

    fasta_path = Path("tests/data") / fasta_name

    record = pk.load_record(fasta_path)

    t0 = time.perf_counter()
    report = pk.analyze(record)
    t1 = time.perf_counter()

    annotations = report["annotations"]
    assert annotations, f"expected at least one annotation for {fasta_name}"
    
    total_s = t1 - t0

    # Persist per-sequence timing for later inspection
    out_path = tmp_path / f"{fasta_name}.timing.json"
    pk.export_json(
        {
            "sequence_id": record.id,
            "length": len(record.seq),
            "file": str(fasta_path),
            "timings": {
                "total_s": total_s,
            },
        },
        out_path,
    )
    assert out_path.exists()

    print(f"{fasta_name}: total {total_s:.3f}s")


def _csv_for_fasta(fasta: Path) -> Path:
    # Expected naming pattern: <stem>-annotations.csv
    return fasta.with_name(f"{fasta.stem}-annotations.csv")


def _paired_fastas_with_csvs() -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for fasta in _list_fasta_files():
        csv_path = _csv_for_fasta(fasta)
        if csv_path.exists():
            pairs.append((fasta.name, csv_path.name))
    return pairs


@pytest.mark.parametrize("fasta_name,csv_name", _paired_fastas_with_csvs())
def test_annotations_overlap_expected_csv(fasta_name: str, csv_name: str) -> None:
    import sys
    import csv

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plasmidkit as pk

    fasta_path = Path("tests/data") / fasta_name
    csv_path = Path("tests/data") / csv_name

    record = pk.load_record(fasta_path)
    anns = pk.annotate(record)

    by_type: dict[str, list] = {}
    for a in anns:
        by_type.setdefault(a.type.lower(), []).append(a)

    expected_rows = []
    with open(csv_path, newline="", encoding="utf8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            expected_rows.append(row)

    # Focus on a subset of high-confidence types present in our simple DB
    expected_focus = [
        r for r in expected_rows if r.get("Type", "").lower() in {"rep_origin", "promoter", "cds"}
    ]

    def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
        return not (a_end < b_start or b_end < a_start)

    hits = 0
    for row in expected_focus:
        row_type = row["Type"].lower()
        try:
            start = int(row["start location"]) - 1  # CSV likely 1-based
            end = int(row["end location"])  # end-inclusive in CSV
            if start > end:
                # Circular wrap; skip for this simple test
                continue
        except Exception:
            continue
            
        # ALLOW MARKERS TO MATCH EXPECTED CDS
        candidates = by_type.get(row_type, [])
        if row_type == "cds":
            candidates.extend(by_type.get("marker", []))
            candidates.extend(by_type.get("resistance_marker", []))
            
        if any(overlaps(a.start, a.end, start, end) for a in candidates):
            hits += 1

    # Require at least 70% of CSV rows in focus to be overlapped by at least one annotation
    required_hits = max(1, int(len(expected_focus) * 0.7))
    assert hits >= required_hits, f"found {hits}/{len(expected_focus)} expected features for {csv_name}, needed at least {required_hits}"


def test_psc101_requires_rep_origin_overlap() -> None:
    import sys
    import csv

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plasmidkit as pk

    fasta_path = Path("tests/data/pSC101.fasta")
    csv_path = Path("tests/data/pSC101-annotations.csv")

    record = pk.load_record(fasta_path)
    anns = pk.annotate(record)

    # Partition annotations by type
    rep_origins = [a for a in anns if a.type.lower() == "rep_origin"]

    # Load expected rep_origin rows
    expected_rows = []
    with open(csv_path, newline="", encoding="utf8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("Type", "").lower() == "rep_origin":
                expected_rows.append(row)

    def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
        return not (a_end < b_start or b_end < a_start)

    # Require at least one expected rep_origin to overlap any detected rep_origin
    hit = False
    for row in expected_rows:
        try:
            start = int(row["start location"]) - 1  # CSV likely 1-based
            end = int(row["end location"])  # end-inclusive in CSV
            if start > end:
                continue
        except Exception:
            continue
        if any(overlaps(a.start, a.end, start, end) for a in rep_origins):
            hit = True
            break

    assert hit, "Expected a pSC101 rep_origin overlapping the CSV interval"