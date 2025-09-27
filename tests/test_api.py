from __future__ import annotations

from pathlib import Path


def test_annotate_and_score(tmp_path: Path) -> None:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import plasmidkit as pk

    record = pk.load_record(Path("tests/data/example.fasta"))
    annotations = pk.annotate(record)
    assert annotations, "expected at least one annotation"
    report = pk.score(record, annotations=annotations)
    assert 0 <= report["total"] <= 100
    exported = tmp_path / "report.json"
    pk.export_json({"annotations": [feat.to_dict() for feat in annotations], "score": report}, exported)
    assert exported.exists()
