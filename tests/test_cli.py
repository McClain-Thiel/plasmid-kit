from pathlib import Path
from typer.testing import CliRunner
from plasmidkit.cli import app
import json

runner = CliRunner()

def test_analyze_cli(tmp_path: Path) -> None:
    # Use absolute path or relative to project root
    fasta_path = Path("tests/data/pUC19.fasta")
    json_path = tmp_path / "output.json"
    
    result = runner.invoke(app, ["analyze", str(fasta_path), "--out-json", str(json_path)])
    
    assert result.exit_code == 0
    assert "sequence_id" in result.stdout
    # pUC19 ID in the file is "Addgene_NGS_Result"
    assert "Addgene_NGS_Result" in result.stdout
    assert json_path.exists()
    
    with open(json_path) as f:
        data = json.load(f)
        assert data["length"] > 0
        assert "annotations" in data

def test_annotate_cli() -> None:
    fasta_path = Path("tests/data/pUC19.fasta")
    
    result = runner.invoke(app, ["annotate", str(fasta_path)])
    assert result.exit_code == 0
    # Annotate outputs a JSON list
    data = json.loads(result.stdout)
    assert isinstance(data, list)
    assert len(data) > 0
    # Check for a known feature
    assert any(f["id"] == "TEM-116" for f in data)

def test_cache_list_cli() -> None:
    result = runner.invoke(app, ["cache", "list"])
    assert result.exit_code == 0
