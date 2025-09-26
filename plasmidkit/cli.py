from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import typer

from . import api
from .cache import manager
from .exporters import export_gff3, export_json, export_minimal_genbank

app = typer.Typer(help="PlasmidKit command line interface")


@app.command()
def annotate(
    input: Path = typer.Argument(..., help="Input FASTA/GenBank file"),
    db: str = typer.Option("engineered-core@1.0.0", help="Database identifier"),
    detectors: Optional[str] = typer.Option(None, help="Comma-separated detector list"),
    out_json: Optional[Path] = typer.Option(None, help="Write annotations+score JSON"),
    out_gff: Optional[Path] = typer.Option(None, help="Write annotations as GFF3"),
    out_gb: Optional[Path] = typer.Option(None, help="Write annotations as minimal GenBank"),
) -> None:
    record = api.load_record(input)
    detector_list = detectors.split(",") if detectors else None
    annotations = api.annotate(record, db=db, detectors=detector_list)
    result = {
        "sequence_id": record.id,
        "length": len(record.seq),
        "annotations": [feature.to_dict() for feature in annotations],
        "db": db,
    }
    if out_json:
        export_json(result, out_json)
    if out_gff:
        export_gff3(record, annotations, out_gff)
    if out_gb:
        export_minimal_genbank(record, annotations, out_gb)
    typer.echo(json.dumps(result, indent=2))


@app.command()
def score(
    input: Path = typer.Argument(..., help="Input FASTA/GenBank file"),
    db: str = typer.Option("engineered-core@1.0.0", help="Database identifier"),
    detectors: Optional[str] = typer.Option(None, help="Comma-separated detector list"),
    out_json: Optional[Path] = typer.Option(None, help="Write annotations+score JSON"),
) -> None:
    record = api.load_record(input)
    detector_list = detectors.split(",") if detectors else None
    annotations = api.annotate(record, db=db, detectors=detector_list)
    score = api.score(record, annotations=annotations, db=db)
    result = {
        "sequence_id": record.id,
        "length": len(record.seq),
        "annotations": [feature.to_dict() for feature in annotations],
        "score": score,
        "db": db,
    }
    if out_json:
        export_json(result, out_json)
    typer.echo(json.dumps(result, indent=2))


@app.command()
def fetch(db: str = typer.Argument(..., help="Database identifier")) -> None:
    manager.ensure_cache_ready()
    manager.get_artifacts(db)
    typer.echo(f"Fetched {db} into {manager.get_cache_dir()}")


@app.command()
def cache(action: str = typer.Argument(..., help="list or purge")) -> None:
    cache_dir = manager.get_cache_dir()
    if action == "list":
        files = sorted(cache_dir.glob("**/*"))
        for file in files:
            typer.echo(str(file))
    elif action == "purge":
        for file in cache_dir.glob("**/*"):
            if file.is_file():
                file.unlink()
        typer.echo("Cache cleared")
    else:
        raise typer.BadParameter("Action must be 'list' or 'purge'")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
