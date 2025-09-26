# PlasmidKit

PlasmidKit is a Python library and CLI for annotating plasmid sequences and estimating a synthesis/assembly "makeability" score. The implementation focuses on fast local execution, a clear scoring breakdown, and extensibility through pluggable detectors and registries.

## Quick start

```bash
uv sync
uv run plasmidkit --help
```

Annotate and score a plasmid sequence in a single step:

```bash
uv run plasmidkit score tests/data/example.fasta --out-json report.json
```

## Features

- Pure Python API with optional extras for accelerated motif scanning.
- Deterministic annotation and scoring with transparent rule contributions.
- Offline-friendly cache layout with HuggingFace-style directories.
- CLI utilities for annotation, scoring, cache management, and data prefetch.
- Simple exporter helpers for JSON, GFF3, and minimal GenBank output.

## Library usage

```python
import plasmidkit as pk

record = pk.load_record("tests/data/example.fasta")
annotations = pk.annotate(record)
score = pk.score(record, annotations=annotations)

pk.export_json({"annotations": [feat.to_dict() for feat in annotations], "score": score}, "report.json")
```

## CLI usage

```
plasmidkit annotate input.fasta --out-json annotations.json
plasmidkit score input.fasta --out-json score.json
plasmidkit fetch engineered-core@1.0.0
plasmidkit cache list
```

## Development

Install dependencies with [`uv`](https://docs.astral.sh/uv/):

```bash
uv sync
uv run --with pytest python -m pytest
uv run ruff check .
```

## Data sources

The repository ships with a lightweight `engineered-core@1.0.0` dataset that covers common replicons, selectable markers, promoters, terminators, restriction sites, and forbidden motifs used by the built-in heuristics. Additional registries can be added at runtime with `plasmidkit.add_registry()`.
