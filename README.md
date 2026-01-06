# PlasmidKit

[![Tests](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/tests.yml/badge.svg)](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/tests.yml)
[![PyPI Publish](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/publish.yml/badge.svg)](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/publish.yml)
[![PyPI version](https://img.shields.io/pypi/v/plasmidkit)](https://pypi.org/project/plasmidkit/)
[![Documentation](https://img.shields.io/badge/docs-live-blue)](https://mcclain-thiel.github.io/plasmid-kit/)

**Fast Python library and CLI for comprehensive plasmid annotation and analysis.**

## Design Philosophy

Most plasmid annotation tools suffer from two major problems:
1. They are just CLI wrappers around "weird" bioinformatics tools that are hard to install and use.
2. They require you to bring your own database, leaving the burden of curation on the user.

**PlasmidKit is different.** It comes with a built-in, highly curated database of over 1,000 origins of replication (DoriC) and 7,000+ resistance markers (AMRFinderPlus). It uses fast, native Python implementations for motif matching and analysis, ensuring it works out of the box with zero configuration.

**Key Features:**
- **Curated Database:** Includes DoriC 12.1 and NCBI AMRFinderPlus data.
- **Intelligent Annotation:** Prioritizes high-confidence markers over generic open reading frames (ORFs).
- **Fast:** Optimized for 2–20 kb engineered plasmids.
- **Easy to Use:** Simple Python API and CLI.

## Quick Start

```bash
# Install with uv
uv sync

# Run the CLI
uv run plasmidkit --help

# Or install with pip
pip install -e .
```

## Example Usage

### Python API

```python
import json
import plasmidkit as pk

# Load a plasmid sequence
record = pk.load_record("tests/data/pUC19.fasta")

# Analyze plasmid (returns a report dictionary)
report = pk.analyze(record)

# Display results
print(f"Sequence: {report['sequence_id']}")
print(f"Length: {report['length']} bp")
print(f"GC Content: {report['gc_content']}%")
print(f"Features: {report['feature_counts']}")

# Show first few annotations
for ann in report['annotations'][:3]:
    print(f"  {ann['type']}: {ann['id']} at {ann['start']}-{ann['end']}")
```

### Command Line

```bash
# Analyze a plasmid and get a report
uv run plasmidkit analyze tests/data/pUC19.fasta

# Save full analysis to JSON
uv run plasmidkit analyze tests/data/pUC19.fasta --out-json report.json

# Annotate and export to GFF3
uv run plasmidkit annotate tests/data/pUC19.fasta --out-gff output.gff
```

## Feature Detection

PlasmidKit identifies:

| Feature Type | Detection Method | Examples |
|-------------|------------------|----------|
| **rep_origin** | Motif matching (DoriC) | ColE1, pMB1, p15A, pSC101, RSF1030 |
| **marker** | Motif matching (AMRFinder) | AmpR (TEM-1), KanR (nptII), CmR (cat) |
| **promoter** | Motif matching | lac, T7, CMV, AmpR promoter |
| **terminator** | Motif matching | rrnB T1, T7 terminator |
| **cds** | Prodigal ORF prediction | Any plausible coding sequence (clobbered by specific markers if overlapping) |

## Data Sources

We curate signatures from public sources, with per-entry citations in `plasmidkit/data/engineered_core_signatures.json`:

- **DoriC 12.1** (Database of Replication Origins):
  - Over 1,000 plasmid origins of replication
  - Source: http://tubic.tju.edu.cn/doric/

- **NCBI AMRFinderPlus**:
  - Comprehensive database of antimicrobial resistance genes
  - Includes over 6,000 resistance markers grouped by gene family
  - Source: https://github.com/ncbi/amr

- **PlasMapper** features API (promoters/terminators/origins)
  - Portal: https://plasmapper.wishartlab.com/search
  - API: `https://plasmapper.ca/api/features`

- **NCBI** ori sequences via query:
  - `origin_of_replication[All Fields] AND (bacteria[filter] AND plasmid[filter])`
  - Citation: `https://www.ncbi.nlm.nih.gov/nuccore/<ACCESSION>`

- **UniProt (Swiss-Prot)** for reviewed markers:
  - Examples: blaTEM-1 (P62593), nptII (P00552)
  - API: `https://rest.uniprot.org/uniprotkb/{accession}`

- **CARD** (Comprehensive Antibiotic Resistance Database):
  - Protein homolog models for bacterial AMR determinants
  - PHM entries: beta-lactamases, aminoglycoside-modifying enzymes, etc.

- **SnapGene** Standard Features export:
  - Engineered backbone motifs: promoters, terminators, origins, markers
  - Short DNA motifs for fast exact/fuzzy matching
  - Citation: `{ "database": "SnapGene", "source": "Standard Features export" }`

- **pLannotate** bundle indices:
  - SnapGene/FPbase/Swiss-Prot indices
  - Rfam models for RNA features

## Testing

PlasmidKit includes comprehensive tests that run automatically via GitHub Actions on every push.

**Run tests locally:**
```bash
# All tests
uv run pytest tests/ -v

# With coverage
uv run pytest tests/ --cov=plasmidkit --cov-report=html
```

## Development

```bash
# Clone the repository
git clone https://github.com/McClain-Thiel/plasmid-kit.git
cd plasmid-kit

# Install with development dependencies
uv sync

# Run tests
uv run pytest

# Format code
uv run black plasmidkit/ tests/
uv run ruff check plasmidkit/ tests/
```

## License

See LICENSE file for details.

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Citation

If you use PlasmidKit in your research, please cite:
```
[Citation information to be added]
```