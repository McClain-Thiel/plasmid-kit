# PlasmidKit

[![Tests](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/tests.yml/badge.svg)](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/tests.yml)
[![PyPI Publish](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/publish.yml/badge.svg)](https://github.com/McClain-Thiel/plasmid-kit/actions/workflows/publish.yml)
[![PyPI version](https://img.shields.io/pypi/v/plasmidkit)](https://pypi.org/project/plasmidkit/)

PlasmidKit is a fast Python library and CLI for annotating plasmid sequences and estimating a synthesis/assembly "makeability" score. It focuses on engineered plasmids (2–10 kb typical), validates that backbone essentials are present, and produces an interpretable quality score.

## Quick Start

```bash
# Install with uv
uv sync

# Run the CLI
uv run plasmidkit --help

# Or install with pip
pip install -e .
```

## What It Evaluates

**Backbone Components:**
- **Origins of replication** (e.g., ColE1/pMB1, p15A, pSC101, RSF)
- **Selectable markers** (e.g., AmpR/blaTEM, KanR/nptII, CmR/cat)
- **Promoters** (e.g., lac, T7, CMV) and **terminators**
- **ORFs** via Prodigal to confirm coding potential exists

**Synthesis & Assembly Hygiene:**
- Length optimization (2–6 kb ideal)
- GC content (45–55% ideal)
- Repeat sequences and palindromes
- Homopolymer runs
- Forbidden motifs (e.g., BsaI, BsmBI, NotI)

## Example Usage

### Python API

```python
import json
import plasmidkit as pk

# Load a plasmid sequence
record = pk.load_record("tests/data/pUC19.fasta")

# Annotate features
annotations = pk.annotate(record)

# Fast annotation (skip ORF prediction)
annotations_fast = pk.annotate(record, skip_prodigal=True)

# Annotate a raw sequence string (must specify is_sequence=True)
raw_seq = "ATCG..."
annotations_raw = pk.annotate(raw_seq, is_sequence=True)

# Calculate quality score
score_report = pk.score(record, annotations=annotations)

# Display results
print(f"Found {len(annotations)} features")
print(f"Overall score: {score_report['total']:.1f}/100")

# Show first few annotations
for ann in annotations[:3]:
    print(f"  {ann.type}: {ann.id} at {ann.start}-{ann.end}")
```

### Real Output (pUC19)

**Annotations (first 5):**
```json
[
  {
    "type": "rep_origin",
    "id": "ColE1",
    "start": 2314,
    "end": 2903,
    "strand": "+",
    "method": "motif_fuzzy",
    "confidence": 1.0
  },
  {
    "type": "rep_origin",
    "id": "pBR322_origin",
    "start": 2298,
    "end": 2917,
    "strand": "-",
    "method": "motif_fuzzy",
    "confidence": 1.0
  },
  {
    "type": "marker",
    "id": "TEM-116",
    "start": 1283,
    "end": 2144,
    "strand": "+",
    "method": "motif_fuzzy",
    "confidence": 1.0
  },
  {
    "type": "promoter",
    "id": "AmpR_promoter_4",
    "start": 1178,
    "end": 1283,
    "strand": "+",
    "method": "motif_fuzzy",
    "confidence": 1.0
  },
  {
    "type": "promoter",
    "id": "lac_promoter",
    "start": 540,
    "end": 569,
    "strand": "-",
    "method": "motif_fuzzy",
    "confidence": 1.0
  }
]
```

**Score Report:**
```json
{
  "total": 22.46,
  "components": {
    "length": 15.0,
    "gc": 10.0,
    "repeats": -0.54,
    "palindromes": -2.0,
    "homopolymers": -0.0,
    "BsaI": -1.0,
    "BsmBI": -1.0,
    "forbidden_motifs": -2.0,
    "ori_recognition": 8.0,
    "marker_recognition": -8.0,
    "promoter_terminator": 4.0,
    "burden": 0.0
  }
}
```

### Command Line

```bash
# Annotate a plasmid
uv run plasmidkit annotate tests/data/pUC19.fasta

# Skip slow ORF prediction
uv run plasmidkit annotate tests/data/pUC19.fasta --skip-prodigal

# Get detailed output with scores
uv run plasmidkit annotate tests/data/pUC19.fasta --output report.json

# View all detected features
uv run plasmidkit annotate tests/data/pUC19.fasta --verbose
```

## How It Works

1. **DNA Motif Matching**: Uses `pyahocorasick` for fast multi-pattern scanning of known sequences
   - Supports fuzzy matching with configurable mismatches
   - Handles circular plasmids correctly (wrap-around search)

2. **ORF Prediction**: Runs Prodigal (`pyrodigal`) to identify coding regions
   - Validates that the plasmid has protein-coding potential
   - No exhaustive protein identification needed

3. **Sequence Analysis**: Evaluates synthesis/assembly properties
   - GC content, length optimization
   - Detects repeats, palindromes, homopolymers
   - Flags forbidden restriction sites

4. **Scoring System**: 
   ```
   Total Score = Synthesis Hygiene + Backbone Recognition + Assembly Burden
   ```
   - **Synthesis**: length (15), GC (10), minus penalties for repeats/palindromes/homopolymers/forbidden motifs
   - **Backbone**: ori (8), marker (6), promoter/terminator (4-7)
   - **Burden**: penalties for high-copy + strong promoter combinations

## Feature Detection

PlasmidKit identifies:

| Feature Type | Detection Method | Examples |
|-------------|------------------|----------|
| **rep_origin** | Motif matching | ColE1, pMB1, p15A, pSC101, RSF1030 |
| **marker** | Motif matching + aliases | AmpR (TEM-1), KanR (nptII), CmR (cat) |
| **promoter** | Motif matching | lac, T7, CMV, AmpR promoter |
| **terminator** | Motif matching | rrnB T1, T7 terminator |
| **cds** | Prodigal ORF prediction | Any plausible coding sequence |

## Data Sources

We curate signatures from public sources, with per-entry citations in `plasmidkit/data/engineered_core_signatures.json`:

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

# Specific test
uv run pytest tests/test_api.py::test_annotate_and_score -v
```

**Test coverage includes:**
- Annotation accuracy across multiple plasmids (pUC19, pSC101, etc.)
- Score calculation and component breakdown
- Feature detection for ori, markers, promoters, terminators
- Edge cases and circular sequence handling

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

## Data Caching

Large database files are not stored in git. They're cached under `plasmidkit/data/_cache/` and ignored by `.gitignore`.

- Default cache directory can be overridden with `PLASMIDKIT_CACHE`
- To prefetch caches for offline use:

```bash
uv run python -m plasmidkit.cli bootstrap --cache-dir plasmidkit/data/_cache
```

This warms up the built-in `engineered-core@1.0.0` database. Optional external indices (BLAST/Rfam/SnapGene/SwissProt) can be placed under the cache dir if available.

## Design Philosophy

**Focus on Backbone Recognition:**
- PlasmidKit targets engineered backbone essentials rather than exhaustive protein identification
- Motif-based detection for known functional elements (fast, interpretable)
- ORF prediction confirms coding potential exists (no protein ID needed)

**CDS Detection Strategy:**
1. **Motif matches** for common selectable markers (presence implies function)
2. **ORF prediction** (via Prodigal/pyrodigal) to confirm plausible coding regions

**Why not exhaustive annotation?**
- Engineered plasmids have predictable backbones
- Exact/fuzzy motif matching is fast and accurate for known parts
- Full protein BLAST is slow and often unnecessary for quality assessment

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
