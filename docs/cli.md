# Command Line Interface

PlasmidKit provides a CLI tool `plasmidkit` for quick analysis and annotation.

## Commands

### `analyze`

Analyze a plasmid sequence and produce a comprehensive report.

```bash
plasmidkit analyze [OPTIONS] INPUT
```

**Arguments:**

*   `INPUT`: Input FASTA/GenBank file path. [Required]

**Options:**

*   `--db TEXT`: Database identifier. [default: engineered-core@1.0.0]
*   `--detectors TEXT`: Comma-separated list of detectors to use (e.g. `ori,marker`).
*   `--skip-prodigal`: Skip slow ORF detection (Prodigal).
*   `--out-json PATH`: Write full analysis report to a JSON file.
*   `--help`: Show this message and exit.

**Example Output:**

```json
{
  "sequence_id": "pUC19",
  "length": 2686,
  "gc_content": 50.63,
  "annotations": [ ... ],
  "feature_counts": { "marker": 1, "rep_origin": 1 },
  "db": "engineered-core@1.0.0"
}
```

### `annotate`

Annotate plasmid features and write to various standard formats.

```bash
plasmidkit annotate [OPTIONS] INPUT
```

**Arguments:**

*   `INPUT`: Input FASTA/GenBank file path. [Required]

**Options:**

*   `--db TEXT`: Database identifier. [default: engineered-core@1.0.0]
*   `--detectors TEXT`: Comma-separated list of detectors to use.
*   `--skip-prodigal`: Skip slow ORF detection.
*   `--out-json PATH`: Write annotations list to JSON.
*   `--out-gff PATH`: Write annotations as GFF3.
*   `--out-gb PATH`: Write annotations as minimal GenBank.
*   `--help`: Show this message and exit.

### `cache`

Manage the internal database cache.

```bash
plasmidkit cache [list|purge]
```

### `fetch`

Fetch a database artifact into the cache.

```bash
plasmidkit fetch DB_NAME
```

### `bootstrap`

Initialize the cache and load default data.

```bash
plasmidkit bootstrap [OPTIONS]
```
