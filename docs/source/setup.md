# Setup

## Install

```bash
pip install genal-python
```

```python
import genal
```

`genal` requires Python ≥ 3.8.

## PLINK 2

Several features require a working **PLINK 2** executable (clumping, PRS, association tests, proxy search, allele-frequency updates).

```python
import genal

# Option 1: let genal download/install PLINK2 (requires internet)
genal.install_plink()

# Option 2: point genal to an existing plink2 executable
genal.set_plink("/path/to/plink2")  # or "/path/to/folder/containing/plink2"
```

## Configuration and cache locations

On first import, `genal` creates a configuration folder under:

- `~/.genal/` (see {py:data}`genal.constants.CONFIG_DIR`)
- `~/.genal/config.json` stores paths (PLINK 2, liftover, default genotype path, reference data folder)

Temporary PLINK outputs are written to:

- `./tmp_GENAL/` (relative to your current working directory)
- You can delete it with {py:func}`genal.delete_tmp`

### Changing where reference files are stored

If you want reference panels to live on a different disk (recommended on shared clusters), set the reference folder:

```python
from genal.tools import set_reference_folder

set_reference_folder("/path/to/genal_reference_data")
```

### Saving a default genotype path

Most methods that take a genotype `path=...` (PRS, extraction, association testing) save that path in `~/.genal/config.json` so you can omit it next time by passing `path=None`.

## Reference data (two different concepts)

`genal` uses reference data in two different ways; the same parameter name (`reference_panel`) is used in multiple methods but the *expected value differs by workflow*.

### 1) Build-only reference (for preprocessing / filling columns)

Used by: {py:meth}`genal.Geno.preprocess_data` and {py:meth}`genal.Geno.get_reference_panel`.

- Purpose: fill/validate **rsIDs, CHR/POS, and alleles**.
- Typical values: `reference_panel="37"` or `reference_panel="38"`.
- Also accepts: a path to a `.bim`/`.pvar` (no extension), or a DataFrame with `CHR, POS, SNP, A1, A2`.

### 2) LD reference panel (for LD operations and frequencies)

Used by: {py:meth}`genal.Geno.clump`, proxy search in {py:meth}`genal.Geno.prs` / {py:meth}`genal.Geno.query_outcome`, and {py:meth}`genal.Geno.update_eaf`.

- Purpose: compute LD and/or allele frequencies from a 1000G-like panel via PLINK.
- Typical values: `reference_panel="EUR_37"` (or `"AFR_38"`, etc.).
- Also accepts: a path to a custom PLINK fileset (bed/bim/fam or pgen/pvar/psam), without extension.
- Built-in panels are downloaded automatically on first use and cached under your configured reference folder.

## Offline usage

You can run fully offline once you have:

- a configured `plink2` path, and
- downloaded (or custom) reference panels already present locally.

Features that may require internet on first use include: `install_plink()`, built-in reference panel downloads, chain file downloads for liftover, gene mapping downloads for `filter_by_gene`, and GWAS Catalog queries.
