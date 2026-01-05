# FAQ / Troubleshooting

## PLINK path is not set

- Run {py:func}`genal.install_plink` (downloads PLINK 2), or
- Run {py:func}`genal.set_plink` to point to an existing `plink2` executable.

## `preprocess_data()` can’t find rsIDs / many SNPs are missing

Most commonly this is a **build mismatch**.

- If your GWAS is GRCh38/hg38, use `reference_panel="38"` for preprocessing.
- If your GWAS is GRCh37/hg19, use `reference_panel="37"`.
- If you need to convert, use {py:meth}`genal.Geno.lift` before preprocessing.

## PRS uses far fewer SNPs than expected

Common causes:

- Variant ID mismatch between summary stats and genotype dataset.
  - If you have `CHR`/`POS`, `genal` can map to genotype IDs before extraction.
  - Otherwise, consider curating SNP IDs or using proxying.
- Missing variants in genotypes (e.g., limited imputation).
  - Use `proxy=True` with an ancestry-matched LD reference panel.

## Clumping returns no SNPs

- Check the `P` column values (numeric, within $[0,1]$).
- Relax `p1` temporarily (e.g., `p1=1e-4`) to validate the pipeline.
- Ensure the chosen LD reference panel matches your variants (build + ancestry).

## MR results look wrong for binary outcomes

- Ensure effects are on the **log-odds** scale.
- If your GWAS reports odds ratios, call `preprocess_data(effect_column="OR")` to log-transform.
- Use `odds=True` in {py:meth}`genal.Geno.MR` only to *format* odds ratios from MR betas.

## Palindromic SNPs are removed/flipped

Control this with `action` in MR/MR-PRESSO:

- `action=1`: keep palindromes (assumes aligned)
- `action=2`: frequency-based (default)
- `action=3`: drop all palindromes

## GWAS Catalog queries fail or time out

- This is expected for many SNPs (catalog coverage is incomplete and rate limits exist).
- Consider increasing `timeout` or limiting `max_associations`.
- For very large SNP lists, annotate a subset.

## `tmp_GENAL/` is large

`genal` does not automatically clean PLINK temporary outputs (useful for debugging).

- Delete it with {py:func}`genal.delete_tmp`, or remove the folder manually.

## PLINK “out of memory”

`genal` passes `--memory` based on `G.ram`. If PLINK reports an OOM error:

```python
G.ram = 25_000  # MB
```

Then rerun the command.
