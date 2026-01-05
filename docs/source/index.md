```{image} Images/genal_logo.png
:alt: genal logo
:width: 320px
:align: center
```

# genal

Polygenic risk scoring (PRS) and Mendelian randomization (MR) in Python.

**Quick links**
- Project docs: https://genal.readthedocs.io
- Paper (citation): https://doi.org/10.1093/bioadv/vbae207
- Source code: https://github.com/CypRiv/genal

```{toctree}
:maxdepth: 2
:caption: Documentation

introduction
setup
concepts
workflows
methods
api
faq
```

## Citation

If you use genal, please cite:

Cyprien A. Rivier, Santiago Clocchiatti-Tuozzo, Shufan Huo, Victor Torres-Lopez, Daniela Renedo, Kevin N. Sheth, Guido J. Falcone, Julian N. Acosta.  
**Genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization**.  
Bioinformatics Advances (2024). DOI: https://doi.org/10.1093/bioadv/vbae207

### Additional method citations

If you use the following tools/methods through `genal`, please also cite the original method/tool papers:

- **MR-PRESSO** (via {py:meth}`genal.Geno.MRpresso`): Marie Verbanck, Chia‑Yen Chen, Benjamin Neale, Ron Do. *Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases*. Nature Genetics (2018). DOI: https://doi.org/10.1038/s41588-018-0099-7. PMID: 29686387.
- **Colocalization (ABF)** (via {py:meth}`genal.Geno.colocalize`): Claudia Giambartolomei et al. *Bayesian test for colocalisation between pairs of genetic association studies using summary statistics*. PLoS Genetics (2014). DOI: https://doi.org/10.1371/journal.pgen.1004383.
- **PLINK 2** (used for clumping/PRS/association tests/allele frequencies): Chang C.C. et al. *Second-generation PLINK: rising to the challenge of larger and richer datasets*. GigaScience (2015). DOI: https://doi.org/10.1186/s13742-015-0047-8.
