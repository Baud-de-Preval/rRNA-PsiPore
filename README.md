# rRNA - PsiPore

Simple Nextflow pipeline to call RNA modifications via dorado basecaller models.

## Overview

This Nextflow pipeline quantifies pseudouridine (Ψ) stoichiometries on human ribosomal RNA from Oxford Nanopore direct RNA sequencing data (pod5). It integrates SeqTagger for barcode demultiplexing and modkit for per-site modification frequency estimation, producing ready-to-use modification tables.

> ⚠️ GPU access is required for the basecalling and demultiplexing steps.

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥ 25.0 | Workflow manager |
| Dorado | = 1.4.0 | Reads basecall |
| Apptainer | ≥ 1.3 | Containerization |
| samtools | = 1.22.1 | BAM processing |
| modkit | = 0.5.1 | Manipulation of bedMethyl files |

Make sure you have [Nextflow](https://www.nextflow.io/docs/latest/install.html), [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation) installed.
 
## Run

Please, note this pipeline was run with the rna004_130bps_sup@v5.1.0 model. You can access more recent models with the latest version of dorado (see https://software-docs.nanoporetech.com/dorado/latest/models/list/).

The sample names can be renamed and adapted in `data/sample_sheet.csv` with the corresponding number of barcode (up to 96).

```bash
nextflow run main.nf -with-conda --multiplex true
```
---

## Citation

---

**Contact**: Baudouin.Seguineau.De.Preval@USherbrooke.ca
