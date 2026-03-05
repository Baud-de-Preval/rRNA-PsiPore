# rRNA - PsiPore

Simple Nextflow pipeline to call RNA modifications via dorado basecaller models.

## Overview

This Nextflow pipeline quantifies pseudouridine (Ψ) stoichiometries on human ribosomal RNA from Oxford Nanopore direct RNA sequencing data (pod5). It integrates SeqTagger for barcode demultiplexing and modkit for per-site modification frequency estimation, producing ready-to-use modification tables.

> ⚠️ GPU access is required for the basecalling and demultiplexing steps.

## Dependencies

Please, note this pipeline was run with the rna004_130bps_sup@v5.1.0 model. 
You can access more recent models with the latest version of dorado (see https://software-docs.nanoporetech.com/dorado/latest/models/list/).

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥ 25.0 | Workflow manager |
| Dorado | = 1.4.0 | Reads basecall |
| Apptainer | ≥ 1.3 | Containerization |
| samtools | = 1.22.1 | BAM processing |
| modkit | = 0.5.1 | Manipulation of bedMethyl files |

## Installation

You will need to install nextflow []

```bash

```

## Run

```bash
nextflow run main.nf -with-conda --multiplex true
```
---

## Citation

If you use this pipeline, please cite:

> Author et al. (2026). *Nanopore Direct RNA Sequencing Enables Reproducible, Site-Resolved Pseudouridine Quantification in Human Ribosomal RNA *. Journal, volume(issue), pages. https://doi.org/xxxxx

Also cite the underlying tool:
- **Seqtagger**: 

---

## Contact

Email: Baudouin.Seguineau.De.Preval@USherbrooke.ca
---
