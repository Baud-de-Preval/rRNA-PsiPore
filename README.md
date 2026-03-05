# rRNA - PsiPore
Simple Nextflow pipeline to call RNA modifications via dorado basecaller models.

---

## Table of Contents

- [Overview] (#Overview)
- [Dependencies] (#dependencies)
- [Installation] (#installation)
- [Run] (#Run)
- [Citation] (#Citation)
- [Contact] (#Contact)
 
---

## Overview

This nextflow pipeline uses tools developped by nanopore sequencing to quantify pseudouridine stoichiometries on human ribosomal RNA.
It uses seqtagger with Apptainer to demultiplex aligned reads and modkit to assess the modification frequencies. GPU access is needed to run the pipeline for the basecalling and demplutiplexing steps.

---

## Dependencies

Please, make sure you have the latest version of dorado (https://github.com/nanoporetech/dorado) to access the up-to-date modification models.

---

## Installation

```bash

```
---

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
