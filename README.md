# MDD Variant Prioritization Polygenic Risk Score Analysis

This repository contains code for a variant prioritization pipeline that uses functional annotations to construct partial polygenic risk scores (PRS) for major depressive disorder. This framework can be extended to integrate cross-trait information to improve prediction by leveraging shared genetic architecture across related psychiatric disorders. 

## Find your way around

| Folder | What it does |
|---|---|
| `prs_construction/` | Prepares GWAS inputs, runs PRS-CS and calculates individual PRS. |
| `elastic_net/` | Combines scores and evaluates the MDD-only, multi-trait and cross-trait models. |
| `figures/` | Plots the original PRS comparisons, SNP-set overlaps and ancestry-transfer results. |
| `config/` and `scripts/` | Set file locations, launch R scripts and run checks. |
| `docs/` | Explains inputs, software, execution order and reproducibility limits. |

The usual analysis order is: prepare inputs, calculate scores, combine models,
evaluate performance, then make figures and tables. You can run the part you
need once its required inputs are available.

## Before you start

You will need R and the packages used by your chosen analysis. The PRS-CS and
LDSC steps also need their Python environments; full genotype analyses generally
require a compute cluster. Details are in the [software guide](docs/SOFTWARE.md).

Scripts are numbered in execution order within each analysis folder. A matching `.R` script and `.sh` launcher perform the same step; use one execution route. See the [run-order guide](docs/WORKFLOW.md) for required arguments, optional steps and setup instructions.


