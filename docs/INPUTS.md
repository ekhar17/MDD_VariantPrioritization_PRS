# Required inputs

These inputs are not distributed in this repository. Set their locations in
the ignored `config/paths.local.sh` file using `config/paths.example.sh`.

## Participant-level inputs

- `MDD_BACKING`: the original `bigsnpr` genotype `.rds` descriptor with its
  corresponding `.bk` file. Retain the original sample and variant order.
- `MDD_SAMPLE`: its UK Biobank `.sample` file. Original readers expect the
  placeholder row and align `ID_1`/`ID_2` with `FID`/`IID`.
- `$MDD_PROJECT/01.TrainTestSplit/data/`: original QC cohort/covariate files.
  Readers use `FID`, `IID`, `MDD`, `pc1`-`pc8`, `sex`, `age`,
  `assessment_center` and `genotype_axiom` as appropriate.
- The original `splits/` subdirectory, with
  `split1_QC_MDD_alltraining.txt`, `split1_QC_MDD_testing.txt`, through split 10.
- Original AFR/SAS files such as
  `AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt` and the
  corresponding SAS file.

Use the original split memberships. This release does not regenerate new QC
cohorts, ancestry assignments or split files. Keep participant-level inputs and
generated participant predictions outside the public upload.

## Discovery statistics and prioritized variant lists

`MDD_SOURCE_DATA` retains the original `original_data/` and
`derived_data/procedure_data/` layout. Input basenames identify the original
discovery versions. The MDD file is
`daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz`.

The original MDD prioritization table is
`derived_data/procedure_data/ipsychMDDGWAS.sumstats.QC.maf01info09_5grp.txt`.
It supplies original coordinates/alleles and `cr95set`, `logbf_over4` and
`twassig` memberships. Related-trait tables use the original `fgannosumstats/`
layout. Published discovery files and the original post-MetaSubtract ANX
statistics must be obtained from their original sources.

Initial fGWAS/TWAS prioritization and the original MetaSubtract code were not
present in the supplied source collection; their saved outputs are inputs to
this core release. This repository does not replace those original analyses.

`MDD_HAPMAP` points to `snpinfo_1kg_hm3`; `MDD_LD_REF` points to the corresponding
European 1000 Genomes LD blocks for PRS-CS.

## Existing result paths

Data/output filenames retain their historical names under
`$MDD_PROJECT/05.PRS_CS_construction/` and `$MDD_PROJECT/06.Elastic_Net/`.
Those external data names are separate from the numbered code in
`prs_construction/` and `elastic_net/`.

The cross-trait data folders retain the original spellings
`muti_trait_MDDcoef_alttraitset` and `muti_trait_MDDset_alttraitcoef`, corresponding
to the code folders `elastic_net/cross_trait_variant_sets/` and
`elastic_net/cross_trait_effect_sizes/`.

Original plotting intermediates `auc_eur_res.csv` and
`nagr2_diff_eur_res.csv` under `05.PRS_CS_construction/data/` are required by
the AUC and risk-metric plots. Their generators were not found. The optional
historical coefficient script also reads the original `OLD` and `P_version_PRSCS`
comparison datasets.
