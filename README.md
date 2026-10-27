# Melon spatial and single-nucleus transcriptomics

Analysis code for *A new cell domain forming a metabolic envelope for seeds in cucurbits*.

Software and package versions are listed in [SOFTWARE_VERSIONS.md](SOFTWARE_VERSIONS.md).

## Data Source

The sequencing data are deposited under NCBI BioProject **PRJNA1511485**.

## Content

```text
01_spatial_transcriptomics/
  00_spaceranger_count.sh
  01_ST_analysis.R
02_single_nucleus/
  00_DNBC4tools_mixseq.sh
  01_snRNA_analysis.R
  02_DEG_GO_Augur/
    01_DEG_identification.R
    02_GO_enrichment.R
    Augur.R
  03_reclustering_states/
    Parenchyma/
      01_reclustering.R
      02_hdWGCNA.R
      MO_monocle.R
      C_monocle.R
      P_monocle.R
      03_prepare_velocity_h5ad.R
      04_scvelo.py
      05_CytoTRACE.R
    OEC/
      01_reclustering.R
      02_hdWGCNA.R
      MO_monocle.R
      C_monocle.R
      P_monocle.R
      03_SCODE.sh
03_integration/
  01_stereoscope.sh
05_cross_species_synteny/
  01_tomato_celltype_correlation.R
06_bulk_RNA-seq/
  01_alignment_quantification.sh
  02_differential_expression.R
```
