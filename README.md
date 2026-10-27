# Melon spatial and single-nucleus transcriptomics code

Code-availability package for the manuscript *A new cell domain forming a metabolic envelope for seeds in cucurbits*.

The repository contains only the core input paths, key parameters and a few essential analysis commands. Figure layout, color definitions, repeated plotting code and large data objects are intentionally omitted.

## Data source

The manuscript reports the raw sequencing and spatial-metabolomics data under NCBI BioProject **PRJNA1511485**. Raw FASTQ files, filtered matrices, Seurat objects, annotation tables and vendor reports are not stored in GitHub. See `00_metadata/data_sources.tsv` and `input_requirements.tsv`.

## Contents

- `01_spatial_transcriptomics/`: data source and unresolved Space Ranger inputs.
- `02_single_nucleus/`: the supplied DNBC4tools command and the minimal Seurat core workflow.
- `03_integration/`: integration method and missing Stereoscope inputs.
- `04_deg_go_augur/`: minimal GO-enrichment command.
- `05_parenchyma_states/` to `11_validation/`: concise method/data-source notes and unresolved-code records.
- `00_metadata/`: provenance, software versions, input requirements and missing-code notes.

The package is deliberately not presented as a complete end-to-end reproduction pipeline. Missing commands are recorded rather than reconstructed.

## Local execution

Replace the example paths in the scripts with the corresponding deposited data paths. The snRNA script expects the `MO`, `C`, `P` and `U` filtered count-matrix directories. The DNBC4tools script requires the study-specific `melon_v4` STAR index and `gene.gtf`.

## License

No license has been selected yet. Choose a license before making the repository public.
