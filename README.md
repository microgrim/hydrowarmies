# Active Autotrophs in Hydrothermal Vents Manuscript Analysis

This is a cleaned version of the repository containing only the essential code and outputs referenced in the Elkassas, Fortunato et al. manuscript.

## Repository Structure

- `data/` - Essential input data and metadata
- `output/` - Analysis results and data exports  
- `figures/` - Manuscript-ready figures
- `*.R` - Analysis scripts (numbered pipeline)

## Analysis Pipeline

Run scripts in numerical order:

**Environment Setup:**
1. `00_custom_functions.R` - Load custom functions
2. `00_resource_allocation.R` - Set up computational resources

**Main Analysis Pipeline:**
3. `01_DifferentialExpression.R` - Differential expression analysis (originally: 01_metat_coa_deseq.R)
4. `02_SDE_VennAnalysis.R` - SDE analysis and Venn diagrams (originally: 02_metat_coa_sde_venn.R)
5. `03_SDE_Visualization.R` - SDE visualization (originally: 03_metat_coa_sde_plots.R)
6. `04_CollateBins.R` - Collate genomic bins (originally: 04_metat_collate_bins.R)
7. `05_MAG_Transcriptomics.R` - MAG transcriptomics analysis (originally: 05_metat_mag_tx.R)
8. `06_GeneTaxonomy.R` - Gene taxonomy analysis (originally: 06_metat_gene_taxonomy.R)
9. `07_BinAbundances.R` - Bin abundance analysis (originally: 07_mag_bin_abundances.R)
10. `08_BinExpression.R` - Bin expression analysis (originally: 08_mag_bin_expression.R)
11. `09_ActiveMAG_Expression.R` - Active MAG expression (originally: 09_active_mag_expression.R)
12. `10_FinalFigures.R` - Generate final figures (originally: 10_coa_mag_bin_figures.R)

## Note

This repository has been cleaned of hard-coded paths. You may need to modify file paths in scripts to match your local directory structure.

