# RNA-Seq Trajectory & Expression Analysis Pipeline

An end-to-end, modular R and Bash pipeline designed for comprehensive transcriptomic data analysis, gene expression trajectory clustering, and interactive visualization. 

This repository is a workflow of end-to-end analysis of a complex, multi time point, multi-condition RNASeq experiment. It covers everything from QC of raw data through differential expression tests with DESeq2 or sleuth, utilizes fuzzy clustering of trajectory-based DE data and Mahalanobis outliers, through functional enrichment and gene set enrichment analysis, and a deployed interactive Shiny application for exploratory data analysis.

## Pipeline Workflow

1. **QC & Quantification (`01`):** Shell-based processing and quality control for raw sequencing reads.
2. **Differential Expression (`02a`/`2b`):** Robust statistical modeling using `DESeq2` or `sleuth` to identify significant transcriptional changes.
3. **Clustering & Trajectory Analysis (`03` - `05`):** Grouping genes by expression patterns over trajectories, utilizing distance metrics like Mahalanobis distance to model relationships and track trends.
4. **Functional Enrichment (`04a` - `04b`, `06`):** Automated cluster profiling, Gene Ontology (GO), COG, and KEGG pathway enrichment frameworks.
5. **Advanced Visualization & Interactive Dashboard (`07` - `10`):** Additional, publication-ready static plots, custom faceting, interactive `plotly` implementations, and a dedicated `Shiny` web application for dynamic exploration.

## Example Figures
#### Volcano:
<img src=https://github.com/kiragoff/TranscriptomicTrajectoryToolkit/blob/main/figures/volcano-res_24h_exposure-qval_0.05-lfc_1.png>

#### Fuzzy cluster Mahalanobis outliers:
<img src = https://github.com/kiragoff/TranscriptomicTrajectoryToolkit/blob/main/figures/mahalanobis_0.95.png>

#### COG category enrichment hierarchical dotplot:
<img src = https://github.com/kiragoff/TranscriptomicTrajectoryToolkit/blob/main/figures/COG_cat_enrichment_m0.25-resc1.0_c9_hierarchical_dotplot.png>

#### Heatmap of enrichment of cyclic and branched chain amino acid GO processes:
<img src = https://github.com/kiragoff/TranscriptomicTrajectoryToolkit/blob/main/figures/figure_1_bcaa_and_cyclic-1.png>

#### cnetplot of enriched GO processes in a single cluster:
<img src = https://github.com/kiragoff/TranscriptomicTrajectoryToolkit/blob/main/figures/Cluster_8_BP_cnetplot.png>
