# KMT2A-r Leukemia Paper

Notes for R codes for manuscript "**Single-cell multiomics reveals increased plasticity, resistant populations and stem-cell-like blasts in KMT2A-rearranged leukemia**"

## Process/integrate scRNA-Seq data
- script_process_scRNA_single_sample.R  -- process each sample from cellranger output
- integrate_18MLLr_scRNA.R  -- integrate scRNA-Seq data for all patient samples
- integrate_hd_scRNA.R  -- integrate scRNA-Seq data for all healthy donor samples

## Process/integrate scATAC-Seq data
- integrate_18MLLr_scATAC.R  -- integrate scRNA-Seq data for all patient samples
- integrate_hd_scATAC.R -- integrate scATAC-Seq data for all healthy donor samples

## Label transfering 
- transfer_lable.R -- Label transfering form scRNA-Seq data to scATAC-Seq data to annotate HD scATAC-Seq data

## Projection 
- project2norm_rna.R. -- project patient cells on top of umap of normal reference (HD) using scRNA-Seq data
- project2norm_atac.R -- project patient cells on top of umap of normal reference (HD) using scATAC-Seq data

## TF enrichment analysis
- script_chromVAR_allPeaks.R -- run chromVAR given a seurat object
- compTFenrich_18MLLr.R -- comparing TF enrichment between different groups using saved chromVAR outputs/objects

## Analysis of snmC-Seq2 data
- integrate_snmc.R -- integrate all snmc-seq2 data
- coembedding_snmc.R -- coembedding snmc-seq2 with corresponding scATAC-Seq data
- pseudo_bulk_snmc.R -- compare differential methylated regions between young vs old group by developmental stage

## Cell-Cell-Communication analysis
- prepInputs4CellphoneDB.R -- prepare inputs for running cellphoneDB
- compare_cell_interactions_HSPC1vsBlasts.R -- compare ligand-receptor pairs between HSPC1vsImmune and BlastvsImmnue


## Lineage switch analysis
- lineage_swithc_example.R  -- analysis of 2 patients whose leukemia underwent lineage switch upon treatment
- lineage_switch_analysis_updated.R -- overall analysis of lineage switch


## Long read analysis
- filter_alignqc_res.R -- filter alignqc result for downstream use
- sum_longread_final.R --- summarize detection rate using longGF & alignqc overlapped result


## Reference
Chen, C., Yu, W., Alikarami, F., Qiu, Q., Chen, C. H., Flournoy, J., ... & Tan, K. (2020).   [Single-cell multiomics reveals increased plasticity, resistant populations and stem-cell-like blasts in KMT2A-rearranged leukemia](https://www.biorxiv.org/content/10.1101/2020.12.06.413930v1)

