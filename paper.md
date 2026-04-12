---
title: 'SPADE: An R package for contamination-aware Spatial Proximity Analysis of Differential Expression'
tags:
  - R
  - spatial transcriptomics
  - bioinformatics
  - single-cell analysis
  - differential expression
authors:
  - name: Mingke Wu
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Immunology, School of Basic Medical Sciences, Peking University, Beijing, 100191, China
    index: 1
date: 2026-04-12
bibliography: paper.bib
---

# Summary

Spatial transcriptomics has offered us a remarkable tool to investigate the intricate cellular ecosystems within tissues, capturing the crucial gene expression at spatial level and providing unique insights about cell-cell interaction. Yet, despite the growing power of these technologies, a critical gap remains. Most existing tools are unable to systematically analyze gene expression in one cell type that is far or near another reference cell type, and the inherent complexity of spatial transcriptomics technology further compound this challenge: each spatial spot often contains mixed cell types regardless of whether we use reference-based label transfer or other deconvolution methods for annotation, leading to significant gene expression contamination from neighboring cells. This contamination will cover up true spatial expression patterns and makes it difficult to accurately identify genes that are up- or down-regulated with increasing distance from the reference cells. Consequently, researchers are unable to answer how nearby cells respond to a tumor boundary, an immune infiltration, or a stem cell niche from the level of spatial transcriptomics.

We developed SPADE (Spatial Proximity Analysis of Differential Expression), an R package built on the Seurat framework, to address this gap. SPADE calculates distances between query cells and reference cells, classifies cells into different states based on user-defined distance thresholds, and performs differential expression analysis between different states. Importantly, the package implements a control cell filtering strategy to remove contaminated signals and integrates functional enrichment analysis for biological interpretation.

Using a lung cancer Visium dataset, we demonstrated SPADE’s ability in discovering spatial-related different expressed genes while removing contaminated signals, and we successfully identified an antiviral activation and metabolic suppression program in epi-like tumor cells near myeloid cells.

SPADE provides a straightforward and reproducible protocol for distance-dependent differential expression analysis in spatial transcriptomics data. The package is open-source, freely available at https://github.com/sorrymaker03/SPADE.

Keywords: spatial transcriptomics; differential expression analysis; R package; cell-cell interaction

# Statement of need

Spatial transcriptomics technology has revolutionized biomedical research by integrating gene expression profiling with spatial localization and enabling the decoding of cellular communication networks within their surrounding microenvironments, offering an unprecedented resolution to investigate tumor invasion and immune infiltration[1][2].

Despite these advances, a critical analytical gap remains both in the technology and analysis level. Most existing computational tools lack the capacity to systematically quantify how gene expression changes in various distances[3][4]. For example, understanding how stromal cells transcriptionally adapt at 50 μm versus 200 μm from a tumor boundary is fundamental to studying metastatic niches. The challenge is further compounded by the inherent technical limitations of Spatial technologies. Unlike single-cell resolution, most spatial spots capture 5–20 cells per capture area[5]. In addition, most existing annotation methods, including reference-based label transfer or deconvolution approaches, will regard each spatial spot as a mixture of multiple cell types, even when one cell type constitutes most of the spot[6][7]. This leads to pervasive gene expression contamination, where signals from adjacent cells influence the expression of an ideal cell type. Consequently, when attempting to compare spatially resolved differential expression, we couldn’t observe the real biological shifts as response to distance[8][9].

To address this, we developed SPADE (Spatial Proximity Analysis of Differential Expression), an R package that integrates with Seurat objects to identify differentially expressed genes between cells located near versus far from a reference cell type, while filtering out potential contamination from adjacent cells. Here, we demonstrate SPADE’s capability using a human lung cancer spatial transcriptomics dataset, revealing that tumor cells proximal to myeloid cells exhibit a antiviral activation and metabolic quiescence program.

# Software design

## SPADE: Spatial Proximity Analysis of Differential Expression

We developed SPADE, an R package designed to identify truly differentially expressed genes between cells in close proximity versus distant from a reference cell type. SPADE integrates seamlessly with Seurat objects and leverages both spatial coordinates and cell-type annotations to perform systematic comparisons.

Cell selection and distance calculation: For each pair of target cell types (cell A vs cell B), SPADE identifies cells of each type based on metadata annotations. Euclidean distances from each cell A to the nearest cell B are computed using the get.knnx function from the FNN package[10].

Grouping by spatial proximity: Cell A are first classified into near, mid, and far groups based on user-defined distance thresholds relative to Cell B. Differential expressed genes between the near and far groups will be initially identified using Seurat’s FindMarkers function (log fold-change threshold = 0.25, minimum detection fraction = 10%)[11]. To remove potential confounding effects from neighboring Cell B, SPADE selects the closest Cell B for each Cell A as a control cluster. Genes that are differentially expressed between this control cluster and the far group of Cell A are then excluded from the initial near-vs-far list, effectively filtering signals that may result from contamination from Cell B genes. This two-step filtering ensures that the final list of differential genes reflects true spatially associated differences.

Functional enrichment analysis: Significantly upregulated and downregulated genes (adjusted p-value < 0.05) are mapped from gene symbols to Entrez IDs using clusterProfiler::bitr. Gene Ontology (GO) biological process enrichment is performed with clusterProfiler::enrichGO. Gene set enrichment analysis (GSEA) is performed on the ranked gene list using clusterProfiler::gseGO[12].

## Visualization

Spatial distribution of distances and groups (near, mid, and far) are visualized using SpatialFeaturePlot and SpatialDimPlot in Seurat[6]. Color gradients indicate distance from the reference cell type, and distinct colors highlight spatial proximity groups.

## Implementation and usage

SPADE is implemented in R (version ≥ 4.2) and depends on Seurat, FNN, cluster-Profiler, and org.Hs.eg.db. The package provides a single wrapper function, run_spade_full, that automates distance calculation, differential expression, functional enrichment, and spatial visualization. SPADE is open-source and available at https://github.com/sorrymaker03/SPADE.

The main function run_spade_full requires a spatial transcriptomics object (e.g., Seurat object with spatial coordinates), along with user-defined cell type groups and distance thresholds. In this workflow, SPADE automatically computes distances between specified cell populations, classifies cells into spatial proximity groups, performs differential expression analysis, applies control-based filtering, and outputs both table results and spatial visualizations.

## Spatial Transcriptomics Demo Data

Demo spatial transcriptomics dataset is derived from a human lung cancer sample (Neuroendocrine Carcinoma), provided by 10x Genomics (https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard). Spatial transcriptomics datasets were preprocessed using standard pipelines.

## Results

We successfully demonstrate SPADE’s ability in a spatial transcriptomics dataset derived from a human lung cancer sample. To explore the impact of spatially adjacent myeloid cells on tumor cells, we utilized SPADE to analyze epi-like lung tumor cells located proximal and distal to myeloid cells.

The gene expression in nearby spatial clusters would significantly contaminate their adjacent cells, as epi-like tumor cells near myeloid cells have a high expression of LYZ and IFIT3, which is the marker gene of myeloid cells. After SPADE’s process, the results of differential expressed genes, GO enrichment, and GSEA analysis collectively revealed a highly consistent and functionally transcriptional reprogramming features.

The most significant finding is that tumor cells located close to myeloid cells initiated a strong antiviral immune program. Except for two related genes IFI16 and LGALS3BP were significantly upregulated, enrichment analysis clearly showed that the entire interferon response pathway (including type I interferon signaling and response to virus) was systematically activated in proximal tumor cells. Simultaneously, pathways related to antigen processing and presentation were also significantly enriched.

In parallel with this, proximal tumor cells almost completely down-regulate energy metabolism and stress response programs. Enrichment analysis revealed that energy metabolism pathways such as oxidative phosphorylation and glycolysis were highly enriched in tumor cells distant from myeloid cells, while they were significantly inhibited in proximal tumor cells.

In summary, lung tumor cells located near myeloid cells undergo a state shift to metabolically quiescent and immune alert.

![Figure 1](./figure1.png)

Figure 1. (A) Volcano plot indicates that differential expressed genes would be contaminated by adjacent cells. (B) Spatial plot shows the demo distance and group information in SPADE. (C) Volcano plot demonstrates the true differential expressed genes between tumor cells close to or far from myeloid cells. (D) GO analysis shows the up and down regulated terms. (E) GSEA analysis shows the up and down regulated terms.

# Research impact statement

In this study, we introduced SPADE as a powerful tool for spatial proximity analysis of differential expression. Unlike conventional approaches that simply compare predefined anatomical regions, SPADE will automatedly calculate distance to a reference cell type and applies a control-based filtering step to remove genes potentially contributed by adjacent cells.

Our demonstration confirms that SPADE effectively eliminates myeloid marker genes (e.g., LYZ, IFIT3) from the final differential list.

Several limitations should be acknowledged. The current implementation uses Euclidean distance based on spatial coordinates, which assumes isotropic tissue architecture.

Nevertheless, SPADE provides a straightforward and effective solution for hypothesis-driven spatial analysis.


# AI usage disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions from Mingke Wu during the genesis of this project.

# References
1.  Jin, Y.; Zuo, Y.; Li, G.; Liu, W.; Pan, Y.; Fan, T.; Fu, X.; Yao, X.; Peng, Y. Advances in Spatial Transcriptomics and Its Applica-tions in Cancer Research. Mol Cancer 2024, 23, 129, doi:10.1186/s12943-024-02040-9.
2. 	Hunter, M.V.; Moncada, R.; Weiss, J.M.; Yanai, I.; White, R.M. Spatially Resolved Transcriptomics Reveals the Architecture of the Tumor-Microenvironment Interface. Nat Commun 2021, 12, 6278, doi:10.1038/s41467-021-26614-z.
3. 	Wang, Q.; Zhu, H.; Deng, L.; Xu, S.; Xie, W.; Li, M.; Wang, R.; Tie, L.; Zhan, L.; Yu, G. Spatial Transcriptomics: Biotechnologies, Computational Tools, and Neuroscience Applications. Small Methods 2025, 9, 2401107, doi:10.1002/smtd.202401107.
4. 	Long, K.A.; Lister, A.; Jones, M.R.W.; Adamski, N.M.; Ellis, R.E.; Chedid, C.; Carpenter, S.J.; Liu, X.; Backhaus, A.E.; Goldson, A.; et al. Spatial Transcriptomics Reveals Expression Gradients in Developing Wheat Inflorescences at Cellular Resolution. The Plant Cell 2026, 38, koaf282, doi:10.1093/plcell/koaf282.
5. 	Ståhl, P.L.; Salmén, F.; Vickovic, S.; Lundmark, A.; Navarro, J.F.; Magnusson, J.; Giacomello, S.; Asp, M.; Westholm, J.O.; Huss, M.; et al. Visualization and Analysis of Gene Expression in Tissue Sections by Spatial Transcriptomics. Science 2016, 353, 78–82, doi:10.1126/science.aaf2403.
6. 	Hao, Y.; Stuart, T.; Kowalski, M.H.; Choudhary, S.; Hoffman, P.; Hartman, A.; Srivastava, A.; Molla, G.; Madad, S.; Fernan-dez-Granda, C.; et al. Dictionary Learning for Integrative, Multimodal and Scalable Single-Cell Analysis. Nat Biotechnol 2024, 42, 293–304, doi:10.1038/s41587-023-01767-y.
7. 	Ma, Y.; Zhou, X. Spatially Informed Cell-Type Deconvolution for Spatial Transcriptomics. Nat Biotechnol 2022, 40, 1349–1359, doi:10.1038/s41587-022-01273-7.
8. 	Bergenstråhle, L.; He, B.; Bergenstråhle, J.; Abalo, X.; Mirzazadeh, R.; Thrane, K.; Ji, A.L.; Andersson, A.; Larsson, L.; Stak-enborg, N.; et al. Super-Resolved Spatial Transcriptomics by Deep Data Fusion. Nat Biotechnol 2022, 40, 476–479, doi:10.1038/s41587-021-01075-3.
9. 	Yan, L.; Sun, X. Benchmarking and Integration of Methods for Deconvoluting Spatial Transcriptomic Data. Bioinformatics 2023, 39, btac805, doi:10.1093/bioinformatics/btac805.
10. Beygelzimer, A.; Kakade, S.; Langford, J. Cover Trees for Nearest Neighbor. In Proceedings of the Proceedings of the 23rd international conference on Machine learning  - ICML ’06; ACM Press: Pittsburgh, Pennsylvania, 2006; pp. 97–104.
11. Satija, R.; Farrell, J.A.; Gennert, D.; Schier, A.F.; Regev, A. Spatial Reconstruction of Single-Cell Gene Expression Data. Nat Biotechnol 2015, 33, 495–502, doi:10.1038/nbt.3192.
12. Yu, G.; Wang, L.-G.; Han, Y.; He, Q.-Y. clusterProfiler: An R Package for Comparing Biological Themes among Gene Clusters. OMICS 2012, 16, 284–287, doi:10.1089/omi.2011.0118.
