from gallery import *
import os

# scRNA-seq pipeline in the same look as the W3/W4 pipelines (helpers from g13_pipelines.py).
# Steps follow the lab pipeline: Cell Ranger → CellBender → cell QC → scDblFinder → HVG/PCA →
# scVI/Harmony + clustering → annotation → markers / pseudobulk DE → downstream.
exec(open(os.path.join(os.path.dirname(__file__), "g13_pipelines.py")).read().split("# a) ATAC-seq")[0])

three_band("sc_pipeline.png", "Pipeline · scRNA-seq (lab example)",
    "scRNA-seq: from reads to cell types and their changes",
    pre=[("Align and count UMIs", "Cell Ranger / BD Rhapsody\n(+ demultiplexing)", ".h5"),
         ("Call cells, remove ambient RNA", "CellBender (10x data)", ".h5"),
         ("Primary cell QC", "UMI · genes · mito % · scanpy / Seurat", ".h5ad"),
         ("Remove doublets", "scDblFinder", None),
         ("Normalise, HVG, PCA", "cells × PCs", None),
         ("Integrate, cluster", "scVI / Harmony · Leiden", None),
         ("Annotate", "marker genes · scANVI\n(+ secondary QC)", "labels")],
    qc=[("Reads and alignment", "valid barcodes · saturation\nreads mapped to the transcriptome"),
        ("Knee plot", "cells vs. empty droplets"),
        ("Per-cell metrics", "UMI · genes · mito %"),
        ("Ambient signal", "one cell type's markers in every cluster"),
        ("Doublet score", "cells with two types' markers"),
        ("Integration check", "colour by batch and by cell type")],
    test=[("Marker genes", "Wilcoxon, within one dataset", "table"),
          ("Pseudobulk DE", "donor × cell type sums\nDESeq2 / edgeR", "DEGs"),
          ("Differential abundance", "cell-type proportions\nscCODA · miloR", None)],
    res="Cell types, and how they change between conditions",
    res_sub="→ trajectory (CellRank) · cell–cell communication (CellPhoneDB) · regulons (pySCENIC)",
    sub_pre="reads → annotated cells", sub_qc="what to check at each step",
    sub_test="n = donors, not cells")
