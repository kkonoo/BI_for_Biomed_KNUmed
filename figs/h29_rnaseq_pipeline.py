from gallery import *
import os

# Bulk RNA-seq pipeline in the same look as the W3 pipelines: reuse the helpers in g13_pipelines.py
# (everything before its first figure). Steps follow the lab pipeline: FastQC → TrimGalore → salmon →
# tximport → DESeq2 → clusterProfiler / WGCNA.
exec(open(os.path.join(os.path.dirname(__file__), "g13_pipelines.py")).read().split("# a) ATAC-seq")[0])

three_band("rnaseq_pipeline.png", "Pipeline · bulk RNA-seq (lab example)",
    "Bulk RNA-seq: from reads to differentially expressed genes",
    pre=[("Sequencer output", "reads per sample", "FASTQ"),
         ("Trim reads", "adapters, low-quality bases\nTrimGalore", "FASTQ"),
         ("Quantify transcripts", "pseudo-alignment · salmon\n(or STAR + featureCounts)", ".sf"),
         ("Sum to genes", "tximport", "counts"),
         ("Normalise, transform", "size factors (DESeq2)\nVST for plots only", None)],
    qc=[("Raw reads", "base quality · adapters · GC\nFastQC / MultiQC"),
        ("Read alignment", "mapping rate · rRNA %\ngene-body coverage"),
        ("Count matrix", "library size · genes detected"),
        ("Samples (EDA)", "PCA on VST · batch · outliers")],
    test=[("DESeq2 model", "negative binomial GLM\n~ batch + condition\ndispersion shrinkage", None),
          ("Wald test / LRT", "log₂ fold change · padj", "table")],
    res="Genes that differ between groups",
    res_sub="→ ORA / GSEA (clusterProfiler) · WGCNA · deconvolution",
    sub_pre="reads → a genes × samples count matrix", sub_qc="is every sample usable?",
    sub_test="which genes differ between groups?")
