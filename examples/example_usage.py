import pandas as pd
from genomic_helper.background import compute_background_peaks

# Example mock dataframe
df = pd.DataFrame({
    "prop_shared_cells": [0.1, 0.2, 0.05, 0.3],
    "GC": [0.4, 0.5, 0.45, 0.6],
    "pct_dropout": [0.7, 0.8, 0.6, 0.9],
    "Median Rank": [5, 10, 2, 7],
    "IQR Rank": [1, 3, 2, 4],
    "overlapping": ["geneA", "geneB", "geneA", "geneC"],
    "Chromosome": ["chr1", "chr2", "chr1", "chr3"],
    "Start": [100, 200, 300, 400],
    "End": [150, 250, 350, 450]
})

result = compute_background_peaks(
    df=df,
    feature_cols=["prop_shared_cells", "GC", "pct_dropout", "Median Rank", "IQR Rank"],
    group_col="overlapping",
    chrom_col="Chromosome",
    k=2
)

print(result)
