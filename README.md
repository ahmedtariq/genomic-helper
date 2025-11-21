# genomic-helper

A lightweight Python package for genomic background peak selection and feature-based similarity matching.  
Designed to be fully generic — works for any numeric feature set, not only genomics.

## Features
✔ Generic: accept any numeric features  
✔ Normalization included  
✔ Efficient computation of background peaks  
✔ Excludes same chromosome + same group  
✔ Uses sklearn pairwise distances  
✔ pyproject.toml modern packaging  

---

## Installation

Clone the repo:
```bash
git clone https://github.com/yourusername/genomic-helper.git
cd genomic-helper
pip install .
```
---

## Example

```python
from genomic_helper.background import compute_background_peaks

df_out = compute_background_peaks(
    df=result_df,
    feature_cols=["prop_shared_cells", "GC", "pct_dropout", "Median Rank", "IQR Rank"],
    group_col="overlapping",
    chrom_col="Chromosome",
    k=5
)
```