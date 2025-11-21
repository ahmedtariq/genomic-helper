import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from tqdm import tqdm


def compute_background_peaks(
    df,
    feature_cols,
    group_col="overlapping",
    chrom_col="Chromosome",
    k=5,
    distance_metric="euclidean",
):
    """
    Compute background peaks for each group using nearest-neighbor similarity
    in the given feature space, excluding peaks from the same chromosome
    and same group.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing features and annotations.
    feature_cols : list of str
        Numeric feature columns used for similarity.
    group_col : str
        Column representing groups (e.g., genes).
    chrom_col : str
        Column representing chromosomes.
    k : int
        Number of nearest neighbors to include.
    distance_metric : str
        Metric used by sklearn.pairwise_distances.

    Returns
    -------
    pd.DataFrame
        Combined original dataframe with appended background rows.
    """

    df = df.copy()

    # Normalize features
    norm_df = (df[feature_cols] - df[feature_cols].min()) / (
        df[feature_cols].max() - df[feature_cols].min()
    )

    background_rows = []

    for group in tqdm(df[group_col].unique(), desc="Computing background"):

        chrom = df.loc[df[group_col] == group, chrom_col].iloc[0]

        ref = norm_df[df[group_col] == group]

        candidates_mask = (df[group_col] != group) & (df[chrom_col] != chrom)
        candidates = norm_df[candidates_mask]

        candidates = candidates.loc[~candidates.index.duplicated(keep="first")]

        if candidates.empty:
            continue

        distances = pairwise_distances(ref, candidates, metric=distance_metric)

        mean_dist = pd.Series(distances.mean(axis=0), index=candidates.index)

        nearest = mean_dist.sort_values().iloc[:k].index

        bk = df.loc[nearest].copy()
        bk[group_col] = group
        bk["distance"] = np.nan

        background_rows.append(bk)

    background_df = pd.concat(background_rows, ignore_index=False)
    return pd.concat([df, background_df], ignore_index=False)