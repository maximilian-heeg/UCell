import pandas as pd
import anndata as ad
from typing import List, Dict


def calc_scores(
    adata: ad.AnnData,
    signatures: Dict[str, List[str]],
    prefix="UCell_",
    maxRank=100,
    seed=None,
) -> pd.DataFrame:
    """
    Calculate scores for different gene signatures on single-cell expression data.

    Args:
        adata (ad.AnnData): An AnnData object representing single-cell expression data.
        signatures (Dict[str, List[str]): A dictionary of gene signatures where keys are signature names and values are lists of gene names. Genes can have a direction as indicated by the trailing minus or plus sign.
        prefix (str, optional): Prefix for the column names in the output DataFrame (default is "UCell\_").
        maxRank (int, optional): Maximum rank value (default is 100).
        seed (int, optional): Seed for random number generation (default is None).

    Returns:
        pd.DataFrame: A DataFrame with calculated scores for each gene signature.
    """
    from scipy.sparse import csr_matrix

    if adata.raw:
        x = adata.raw.X
        columns = adata.raw.var_names
    else:
        x = adata.X
        columns = adata.var_names

    if isinstance(x, csr_matrix):
        x = x.todense()
    m = pd.DataFrame(x, columns=columns)

    signatures = _check_signatures(signatures=signatures, indices=m.columns)

    m = create_rankings(pd.DataFrame(m), seed=seed)

    scores: Dict[str, pd.Series] = {}
    for k, v in signatures.items():
        scores[prefix + k] = _calc_score(m, v, maxRank=maxRank)

    scores: pd.DataFrame = pd.DataFrame(scores)
    scores[scores < 0] = 0
    scores.index = adata.obs.index
    return scores


def add_scores(
    adata: ad.AnnData,
    signatures: Dict[str, List[str]],
    prefix="UCell_",
    maxRank=100,
    seed=None,
):
    """
    Calculate scores for different gene signatures on single-cell expression data and add it to AnnData obs.

    Args:
        adata (ad.AnnData): An AnnData object representing single-cell expression data.
        signatures (Dict[str, List[str]): A dictionary of gene signatures where keys are signature names and values are lists of gene names. Genes can have a direction as indicated by the trailing minus or plus sign.
        prefix (str, optional): Prefix for the column names in the output DataFrame (default is "UCell\_").
        maxRank (int, optional): Maximum rank value (default is 100).
        seed (int, optional): Seed for random number generation (default is None).
    """
    scores = calc_scores(
        adata=adata, signatures=signatures, prefix=prefix, maxRank=maxRank, seed=seed
    )
    _add_to_adata(adata=adata, scores=scores)


def _check_signatures(
    signatures: Dict[str, List[str]], indices: List[str]
) -> Dict[str, List[str]]:
    """
    Filter a dictionary of gene signatures by checking the presence of genes in a list of indices.

    Args:
        signatures (Dict[str, List[str]]): A dictionary where keys are signature names and values are lists of genes.
        indices (List[str]): A list of genes to check against the signature genes.

    Returns:
        Dict[str, List[str]]: A filtered dictionary with the same keys as the input 'signatures',
        where the values are lists of genes that were found in 'indices'.

    Warnings:
        If some genes in a signature are not found in 'indices', a warning message is issued indicating
        which genes are missing and that they will be removed from the list.
    """
    from itertools import compress
    import warnings

    filtered_list: Dict[str, List[str]] = {}

    for k, v in signatures.items():
        stripped = [gene.strip(r"\+|\-") for gene in v]
        exist = [gene in indices for gene in stripped]
        exist = list(compress(v, exist))
        missing = list(set(v) - set(exist))

        if len(missing):
            warnings.warn(
                f"Some genes were not found in signature: {k}: {missing}. Missing genes will be removed from the list.",
                stacklevel=2,
            )

        filtered_list[k] = exist

    return filtered_list


def _calc_score(m: pd.DataFrame, signature: List[str], maxRank=100) -> pd.Series:
    """
    Internal helper function that splits the gene list is up and downregulated genes.
    Then a score for each signature is calculated and the difference returned.

    Args:
        m (pd.DataFrame): A DataFrame containing rankings of genes in single-cell expression data.
        signature (List[str]): A list of gene names representing a gene signature.
        maxRank (int, optional): Maximum rank value (default is 100).

    Returns:
        pd.Series: A Series containing calculated scores for the gene signature.
    """
    sig_neg = [m.strip(r"\+|\-") for m in signature if m.endswith("-")]
    sig_pos = [m.strip(r"\+|\-") for m in signature if not m.endswith("-")]

    return _u_stat(ranks=m[sig_pos], maxRank=maxRank) - _u_stat(
        ranks=m[sig_neg], maxRank=maxRank
    )


def create_rankings(ex_mtx: pd.DataFrame, seed=None) -> pd.DataFrame:
    """
    Create a rankings dataframe from a single cell expression profile dataframe.

    Args:
        ex_mtx (pd.DataFrame): The expression profile matrix. Rows correspond to different cells, columns to different genes (n_cells x n_genes).
        seed (int, optional): Seed for random number generation (default is None).

    Returns:
        pd.DataFrame: A DataFrame with gene rankings for each cell.
    """

    return (
        ex_mtx.sample(frac=1.0, replace=False, axis=1, random_state=seed)
        .rank(axis=1, ascending=False, method="first", na_option="bottom")
        .astype("uint32")
    )


def _u_stat(ranks: pd.DataFrame, maxRank=100) -> pd.Series:
    """
    Calculate the U-statistic for a set of gene rankings.

    Args:
        ranks: A DataFrame with gene rankings.
        maxRank (int, optional): Maximum rank value (default is 100).

    Returns:
        pd.Series: The calculated U-statistic for the input rankings.
    """
    ranks = ranks.copy()
    ranks[ranks > maxRank] = maxRank + 1
    rank_sum = ranks.sum(axis=1)
    len_sig = len(ranks.columns)
    if len_sig == 0:
        return 0
    u_value = rank_sum - (len_sig * (len_sig + 1)) / 2
    auc = 1 - u_value / (len_sig * maxRank)
    return auc


def _add_to_adata(adata: ad.AnnData, scores: pd.DataFrame):
    scores.index = adata.obs.index
    for signature in scores.columns:
        adata.obs[[signature]] = scores[[signature]]
