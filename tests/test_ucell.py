from ucell import ucell
import scanpy as sc
import pytest
import pandas as pd


# functions to create fake data
def rank_data() -> pd.DataFrame:
    ranks = pd.DataFrame(
        {"A": [1, 2, 95], "B": [4, 70, 74], "C": [7, 75, 110]},
        index=["cell1", "cell2", "cell3"],
    )
    return ranks


def create_anndata() -> sc.AnnData:
    adata = sc.AnnData(X=rank_data())
    return adata


# tests
def test_u_score():
    ranks = rank_data()
    scores = ucell.__u_stat(ranks)
    assert list(scores) == [0.98, 0.53, 0.12]


def test_u_score_maxRank():
    ranks = rank_data()
    scores = ucell.__u_stat(ranks, maxRank=50)
    assert list(scores) == pytest.approx([0.96, 0.346, 0.02], 0.01)


def test_u_score_maxRank_2():
    ranks = rank_data()
    scores = ucell.__u_stat(ranks, maxRank=500)
    assert list(scores) == pytest.approx([0.996, 0.906, 0.818], 0.01)


def test_calc_score():
    ranks = rank_data()
    ranks = ucell.create_rankings(ranks)
    signatures = {"Sig1": ["A", "C"], "Sig2": ["A-", "C"]}
    scores = ucell.__calc_score(ranks, signatures["Sig1"])
    assert list(scores) == pytest.approx([0.995, 0.995, 1.0], 0.01)
    scores = ucell.__calc_score(ranks, signatures["Sig2"])
    assert list(scores) == pytest.approx([0.02, 0.02, 0.01], 0.01)


def test_create_rankings():
    ranks = rank_data()
    ranks = ucell.create_rankings(ranks)
    # order of the columns can change
    ranks = ranks[["A", "B", "C"]]
    expected = pd.DataFrame(
        {"A": [3, 3, 2], "B": [2, 2, 3], "C": [1, 1, 1]},
        index=["cell1", "cell2", "cell3"],
    )
    assert ranks.equals(expected.astype("uint32"))


def test_warning():
    adata = create_anndata()
    # Make sure this triggers a warning as 'D' does not exist
    with pytest.warns(UserWarning):
        scores = ucell.calc_scores(adata, {"Sig1": ["A", "B", "D"]})
    assert list(scores.UCell_Sig1) == pytest.approx([0.99, 0.99, 0.99], 0.01)


def test_calc_scores():
    adata = create_anndata()
    scores = ucell.calc_scores(adata, {"Sig1": ["A", "B"]}, maxRank=2)
    assert list(scores.UCell_Sig1) == pytest.approx([0.5, 0.5, 0.5], 0.01)


def test_check_signatures():
    signatures = {"Sig1": ["A", "B", "D"], "Sig2": ["C", "B", "D", "E"]}
    with pytest.warns(UserWarning):
        signatures = ucell.__check_signatures(signatures, ["A", "B", "E"])
    expected = {"Sig1": ["A", "B"], "Sig2": ["B", "E"]}
    assert signatures == expected


def test_add_to_anndata():
    adata = create_anndata()
    scores = ucell.calc_scores(adata, {"Sig1": ["A", "B"]}, maxRank=2)
    ucell.__add_to_adata(adata, scores)
    assert adata.obs.columns.item() == "UCell_Sig1"


def test_add_scores():
    adata = create_anndata()
    ucell.add_scores(adata, {"Sig1": ["A", "B"]}, maxRank=2)
    assert adata.obs.columns.item() == "UCell_Sig1"
    assert list(adata.obs["UCell_Sig1"]) == [0.5, 0.5, 0.5]


def test_raw_and_sparse():
    from scipy.sparse import csr_matrix

    adata = create_anndata()
    adata.X = csr_matrix(adata.X)
    adata.raw = adata
    adata.X = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    scores = ucell.calc_scores(adata, {"Sig1": ["A", "B"]}, maxRank=2)
    assert list(scores.UCell_Sig1) == pytest.approx([0.5, 0.5, 0.5], 0.01)
