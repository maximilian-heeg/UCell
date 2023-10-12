from ucell import ucell
import scanpy as sc
import pytest


def test_basic():
    adata = sc.datasets.pbmc3k_processed()
    adata = adata[0:10,]

    signatures = {
        "CD8 cell": ['CD3E', 'CD4-', 'CD8A+', 'CD19-', 'IL7R-', 'GZMB', 'NCAM1-', 'FCGR3A-']
    }


    assert ucell.calc_scores(adata, signatures).values[4].item() == pytest.approx(0.29, 0.1)
