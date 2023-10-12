# ucell

[![codecov](https://codecov.io/gh/maximilian-heeg/UCell/graph/badge.svg?token=YH55RI83BQ)](https://codecov.io/gh/maximilian-heeg/UCell)
![CI](https://github.com/maximilian-heeg/UCell/actions/workflows/ci.yml/badge.svg)

calculation of scores for different gene signatures on single-cell expression data

## Installation

```bash
$ pip install git+https://github.com/maximilian-heeg/UCell.git
```

## Usage

```python
import ucell
import scanpy as sc
# Load sample data
adata = sc.datasets.pbmc3k_processed()

# Calculate and add signatures
signatures = {
    "T cell": ['CD3E', 'CD4', 'CD8A', 'CD19-'],
    "CD4 cell": ['CD3E', 'CD4', 'CD8A-', 'CD19-', 'IL7R'],
    "CD8 cell": ['CD3E', 'CD4-', 'CD8A+', 'CD19-', 'IL7R-', 'GZMB', 'NCAM1-', 'FCGR3A-'],
    "B cells": ['CD19', 'CD3E-', 'FAKE GENE', 'CD37', 'CD8A-', 'CD4-', 'ITGAX-']
}

ucell.add_scores(adata, signatures=signatures, maxRank=1000)

# Plot UMAPs with signatures
sc.pl.umap(
    adata,
    color=["UCell_" + k for k in signatures.keys()],
    ncols=2
)
```

See a more detailed example in the [documentation](https://ucell.heeg.io/).

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`ucell` was created by Maximilian Heeg. It is licensed under the terms of the MIT license.


