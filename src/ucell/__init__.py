"""
The UCell package is a designed to facilitate
the calculation of scores for different gene signatures on single-cell 
expression data. It provides a set of functions and utilities for working
with single-cell gene expression data and gene signatures
"""
__version__ = "0.1"

from .ucell import (
    calc_scores,
    add_scores,
    create_rankings
)
