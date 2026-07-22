"""
Load Reactome gene-pathway mapping and compute pathway-level score components.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REACTOME_COLS = ["gene_id", "reactome_id", "url", "pathway_name", "evidence", "species"]


def load_reactome_map(path: Path, species: str = "Homo sapiens") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=REACTOME_COLS)
    return df[df["species"] == species][["gene_id", "reactome_id", "pathway_name"]].copy()


def pathway_score_components(
    score_matrix: pd.DataFrame,
    v2g: pd.DataFrame,
    g2p: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute pathway-level score components.

    pathway_score_i,k = Σ genotype_ij × beta_j  for variants j in pathway k

    Parameters
    ----------
    score_matrix : individual × variant weight-applied matrix
    v2g : variant-to-gene mapping  [variant_id, gene_id]
    g2p : gene-to-pathway mapping  [gene_id, reactome_id, pathway_name]

    Returns
    -------
    pd.DataFrame: individual × pathway component scores (long format)
    """
    merged = v2g.merge(g2p, on="gene_id", how="inner")
    pathways = merged["reactome_id"].unique()

    records = []
    for pathway_id in pathways:
        variants_in_path = merged.loc[
            merged["reactome_id"] == pathway_id, "variant_id"
        ].unique()
        cols_present = [v for v in variants_in_path if v in score_matrix.columns]
        if not cols_present:
            continue
        pathway_name = merged.loc[
            merged["reactome_id"] == pathway_id, "pathway_name"
        ].iloc[0]
        component = score_matrix[cols_present].sum(axis=1)
        for sample_id, val in component.items():
            records.append(
                {
                    "sample": sample_id,
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                    "pathway_score": val,
                }
            )
    return pd.DataFrame(records)
