"""
Pandera schema definitions for all major pipeline tables.
"""

import pandera as pa
from pandera.typing import Series


class HarmonizedScoreFileSchema(pa.DataFrameModel):
    chr: Series[str]
    pos: Series[int]
    effect_allele: Series[str]
    other_allele: Series[str]
    effect_weight: Series[float]
    variant_id: Series[str]
    pgs_id: Series[str]
    trait: Series[str]

    class Config:
        coerce = True
        strict = False


class SampleMetadataSchema(pa.DataFrameModel):
    sample: Series[str]
    population: Series[str]
    super_population: Series[str]
    gender: Series[str]

    class Config:
        coerce = True
        strict = False


class AncestryPCsSchema(pa.DataFrameModel):
    sample: Series[str]
    PC1: Series[float]
    PC2: Series[float]
    PC3: Series[float]

    class Config:
        coerce = True
        strict = False


class PopulationShiftSchema(pa.DataFrameModel):
    pgs_id: Series[str]
    population: Series[str]
    super_population: Series[str]
    mean_score: Series[float]
    sd_score: Series[float]
    standardized_shift: Series[float]
    ci_lower: Series[float]
    ci_upper: Series[float]

    class Config:
        coerce = True
        strict = False


class PathwayComponentSchema(pa.DataFrameModel):
    sample: Series[str]
    pgs_id: Series[str]
    pathway_id: Series[str]
    pathway_name: Series[str]
    pathway_score: Series[float]

    class Config:
        coerce = True
        strict = False
