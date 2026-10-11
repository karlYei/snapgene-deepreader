"""Enhanced SnapGene .dna reader."""

from .snapgene_reader import (
    snapgene_file_to_dict,
    snapgene_file_to_gbk,
    snapgene_file_to_seqrecord,
)

__version__ = "0.3.0"

__all__ = [
    "snapgene_file_to_dict",
    "snapgene_file_to_gbk",
    "snapgene_file_to_seqrecord",
]
