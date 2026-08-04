# snapgene-deepreader

[English](README.md) | [中文](README.zh-CN.md)

`snapgene-deepreader` is an enhanced SnapGene `.dna` parser based on the
MIT-licensed SnapGeneReader project. It keeps the original high-level API while
adding extraction of primer blocks and SnapGene alignable sequence blocks.

## Origin and development

This project is derived from the open-source
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)
project by Edinburgh Genome Foundry, which is released under the MIT License.

All changes that make this a separate package named `snapgene-deepreader` were
developed entirely by an AI coding agent. The added work includes package
renaming, primer block parsing, alignable sequence metadata parsing, zlib-backed
alignable sequence extraction, tests, and project documentation.

## Install for local development

```bash
/home/karl/miniconda3/bin/python -m pip install -e .
```

## Usage

```python
from snapgene_deepreader import snapgene_file_to_dict

data = snapgene_file_to_dict("example.dna")
print(data["seq"])
print(data["primers"])
print(data["alignable_sequences"])
```

The returned dictionary keeps the original keys from `snapgene-reader`:

- `dna`
- `seq`
- `features`
- `notes`
- `isDNA`
- `exportVersion`
- `importVersion`

It also adds:

- `primers`: parsed from SnapGene block `5`.
- `alignable_sequences`: parsed from blocks `16` and `17`.
- `raw_blocks_summary`: block type, size, and byte offset for debugging.

## Compatibility

The following functions are provided with the same call style as
`snapgene-reader`:

- `snapgene_file_to_dict(filepath=None, fileobject=None)`
- `snapgene_file_to_seqrecord(filepath=None, fileobject=None)`
- `snapgene_file_to_gbk(read_file_object, write_file_object)`

Unknown or unsupported rich blocks are not fatal. Parse failures are reported in
the relevant returned object with a `parse_error` field.
