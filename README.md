# snapgene-deepreader

[English](README.md) | [中文](README.zh-CN.md)

`snapgene-deepreader` is an enhanced SnapGene `.dna` parser based on the
MIT-licensed SnapGeneReader project. It keeps the original high-level API while
adding extraction of primer blocks, SnapGene alignable sequence blocks, and
SnapGene-saved comparison alignment data and construction history.

## Origin and development

This project is derived from the open-source
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)
project by Edinburgh Genome Foundry, which is released under the MIT License.

All changes that make this a separate package named `snapgene-deepreader` were
developed entirely by an AI coding agent. The added work includes package
renaming, primer block parsing, alignable sequence metadata parsing, zlib-backed
alignable sequence extraction, BAM-backed comparison parsing, history parsing
and sequence recovery, tests, and project
documentation.

## Install for local development

```bash
python -m pip install -e .
```

## Usage

```python
from snapgene_deepreader import snapgene_file_to_dict

data = snapgene_file_to_dict("example.dna")
print(data["seq"])
print(data["primers"])
print(data["alignable_sequences"])
print(data["comparison_alignments"])
print(data["history"])
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
- `comparison_alignments`: parsed from SnapGene block `27` when it contains a
  BGZF/GZIP-compressed BAM alignment. Each BAM record includes CIGAR, sequence,
  auxiliary tags, linked alignable sequence metadata, and a `variants` list
  reconstructed against the SnapGene reference sequence.
- `raw_blocks_summary`: block type, size, and byte offset for debugging.

## Saved construction history

Version 0.3.0 adds `data["history"]`, read from the saved HistoryTree (block `7`)
and recovery data (block `11`). It always has these fields, including when no
history is saved:

- `trees`: tree attributes, `root_keys`, original XML, and source block locations.
- `nodes`: one entry per node occurrence, with `key`, original `id`, `name`,
  `operation`, `sequence_length`, `topology`, `strandedness`, `input_keys`,
  `input_summaries`, and unmodified XML `attributes` and `properties`.
- `snapshots`: node IDs, Undo/Redo actions, decoded external fragments,
  saved XML properties (including features when present), and packet summaries.
- `parse_errors`: failed block/packet locations and error messages. Local errors
  are also recorded on snapshots and packets as `parse_error`.

`input_keys` link a product to its inputs. Multiple inputs and repeated IDs are
preserved; use occurrence `key` values for traversal, not IDs or list positions.
Nodes are listed in tree traversal order, not chronological order. Dates, authors,
and other attributes are retained if saved, but never inferred.

```python
history = data["history"]
nodes = {node["key"]: node for node in history["nodes"]}
for tree in history["trees"]:
    for root_key in tree["root_keys"]:
        product = nodes[root_key]
        inputs = [nodes[key] for key in product["input_keys"]]
        print(product["name"], product["operation"], [n["name"] for n in inputs])

for node in history["nodes"]:
    print(node["id"], node["recovery_status"], node.get("recovery_reason"))
    if "sequence" in node:
        print(node["sequence"])
```

The `sequence` field is present only when available. `recovery_status` is
`current` (the document sequence), `recovered`, or `unavailable`.
`sequence_source` identifies `document`, `snapshot`, or `undo`; unavailable
nodes have `recovery_reason`. Conflicting validation results use
`recovery_warning`. `snapshot_index`, when present, links to `snapshots`.

Recovery supports observed two-bit compact DNA run encoding `1`, including
lowercase ranges and external fragments, plus saved INSERT/REMOVE actions
(including replacement as remove then insert) and whole-sequence FLIP.
Recovered lengths are checked against tree metadata; saved Redo actions are
used to verify Undo results when present. History action coordinates are
**zero-based with inclusive removal endpoints** and are retained verbatim.

This covers verified formats, not all SnapGene versions. Other compact encodings,
circular origin-spanning removals, ambiguous IDs, and unsupported actions remain
unavailable with diagnostic information. A multi-input construction is retained
as a graph; it is not automatically reversed into its inputs. XML properties are
preserved separately from the recovered sequence, not promised as a complete
historical SeqRecord. There are limits on decompression, XML nesting, and total
recovered sequence size. No new runtime dependency is required.

The original SeqRecord and GenBank APIs still export the current document.
History represents construction data actually saved in the file, not a complete
log of every edit or save.

## Compatibility

The following functions are provided with the same call style as
`snapgene-reader`:

- `snapgene_file_to_dict(filepath=None, fileobject=None)`
- `snapgene_file_to_seqrecord(filepath=None, fileobject=None)`
- `snapgene_file_to_gbk(read_file_object, write_file_object)`

Unknown or unsupported rich blocks are not fatal. Parse failures are reported in
the relevant returned object with a `parse_error` field.
