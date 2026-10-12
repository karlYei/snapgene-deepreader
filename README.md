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

## Node Annotations, Step Primers And Visibility

Version 0.4.0 exposes the annotations actually saved for each historical node:

- `features`, `primers`, `notes`, `primer_settings`: node DNA annotations. A matched
  snapshot takes precedence over tree annotations; the current document takes
  precedence for the validated current node. `annotation_sources` identifies
  the source of each field. Missing records give empty values, not inherited
  annotations from a different molecule.
- `tree_features`, `tree_primers`: the tree's annotations, kept separately from
  the full snapshot. These may be a reduced selection relevant to a later step.
- `used_primers`: each step's saved `Oligo` name, sequence, description and raw
  attributes. This is distinct from `primers` attached to its output DNA. Binding
  sites are linked to direct input templates by matching saved primer sequences
  (case-insensitive, preferring a matching name). Operation-specific tree sites
  take precedence; snapshot sites are a fallback and may include extra matches.
  Sites retain `template_key`, `source`, `evidence`, `bound_strand`,
  `annealed_bases`, `melting_temperature`, `simplified` and `components`.
- `parameters`, `history_colors`, `inputs`, `consumer_keys`: saved step settings,
  colors and both directions of the dependency graph. Input summaries are paired
  with inputs in saved order only when their counts agree. Other XML stays in
  `properties`; snapshot XML and binary packet summaries remain accessible via
  `history["snapshots"][node["snapshot_index"]]` when that index is present.
- `annotation_errors`: local annotation problems do not discard a node's sequence.

Binding-site `location` and component `hybridized_range` strings are retained.
Recognized ranges also have a `region` with **zero-based inclusive** `start` and
`end`, `coordinate_system`, and `wraps_origin`. Feature `start`/`end` remain
zero-based half-open; raw feature segment XML ranges remain one-based inclusive.
`InputSummary` values are retained verbatim: their meanings depend on the saved
manipulation. A step's oligos are recorded reagents, not necessarily laboratory
stock primers (some operations save virtual oligos such as `Fused.FOR`).

`binding_sites_status="unavailable"` means the file did not provide a matched
saved region. The parser does not invent primer positions, a forward/reverse
role, or a primer pair from names. Missing sequence snapshots also remain
unavailable; saved metadata can still be read.

```python
for node in data["history"]["nodes"]:
    print(node["name"], node["operation"], node["is_visible"])
    print(node.get("sequence"), node["features"], node["primers"])
    for primer in node["used_primers"]:
        print(primer["name"], primer["sequence"], primer["binding_sites"])
```

`hide_history` parses the node's saved `hideHistory` flag. It hides that node's
**earlier input subtree**, not the boundary node itself or sibling branches.
`is_visible` is the effective display state; `hidden_by` lists the downstream
boundary keys causing the hiding. Hidden nodes are still fully parsed and
returned. Trees expose `visible_node_keys`, `hidden_node_keys`, and
`unknown_visibility_node_keys`. Unknown flag values produce an unknown (`None`)
state for affected upstream nodes rather than guessing. Feature `visible` is
separate from history-node visibility.

Multi-attribute feature qualifiers (including `bool`/`text`/`predef` inference
values) no longer fail fixed-length unpacking. Complex values keep all their
fields, repeated values are not collapsed, and original XML mappings are kept in
`qualifier_values`; feature and primer `attributes` preserve extra metadata.
The established simple qualifier representations remain supported. Structured
qualifier values are serialized as JSON text for SeqRecord/GenBank export;
GenBank quote escaping and single-segment coordinate conversion are handled.
Historical annotations do not imply a complete historical SeqRecord: the
existing export APIs still export the current document.

## Compatibility

The following functions are provided with the same call style as
`snapgene-reader`:

- `snapgene_file_to_dict(filepath=None, fileobject=None)`
- `snapgene_file_to_seqrecord(filepath=None, fileobject=None)`
- `snapgene_file_to_gbk(read_file_object, write_file_object)`

Unknown or unsupported rich blocks are not fatal. Parse failures are reported in
the relevant returned object with a `parse_error` field.
