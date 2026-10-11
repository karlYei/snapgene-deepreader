import json
import lzma
import struct
from io import BytesIO
from pathlib import Path

import pytest

from snapgene_deepreader import snapgene_file_to_dict, snapgene_file_to_seqrecord
from snapgene_deepreader import history as history_parser


def block(kind, payload):
    return bytes([kind]) + struct.pack(">I", len(payload)) + payload


def document(sequence, *blocks):
    header = block(9, b"SnapGene" + struct.pack(">HHH", 1, 1, 1))
    return BytesIO(
        header
        + b"".join(blocks)
        + block(0, b"\x03" + sequence.encode())
        + block(6, b"<Notes><Type>Natural</Type></Notes>")
    )


def tree(root_id, sequence_length, inputs, operation="insert"):
    return block(
        7,
        lzma.compress(
            (
                '<HistoryTree extra="kept"><Node ID="%s" name="product" type="DNA" '
                'seqLen="%s" circular="1" operation="%s">%s</Node></HistoryTree>'
                % (root_id, sequence_length, operation, inputs)
            ).encode()
        ),
    )


def node(node_id, length, children=""):
    return (
        '<Node ID="%s" name="source" type="DNA" seqLen="%s" circular="1">%s</Node>'
        % (node_id, length, children)
    )


def manipulation(undo, redo=""):
    return block(
        29,
        lzma.compress(
            (
                "<Manipulation><Undo>%s</Undo><Redo>%s</Redo></Manipulation>"
                % (undo, redo)
            ).encode()
        ),
    )


def snapshot(node_id, *packets):
    return block(11, struct.pack(">I", node_id) + b"".join(packets))


def compact(length, packed, lowercase=()):
    # Fixed packed bytes in tests verify G/A/T/C ordering and right-aligned tails.
    return (
        struct.pack(">IBII", length, 3, 1, len(lowercase))
        + struct.pack(">BI", 1, length)
        + packed
        + b"".join(struct.pack(">II", start, end) for start, end in lowercase)
    )


def test_reads_tree_unknown_properties_and_recovers_multiple_steps():
    properties = '<InputSummary manipulation="insertAt" val1="4"/><Unknown flag="x">hello</Unknown>'
    older = node(0, 5)
    source = node(1, 4, older + '<Features><Feature name="old feature"/></Features>')
    stream = document(
        "ACGTTT",
        snapshot(
            0,
            manipulation(
                '<Action type="INSERT" position="1"><Residues type="GENERIC" residues="C"/></Action>',
                '<Action type="REMOVE" range="1-1"/>',
            ),
        ),
        snapshot(
            1,
            manipulation(
                '<Action type="REMOVE" range="4-5"/>',
                '<Action type="INSERT" position="4"><Residues type="GENERIC" residues="TT"/></Action>',
            ),
            block(
                30,
                lzma.compress(
                    block(
                        10, b'<Features><Feature name="snapshot feature"/></Features>'
                    )
                ),
            ),
        ),
        tree(2, 6, properties + source),
    )
    data = snapgene_file_to_dict(fileobject=stream)
    history = data["history"]
    assert not stream.closed
    assert history["parse_errors"] == []
    assert [n["sequence"] for n in history["nodes"]] == ["ACGTTT", "ACGT", "ACCGT"]
    current, old, oldest = history["nodes"]
    assert current["input_keys"] == [old["key"]]
    assert old["input_keys"] == [oldest["key"]]
    assert current["input_summaries"][0]["attributes"]["val1"] == "4"
    assert current["properties"][1]["text"] == "hello"
    assert history["trees"][0]["attributes"] == {"extra": "kept"}
    assert "HistoryTree" in history["trees"][0]["raw_xml"]
    assert history["snapshots"][1]["properties"][0]["xml"]["tag"] == "Features"
    assert current["recovery_status"] == "current"
    assert oldest["sequence_source"] == "undo"
    json.dumps(history)
    stream.seek(0)
    record = snapgene_file_to_seqrecord(fileobject=stream)
    assert str(record.seq) == "ACGTTT"
    assert record.features == []


@pytest.mark.parametrize(
    "sequence,packed",
    [
        ("G", b"\x00"),
        ("A", b"\x01"),
        ("T", b"\x02"),
        ("C", b"\x03"),
        ("GATC", b"\x1b"),
        ("GATCA", b"\x1b\x01"),
        ("GATCAC", b"\x1b\x07"),
        ("GATCACT", b"\x1b\x1e"),
    ],
)
def test_recovers_compact_sources_with_partial_groups(sequence, packed):
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC",
            tree(2, 4, node(0, len(sequence)) + node(1, 1), "clone"),
            snapshot(0, block(1, compact(len(sequence), packed))),
        )
    )["history"]
    assert history["nodes"][1]["sequence"] == sequence
    assert history["nodes"][1]["sequence_source"] == "snapshot"
    assert history["nodes"][2]["recovery_status"] == "unavailable"
    assert "sequence" not in history["nodes"][2]


def test_recovers_external_fragment_and_lowercase_ranges():
    external = struct.pack(">I", 7) + block(1, compact(2, b"\x07", [(0, 1)]))
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC",
            tree(1, 4, node(0, 6), "remove"),
            snapshot(
                0,
                block(31, external),
                manipulation(
                    '<Action type="INSERT" position="4"><Residues type="EXTERNAL" ID="7" length="2"/></Action>',
                    '<Action type="REMOVE" range="4-5"/>',
                ),
            ),
        )
    )["history"]
    assert history["nodes"][1]["sequence"] == "GATCac"
    assert history["parse_errors"] == []


def test_restores_reverse_complement_when_only_flip_redo_saved():
    history = snapgene_file_to_dict(
        fileobject=document(
            "YACGT",
            tree(1, 5, node(0, 5), "flip"),
            snapshot(0, manipulation("", '<Action type="FLIP"/>')),
        )
    )["history"]
    assert history["nodes"][1]["sequence"] == "ACGTR"


def test_duplicate_node_ids_keep_distinct_occurrences_without_guessing():
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC",
            tree(1, 4, node(0, 4) + node(0, 4), "clone"),
            snapshot(0, block(1, compact(4, b"\x1b"))),
        )
    )["history"]
    assert len({n["key"] for n in history["nodes"]}) == 3
    for n in history["nodes"][1:]:
        assert "sequence" not in n
        assert "ambiguous" in n["recovery_reason"]


def test_duplicate_snapshots_are_not_selected_arbitrarily():
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC",
            tree(1, 4, node(0, 4)),
            snapshot(0, block(1, compact(4, b"\x1b"))),
            snapshot(0, block(1, compact(4, b"\xff"))),
        )
    )["history"]
    assert "ambiguous" in history["nodes"][1]["recovery_reason"]


@pytest.mark.parametrize(
    "payload",
    [
        b"\xfd7zXZ\0broken",
        b"<HistoryTree>",
        b"<Other/>",
        b'<!DOCTYPE HistoryTree [<!ENTITY x "value">]><HistoryTree>&x;</HistoryTree>',
        b"<HistoryTree>" + b"<Node>" * 130 + b"</Node>" * 130 + b"</HistoryTree>",
    ],
)
def test_bad_history_does_not_break_document_reading(payload):
    data = snapgene_file_to_dict(fileobject=document("GATC", block(7, payload)))
    assert data["seq"] == "GATC"
    assert data["history"]["trees"] == []
    error = data["history"]["parse_errors"][0]
    assert error["parse_error"]
    assert error["raw_block_type"] == 7
    assert error["raw_block_size"] == len(payload)
    assert error["raw_block_offset"] == 19


def test_unknown_compact_encoding_keeps_properties_and_packet_errors():
    unsupported = struct.pack(">IBII", 4, 3, 1, 0) + struct.pack(">BI", 3, 4)
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC",
            tree(1, 4, node(0, 4)),
            snapshot(
                0,
                block(1, unsupported),
                block(
                    30, lzma.compress(block(6, b"<Notes><Custom>kept</Custom></Notes>"))
                ),
            ),
        )
    )["history"]
    assert "Unsupported" in history["nodes"][1]["recovery_reason"]
    assert history["snapshots"][0]["properties"][0]["xml"]["tag"] == "Notes"
    assert history["parse_errors"][0]["packet_type"] == 1


@pytest.mark.parametrize(
    "undo,redo,length",
    [
        ('<Action type="REMOVE" range="3-1"/>', "", 2),
        (
            '<Action type="INSERT" position="5"><Residues type="GENERIC" residues="A"/></Action>',
            "",
            5,
        ),
        ('<Action type="ROTATE" position="2"/>', "", 4),
        (
            '<Action type="REMOVE" range="0-0"/>',
            '<Action type="INSERT" position="0"><Residues type="GENERIC" residues="A"/></Action>',
            3,
        ),
        (
            '<Action type="INSERT" position="0"><Residues type="EXTERNAL" ID="9" length="1"/></Action>',
            "",
            5,
        ),
        ('<Action type="REMOVE" range="0-0"/>', "", 100),
    ],
)
def test_invalid_recovery_never_returns_a_partial_sequence(undo, redo, length):
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC", tree(1, 4, node(0, length)), snapshot(0, manipulation(undo, redo))
        )
    )["history"]
    old = history["nodes"][1]
    assert "sequence" not in old
    assert old["recovery_status"] == "unavailable"
    assert old["recovery_reason"] != "No validated recovery data"


def test_truncated_snapshot_is_localized():
    history = snapgene_file_to_dict(
        fileobject=document(
            "GATC", tree(1, 4, node(0, 4)), block(11, b"\0\0\0\0\x1d\0\0\0\xffshort")
        )
    )["history"]
    assert history["nodes"][0]["sequence"] == "GATC"
    assert "Truncated" in history["snapshots"][0]["parse_error"]


def test_decompression_budget_and_trailing_stream_are_checked(monkeypatch):
    payload = b"<HistoryTree/>"
    monkeypatch.setattr(history_parser, "MAX_TOTAL_DECOMPRESSED_SIZE", 4)
    assert history_parser.parse_history(
        [({"type": 7, "size": len(payload), "offset": 19}, payload)]
    )["parse_errors"]
    monkeypatch.setattr(history_parser, "MAX_TOTAL_DECOMPRESSED_SIZE", 1024)
    compressed = lzma.compress(payload) + b"trailing"
    assert history_parser.parse_history(
        [({"type": 7, "size": len(compressed), "offset": 19}, compressed)]
    )["parse_errors"]


def test_files_without_history_have_a_stable_empty_object():
    assert snapgene_file_to_dict(fileobject=document("GATC"))["history"] == {
        "trees": [],
        "nodes": [],
        "snapshots": [],
        "parse_errors": [],
    }


def test_private_history_fixtures_if_present():
    root = Path(__file__).resolve().parents[1]
    fixtures = [
        ("pAQ1ex_MP1.dna", 15, 14),
        ("\u5355\u4ea4\u6362\u6784\u578b1.dna", 3, 3),
    ]
    for filename, count, available in fixtures:
        path = root / filename
        if not path.exists():
            continue
        data = snapgene_file_to_dict(str(path))
        history = data["history"]
        assert len(history["nodes"]) == count
        assert sum("sequence" in n for n in history["nodes"]) == available
        assert history["nodes"][0]["sequence"] == data["seq"]
        for n in history["nodes"]:
            if "sequence" in n:
                assert len(n["sequence"]) == n["sequence_length"]
                assert "recovery_warning" not in n
        if filename == "pAQ1ex_MP1.dna":
            missing = next(n for n in history["nodes"] if n["id"] == "0")
            assert "encoding: 3" in missing["recovery_reason"]
