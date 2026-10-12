import json
import lzma
import struct
from io import BytesIO, StringIO

import pytest
from Bio import SeqIO

from snapgene_deepreader import (
    snapgene_file_to_dict,
    snapgene_file_to_seqrecord,
    snapgene_file_to_gbk,
)


def block(kind, payload):
    return bytes([kind]) + struct.pack(">I", len(payload)) + payload


def document(*packets):
    return BytesIO(
        block(9, b"SnapGene" + struct.pack(">HHH", 1, 1, 1))
        + block(0, b"\x03AAAA")
        + b"".join(packets)
        + block(6, b"<Notes><Type>Natural</Type></Notes>")
    )


def feature(name="test", qualifiers="", extra=""):
    return (
        '<Features><Feature name="%s" type="gene" %s>'
        '<Segment range="1-4" color="#00aa00" type="standard"/>%s</Feature></Features>'
        % (name, extra, qualifiers)
    )


def node(node_id, properties="", inputs="", extra=""):
    return (
        '<Node ID="%s" name="node%s" type="DNA" circular="1" '
        'seqLen="4" operation="amplifyFragment" %s>%s%s</Node>'
        % (node_id, node_id, extra, properties, inputs)
    )


def tree(xml):
    return block(7, lzma.compress(("<HistoryTree>" + xml + "</HistoryTree>").encode()))


def snapshot(node_id, properties):
    return block(
        11,
        struct.pack(">I", node_id)
        + block(1, struct.pack(">IBIIBI", 4, 3, 1, 0, 1, 4) + b"\x55")
        + block(
            30,
            lzma.compress(
                b"".join(block(kind, xml.encode()) for kind, xml in properties)
            ),
        ),
    )


@pytest.mark.parametrize(
    "attrs",
    [
        'bool="1" text="1.05" predef="ab initio prediction"',
        'predef="ab initio prediction" bool="1" text="1.05"',
        'text="1.05" predef="ab initio prediction" bool="1"',
    ],
)
def test_complex_qualifiers_are_lossless_and_attribute_order_independent(attrs):
    qualifiers = (
        '<Q name="inference"><V %s/><V bool="1" text="2.13" predef="ab initio prediction"/></Q>'
        % attrs
    )
    stream = document(block(10, feature(qualifiers=qualifiers).encode()))
    data = snapgene_file_to_dict(fileobject=stream)
    annotation = data["features"][0]
    assert annotation["qualifiers"]["inference"] == [
        {"bool": "1", "text": "1.05", "predef": "ab initio prediction"},
        {"bool": "1", "text": "2.13", "predef": "ab initio prediction"},
    ]
    assert (
        annotation["qualifier_values"][0]["V"][0]["@predef"] == "ab initio prediction"
    )
    stream.seek(0)
    record = snapgene_file_to_seqrecord(fileobject=stream)
    assert json.loads(record.features[0].qualifiers["inference"][0])["bool"] == "1"
    stream.seek(0)
    output = StringIO()
    snapgene_file_to_gbk(stream, output)
    output.seek(0)
    exported = SeqIO.read(output, "genbank")
    assert str(exported.seq) == "AAAA"
    assert exported.features[0].location.start == record.features[0].location.start
    assert exported.features[0].location.end == record.features[0].location.end
    assert json.loads(exported.features[0].qualifiers["inference"][1])["text"] == "2.13"


def test_mixed_and_duplicate_qualifiers_do_not_drop_fields():
    qualifiers = (
        '<Q name="mixed"><V text="simple"/><V text="complex" extra="kept"/></Q>'
        '<Q name="legacy"><V predef="one" int="1"/><V text="two" predef="two"/></Q>'
        '<Q name="duplicates"><V text="one" predef="same"/><V text="two" predef="same"/></Q>'
    )
    data = snapgene_file_to_dict(
        fileobject=document(block(10, feature(qualifiers=qualifiers).encode()))
    )
    q = data["features"][0]["qualifiers"]
    assert q["mixed"] == ["simple", {"text": "complex", "extra": "kept"}]
    assert q["legacy"] == {"one": 1, "two": "two"}
    assert q["duplicates"] == [{"same": "one"}, {"same": "two"}]


def primer_xml(location, sequence="AAAA", name="used", extra=""):
    return (
        '<Primers><HybridizationParams minContinuousMatchLen="8"/>'
        '<Primer name="%s" sequence="%s" recentID="1" custom="kept">'
        '<BindingSite location="%s" boundStrand="0" simplified="1" %s>'
        '<Component bases="AAAA" hybridizedRange="%s"/>'
        "</BindingSite></Primer></Primers>"
        % (name, sequence, location, extra, location)
    )


def test_snapshot_annotations_tree_selection_and_document_root_are_separate():
    source = node(1, feature("tree feature", extra='visible="0"') + primer_xml("3-0"))
    root = node(
        2,
        '<Oligo name="used" sequence="AAAA" description="operation primer"/>'
        '<Parameter name="polymerase" val="blunt" extra="kept"/>'
        '<InputSummary manipulation="amplify" val1="3" val2="0"/>'
        '<HistoryColors><StrandColors type="Input"/></HistoryColors>',
        source,
    )
    data = snapgene_file_to_dict(
        fileobject=document(
            block(10, feature("current feature").encode()),
            block(5, primer_xml("0-3", name="attached only").encode()),
            snapshot(
                1,
                [
                    (10, feature("snapshot feature")),
                    (5, primer_xml("1-2")),
                    (6, "<Notes><Comments>source notes</Comments></Notes>"),
                ],
            ),
            tree(root),
        )
    )
    current, old = data["history"]["nodes"]
    assert current["features"][0]["name"] == "current feature"
    assert current["primers"][0]["name"] == "attached only"
    assert current["annotation_sources"]["features"] == "document"
    assert old["features"][0]["name"] == "snapshot feature"
    assert old["tree_features"][0]["name"] == "tree feature"
    assert old["tree_features"][0]["visible"] is False
    assert old["notes"]["Comments"] == "source notes"
    assert old["annotation_sources"]["primers"] == "snapshot"
    assert old["primers"][0]["attributes"]["custom"] == "kept"
    used = current["used_primers"][0]
    assert used["source"] == "history_tree_oligo"
    assert used["description"] == "operation primer"
    assert used["binding_sites_status"] == "recorded"
    site = used["binding_sites"][0]
    assert site["template_key"] == old["key"]
    assert site["location"] == "3-0"
    assert site["region"] == {
        "start": 3,
        "end": 0,
        "wraps_origin": True,
        "coordinate_system": "zero_based_inclusive",
    }
    assert site["source"] == "history_tree"
    assert site["components"][0]["region"] == site["region"]
    assert current["parameters"][0]["extra"] == "kept"
    assert current["inputs"][0]["summary"]["attributes"]["val1"] == "3"
    assert old["consumer_keys"] == [current["key"]]
    assert current["history_colors"][0]["tag"] == "HistoryColors"
    assert data["history"]["parse_errors"] == []
    json.dumps(data)


def test_missing_used_primer_regions_are_not_fabricated_or_inherited():
    source = node(1, primer_xml("0-3", sequence="CCCC"))
    root = node(2, '<Oligo name="used" sequence="AAAA"/>', source)
    current, old = snapgene_file_to_dict(fileobject=document(tree(root)))["history"][
        "nodes"
    ]
    assert old["used_primers"] == []
    assert current["used_primers"][0]["binding_sites"] == []
    assert current["used_primers"][0]["binding_sites_status"] == "unavailable"
    assert current["primers"] == []


def test_used_primer_regions_fall_back_to_matched_snapshot():
    root = node(2, '<Oligo name="alias" sequence="aaaa"/>', node(1))
    data = snapgene_file_to_dict(
        fileobject=document(tree(root), snapshot(1, [(5, primer_xml("0-3"))]))
    )
    used = data["history"]["nodes"][0]["used_primers"][0]
    assert used["binding_sites"][0]["source"] == "snapshot"
    assert used["binding_sites"][0]["matched_primer_name"] == "used"


def test_hidden_history_is_branch_local_and_keeps_boundary_node_visible():
    left = node(
        1,
        inputs=node(0, inputs=node(4), extra='hideHistory="1"'),
        extra='hideHistory="1"',
    )
    right = node(3, inputs=node(0), extra='hideHistory="0"')
    data = snapgene_file_to_dict(
        fileobject=document(tree(node(2, inputs=left + right)))
    )
    root, left, older, oldest, right, right_old = data["history"]["nodes"]
    assert root["is_visible"] is True
    assert left["hide_history"] is True
    assert left["is_visible"] is True
    assert older["is_visible"] is False
    assert older["hidden_by"] == [left["key"]]
    assert oldest["hidden_by"] == [left["key"], older["key"]]
    assert right["is_visible"] is True
    assert right_old["is_visible"] is True
    assert data["history"]["trees"][0]["visible_node_keys"] == [
        root["key"],
        left["key"],
        right["key"],
        right_old["key"],
    ]
    assert data["history"]["trees"][0]["hidden_node_keys"] == [
        older["key"],
        oldest["key"],
    ]


def test_unknown_visibility_flag_and_local_annotation_errors_are_inspectable():
    root = node(
        2,
        '<Features><Feature name="incomplete"/></Features>',
        node(1),
        extra='hideHistory="future_value"',
    )
    data = snapgene_file_to_dict(fileobject=document(tree(root)))
    current, old = data["history"]["nodes"]
    assert current["sequence"] == "AAAA"
    assert current["tree_features"][0]["parse_error"]
    assert current["annotation_errors"]
    assert current["hide_history"] is None
    assert current["is_visible"] is True
    assert old["is_visible"] is None
    assert old["visibility_unknown_by"] == [current["key"]]
    assert data["history"]["parse_errors"] == []


def test_empty_multi_root_tree_and_unmatched_input_summaries():
    empty = snapgene_file_to_dict(fileobject=document(tree("")))["history"]
    assert empty["trees"][0]["visible_node_keys"] == []
    roots = node(1, '<InputSummary manipulation="amplify"/>', node(2) + node(3)) + node(
        4
    )
    history = snapgene_file_to_dict(fileobject=document(tree(roots)))["history"]
    assert len(history["trees"][0]["visible_node_keys"]) == 4
    assert all(i["summary"] is None for i in history["nodes"][0]["inputs"])


def test_empty_primer_records_are_known_empty_not_parse_failures():
    data = snapgene_file_to_dict(
        fileobject=document(
            block(5, b"<Primers/>"),
            tree(node(2, "<Primers/>", node(1))),
            snapshot(1, [(5, "<Primers/>")]),
        )
    )
    assert data["primers"] == []
    current, old = data["history"]["nodes"]
    assert current["tree_primers"] == []
    assert old["primers"] == []
    assert old["annotation_sources"]["primers"] == "snapshot"
    assert not any(n["annotation_errors"] for n in data["history"]["nodes"])


def test_visibility_and_annotations_do_not_leak_between_trees():
    history = snapgene_file_to_dict(
        fileobject=document(
            tree(node(2, inputs=node(1, feature("first")), extra='hideHistory="1"')),
            tree(node(2, inputs=node(1, feature("second")))),
        )
    )["history"]
    first, first_input, second, second_input = history["nodes"]
    assert first_input["hidden_by"] == [first["key"]]
    assert second_input["is_visible"] is True
    assert second["hide_history"] is False
    assert first_input["features"][0]["name"] == "first"
    assert second_input["features"][0]["name"] == "second"
    assert history["trees"][0]["visible_node_keys"] == [first["key"]]
    assert history["trees"][1]["visible_node_keys"] == [
        second["key"],
        second_input["key"],
    ]
