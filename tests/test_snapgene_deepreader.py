import struct
import zlib
from io import StringIO
from pathlib import Path

from snapgene_deepreader import (
    snapgene_file_to_dict,
    snapgene_file_to_gbk,
    snapgene_file_to_seqrecord,
)


ROOT = Path(__file__).resolve().parents[1]
PAQ1EX = ROOT / "pAQ1ex_MP1.dna"
SINGLE_EXCHANGE = ROOT / "单交换构型1.dna"


def _block(block_type, content):
    return bytes([block_type]) + struct.pack(">I", len(content)) + content


def _snapgene_header():
    return b"\x09" + struct.pack(">I", 14) + b"SnapGene" + struct.pack(">HHH", 1, 1, 1)


def _write_minimal_file(path, include_rich_blocks=False):
    parts = [
        _snapgene_header(),
        _block(0, b"\x01ACGTACGT"),
        _block(6, b"<Notes><Type>Natural</Type></Notes>"),
    ]
    if include_rich_blocks:
        primers = (
            b'<?xml version="1.0"?><Primers nextValidID="2">'
            b'<Primer recentID="1" name="test-primer" sequence="ACGTAC" '
            b'description="&lt;html&gt;&lt;body&gt;demo&lt;/body&gt;&lt;/html&gt;" '
            b'dateAdded="2026-08-04T00:00:00Z">'
            b'<BindingSite location="1-6" boundStrand="0" annealedBases="ACGTAC" '
            b'meltingTemperature="42"><Component hybridizedRange="1-6" bases="ACGTAC"/>'
            b"</BindingSite></Primer></Primers>"
        )
        align_metadata = (
            b'<AlignableSequences trimStringency="Lowest">'
            b'<Sequence name="read-1" ID="0" sortOrder="1" trimmedRange="0..8"/>'
            b"</AlignableSequences>"
        )
        compressed = zlib.compress(b"\x00ACGTACGT")
        align_sequence = struct.pack(">I", 0) + (b"\x00" * 32) + compressed
        parts.extend([_block(16, align_sequence), _block(17, align_metadata), _block(5, primers)])
    path.write_bytes(b"".join(parts))


def test_reads_original_fields_and_rich_blocks(tmp_path):
    file_path = tmp_path / "rich.dna"
    _write_minimal_file(file_path, include_rich_blocks=True)
    data = snapgene_file_to_dict(str(file_path))

    assert data["seq"] == "ACGTACGT"
    assert data["dna"]["topology"] == "circular"
    assert data["features"] == []
    assert len(data["primers"]) == 1
    assert len(data["alignable_sequences"]) == 1
    assert data["raw_blocks_summary"]


def test_parses_primers(tmp_path):
    file_path = tmp_path / "rich.dna"
    _write_minimal_file(file_path, include_rich_blocks=True)
    data = snapgene_file_to_dict(str(file_path))

    primer = data["primers"][0]
    assert primer["name"] == "test-primer"
    assert primer["sequence"] == "ACGTAC"
    assert primer["binding_sites"]
    assert primer["binding_sites"][0]["location"] == "1-6"
    assert primer["binding_sites"][0]["components"]


def test_parses_alignable_sequences_with_decompressed_bases(tmp_path):
    file_path = tmp_path / "rich.dna"
    _write_minimal_file(file_path, include_rich_blocks=True)
    data = snapgene_file_to_dict(str(file_path))

    alignable = data["alignable_sequences"][0]
    assert alignable["name"] == "read-1"
    assert alignable["id"] == "0"
    assert alignable["trimmed_range"] == "0..8"
    assert alignable["compressed_offset"] == 36
    assert alignable["length"] == 8
    assert alignable["sequence"] == "ACGTACGT"
    assert set(alignable["sequence"]).issubset(set("ACGTUNacgtun-"))


def test_local_paq1ex_fixture_if_present():
    if not PAQ1EX.exists():
        return
    data = snapgene_file_to_dict(str(PAQ1EX))

    assert len(data["seq"]) == 7134
    assert data["dna"]["topology"] == "circular"
    assert len(data["features"]) == 11
    assert len(data["primers"]) == 11
    assert len(data["alignable_sequences"]) == 6


def test_local_second_fixture_if_present():
    if not SINGLE_EXCHANGE.exists():
        return
    data = snapgene_file_to_dict(str(SINGLE_EXCHANGE))

    assert len(data["seq"]) == 19787
    assert data["dna"]["topology"] == "linear"
    assert len(data["features"]) == 56


def test_files_without_rich_blocks_get_empty_lists(tmp_path):
    file_path = tmp_path / "minimal.dna"
    _write_minimal_file(file_path)

    data = snapgene_file_to_dict(str(file_path))

    assert data["seq"] == "ACGTACGT"
    assert data["primers"] == []
    assert data["alignable_sequences"] == []


def test_seqrecord_compatibility(tmp_path):
    file_path = tmp_path / "minimal.dna"
    _write_minimal_file(file_path)
    record = snapgene_file_to_seqrecord(str(file_path))

    assert len(record.seq) == 8
    assert record.annotations["topology"] == "circular"
    assert len(record.features) == 0


def test_gbk_export_compatibility(tmp_path):
    file_path = tmp_path / "minimal.dna"
    _write_minimal_file(file_path)
    output = StringIO()
    with file_path.open("rb") as handle:
        snapgene_file_to_gbk(handle, output)

    gbk = output.getvalue()
    assert gbk.startswith("LOCUS")
    assert "FEATURES" in gbk
    assert "ORIGIN" in gbk
    assert gbk.rstrip().endswith("//")
