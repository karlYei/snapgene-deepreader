import gzip
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


def _encode_bam_sequence(sequence):
    codes = {
        "=": 0,
        "A": 1,
        "C": 2,
        "M": 3,
        "G": 4,
        "R": 5,
        "S": 6,
        "V": 7,
        "T": 8,
        "W": 9,
        "Y": 10,
        "H": 11,
        "K": 12,
        "D": 13,
        "B": 14,
        "N": 15,
    }
    encoded = bytearray()
    for index in range(0, len(sequence), 2):
        high = codes[sequence[index]]
        low = codes[sequence[index + 1]] if index + 1 < len(sequence) else 0
        encoded.append((high << 4) | low)
    return bytes(encoded)


def _minimal_bam():
    header = b"@HD\tVN:1.6\tSO:unsorted\n"
    reference_name = b"reference\x00"
    bam = bytearray()
    bam.extend(b"BAM\1")
    bam.extend(struct.pack("<i", len(header)))
    bam.extend(header)
    bam.extend(struct.pack("<i", 1))
    bam.extend(struct.pack("<i", len(reference_name)))
    bam.extend(reference_name)
    bam.extend(struct.pack("<i", 8))

    qname = b"0\x00"
    cigar_words = [
        (4 << 4) | 0,  # 4M
        (1 << 4) | 1,  # 1I
        (4 << 4) | 0,  # 4M
    ]
    read_sequence = "ACGTTACGT"
    qualities = bytes([255] * len(read_sequence))
    aux = b"NMi" + struct.pack("<i", 1)
    record = bytearray()
    record.extend(
        struct.pack(
            "<iiBBHHHiiii",
            0,
            0,
            len(qname),
            60,
            0,
            len(cigar_words),
            0,
            len(read_sequence),
            -1,
            -1,
            0,
        )
    )
    record.extend(qname)
    record.extend(struct.pack("<" + "I" * len(cigar_words), *cigar_words))
    record.extend(_encode_bam_sequence(read_sequence))
    record.extend(qualities)
    record.extend(aux)
    bam.extend(struct.pack("<i", len(record)))
    bam.extend(record)
    return bytes(bam)


def _write_minimal_file(path, include_rich_blocks=False, include_comparison=False):
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
    if include_comparison:
        parts.append(_block(27, gzip.compress(_minimal_bam())))
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


def test_parses_comparison_alignment_bam_records(tmp_path):
    file_path = tmp_path / "comparison.dna"
    _write_minimal_file(file_path, include_rich_blocks=True, include_comparison=True)
    data = snapgene_file_to_dict(str(file_path))

    assert len(data["comparison_alignments"]) == 1
    comparison = data["comparison_alignments"][0]
    assert comparison["format"] == "BAM"
    assert comparison["references"] == [{"id": 0, "name": "reference", "length": 8}]
    assert len(comparison["records"]) == 1

    record = comparison["records"][0]
    assert record["qname"] == "0"
    assert record["alignable_sequence_name"] == "read-1"
    assert record["ref_name"] == "reference"
    assert record["start"] == 1
    assert record["cigar"] == "4M1I4M"
    assert record["sequence"] == "ACGTTACGT"
    assert record["aux"]["NM"] == 1
    assert record["variants"] == [
        {
            "type": "insertion",
            "after_ref_position": 4,
            "read_position": 5,
            "inserted_bases": "T",
            "length": 1,
        }
    ]


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
    assert data["comparison_alignments"] == []


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
