"""Enhanced SnapGene reader.

This module keeps the public API of ``snapgene-reader`` and adds support for
primer and alignable sequence blocks observed in SnapGene ``*.dna`` files.
"""

import gzip
import re
import struct
import zlib

import html2text
import xmltodict
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .history import parse_history

try:
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    has_dna_alphabet = False


HTML_PARSER = html2text.HTML2Text()
HTML_PARSER.ignore_emphasis = True
HTML_PARSER.ignore_links = True
HTML_PARSER.body_width = 0
HTML_PARSER.single_line_break = True

ZLIB_HEADER_RE = re.compile(b"\x78[\x01\x5e\x9c\xda]")
DNA_BASES = set("ACGTUNacgtun-")
BAM_MAGIC = b"BAM\1"
GZIP_MAGIC = b"\x1f\x8b"
CIGAR_OPS = "MIDNSHP=XB"
BAM_SEQ_CODES = "=ACMGRSVTWYHKDBN"


def parse(val):
    """Parse SnapGene HTML-ish text into plain text."""
    if isinstance(val, str):
        return HTML_PARSER.handle(val).strip().replace("\n", " ").replace('"', "'")
    return val


def parse_dict(obj):
    """Recursively parse text values in dictionaries and lists."""
    if isinstance(obj, dict):
        for key in obj:
            if isinstance(obj[key], str):
                obj[key] = parse(obj[key])
            elif isinstance(obj[key], (dict, list)):
                parse_dict(obj[key])
    elif isinstance(obj, list):
        for item in obj:
            parse_dict(item)
    return obj


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _attr(obj, name, default=None):
    return obj.get("@" + name, default)


def _int_attr(obj, name, default=None):
    value = _attr(obj, name, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _parse_component(component):
    return {
        "hybridized_range": _attr(component, "hybridizedRange"),
        "bases": _attr(component, "bases"),
    }


def _parse_binding_site(site):
    components = [_parse_component(c) for c in _as_list(site.get("Component"))]
    return {
        "location": _attr(site, "location"),
        "bound_strand": _int_attr(site, "boundStrand"),
        "annealed_bases": _attr(site, "annealedBases"),
        "melting_temperature": _int_attr(site, "meltingTemperature"),
        "components": components,
    }


def _parse_primers_block(block_content, block_summary):
    try:
        parsed = parse_dict(xmltodict.parse(block_content.decode("utf-8")))
        primers_root = parsed.get("Primers", {})
        primers = []
        for primer in _as_list(primers_root.get("Primer")):
            binding_sites = [
                _parse_binding_site(site)
                for site in _as_list(primer.get("BindingSite"))
            ]
            primers.append(
                {
                    "id": _int_attr(primer, "recentID"),
                    "name": _attr(primer, "name"),
                    "sequence": _attr(primer, "sequence"),
                    "description": _attr(primer, "description", ""),
                    "date_added": _attr(primer, "dateAdded"),
                    "binding_sites": binding_sites,
                    "raw_block_type": block_summary["type"],
                    "raw_block_size": block_summary["size"],
                    "raw_block_offset": block_summary["offset"],
                }
            )
        return primers
    except Exception as error:
        return [
            {
                "parse_error": str(error),
                "raw_block_type": block_summary["type"],
                "raw_block_size": block_summary["size"],
                "raw_block_offset": block_summary["offset"],
            }
        ]


def _parse_alignable_metadata_block(block_content, block_summary):
    try:
        parsed = parse_dict(xmltodict.parse(block_content.decode("utf-8")))
        root = parsed.get("AlignableSequences", {})
        metadata = {}
        for sequence in _as_list(root.get("Sequence")):
            sequence_id = _attr(sequence, "ID")
            entry = {
                "id": sequence_id,
                "name": _attr(sequence, "name"),
                "sort_order": _int_attr(sequence, "sortOrder"),
                "use": _int_attr(sequence, "use", 1),
                "trimmed_range": _attr(sequence, "trimmedRange"),
                "raw_block_type": block_summary["type"],
                "raw_block_size": block_summary["size"],
                "raw_block_offset": block_summary["offset"],
            }
            metadata[sequence_id] = entry
        return metadata
    except Exception as error:
        return {
            None: {
                "parse_error": str(error),
                "raw_block_type": block_summary["type"],
                "raw_block_size": block_summary["size"],
                "raw_block_offset": block_summary["offset"],
            }
        }


def _zlib_offsets(block_content):
    offsets = []
    if len(block_content) > 36 and block_content[36:38] in (
        b"\x78\x01",
        b"\x78\x5e",
        b"\x78\x9c",
        b"\x78\xda",
    ):
        offsets.append(36)
    for match in ZLIB_HEADER_RE.finditer(block_content):
        offset = match.start()
        if offset not in offsets:
            offsets.append(offset)
    return offsets


def _clean_decompressed_sequence(decompressed):
    text = decompressed.decode("ascii", errors="ignore")
    while text and text[0] not in DNA_BASES:
        text = text[1:]
    return "".join(base for base in text if base in DNA_BASES)


def _parse_alignable_sequence_block(block_content, block_summary):
    entry = {
        "id": str(struct.unpack(">I", block_content[:4])[0])
        if len(block_content) >= 4
        else None,
        "raw_block_type": block_summary["type"],
        "raw_block_size": block_summary["size"],
        "raw_block_offset": block_summary["offset"],
    }

    for offset in _zlib_offsets(block_content):
        try:
            decompressed = zlib.decompress(block_content[offset:])
        except zlib.error:
            continue
        sequence = _clean_decompressed_sequence(decompressed)
        entry.update(
            {
                "sequence": sequence,
                "length": len(sequence),
                "raw_length": len(decompressed),
                "compressed_offset": offset,
            }
        )
        return entry

    entry["parse_error"] = "Could not locate a decompressible zlib sequence stream."
    return entry


def _parse_features_block(block_content, data):
    strand_dict = {"0": ".", "1": "+", "2": "-", "3": "="}
    format_dict = {"@text": parse, "@int": int}
    features_data = xmltodict.parse(block_content)
    features = _as_list(features_data.get("Features", {}).get("Feature"))
    for feature in features:
        segments = _as_list(feature["Segment"])
        segments_ranges = [
            sorted([int(e) for e in segment["@range"].split("-")])
            for segment in segments
        ]
        qualifiers = _as_list(feature.get("Q"))
        parsed_qualifiers = {}
        for qualifier in qualifiers:
            if qualifier["V"] is None:
                continue
            if isinstance(qualifier["V"], list):
                if len(qualifier["V"][0].items()) == 1:
                    parsed_qualifiers[qualifier["@name"]] = l_v = []
                    for e_v in qualifier["V"]:
                        fmt, value = e_v.popitem()
                        fmt = format_dict.get(fmt, parse)
                        l_v.append(fmt(value))
                else:
                    parsed_qualifiers[qualifier["@name"]] = d_v = {}
                    for e_v in qualifier["V"]:
                        (fmt1, value1), (_, value2) = e_v.items()
                        fmt = format_dict.get(fmt1, parse)
                        d_v[value2] = fmt(value1)
            else:
                fmt, value = qualifier["V"].popitem()
                fmt = format_dict.get(fmt, parse)
                parsed_qualifiers[qualifier["@name"]] = fmt(value)

        if "label" not in parsed_qualifiers:
            parsed_qualifiers["label"] = feature["@name"]
        if "note" not in parsed_qualifiers:
            parsed_qualifiers["note"] = []
        if not isinstance(parsed_qualifiers["note"], list):
            parsed_qualifiers["note"] = [parsed_qualifiers["note"]]
        color = segments[0]["@color"]
        parsed_qualifiers["note"].append("color: " + color)

        data["features"].append(
            {
                "start": min([start - 1 for (start, end) in segments_ranges]),
                "end": max([end for (start, end) in segments_ranges]),
                "strand": strand_dict[feature.get("@directionality", "0")],
                "type": feature["@type"],
                "name": feature["@name"],
                "color": color,
                "textColor": "black",
                "segments": segments,
                "row": 0,
                "isOrf": False,
                "qualifiers": parsed_qualifiers,
            }
        )


def _merge_alignable_sequences(sequence_blocks, metadata_by_id):
    merged = []
    seen = set()
    for block_entry in sequence_blocks:
        sequence_id = block_entry.get("id")
        metadata = metadata_by_id.get(sequence_id, {})
        entry = dict(metadata)
        entry.update(block_entry)
        merged.append(entry)
        seen.add(sequence_id)
    for sequence_id, metadata in metadata_by_id.items():
        if sequence_id not in seen:
            merged.append(metadata)
    return merged


def _read_cstring(blob, pos):
    end = blob.index(0, pos)
    return blob[pos:end].decode("utf-8", errors="replace"), end + 1


def _decode_bam_sequence(blob, length):
    bases = []
    for byte in blob:
        bases.append(BAM_SEQ_CODES[byte >> 4])
        if len(bases) == length:
            break
        bases.append(BAM_SEQ_CODES[byte & 0xF])
        if len(bases) == length:
            break
    return "".join(bases)


def _decode_bam_aux(blob):
    tags = {}
    pos = 0
    while pos + 3 <= len(blob):
        tag = blob[pos : pos + 2].decode("ascii", errors="replace")
        typ = chr(blob[pos + 2])
        pos += 3
        try:
            if typ == "A":
                value = chr(blob[pos])
                pos += 1
            elif typ in "cC":
                value = struct.unpack("<b" if typ == "c" else "<B", blob[pos : pos + 1])[0]
                pos += 1
            elif typ in "sS":
                value = struct.unpack("<h" if typ == "s" else "<H", blob[pos : pos + 2])[0]
                pos += 2
            elif typ in "iI":
                value = struct.unpack("<i" if typ == "i" else "<I", blob[pos : pos + 4])[0]
                pos += 4
            elif typ == "f":
                value = struct.unpack("<f", blob[pos : pos + 4])[0]
                pos += 4
            elif typ in "ZH":
                end = blob.index(0, pos)
                value = blob[pos:end].decode("ascii", errors="replace")
                pos = end + 1
            elif typ == "B":
                subtype = chr(blob[pos])
                count = struct.unpack("<I", blob[pos + 1 : pos + 5])[0]
                pos += 5
                sizes = {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}
                size = sizes[subtype]
                raw = blob[pos : pos + count * size]
                pos += len(raw)
                value = {"subtype": subtype, "count": count}
            else:
                value = {"unsupported_type": typ}
                break
        except Exception as error:
            value = {"parse_error": str(error), "type": typ}
            tags[tag] = value
            break
        tags[tag] = value
    return tags


def _cigar_to_string(cigar_ops):
    return "".join(f"{length}{op}" for length, op in cigar_ops)


def _variant_ref_position(ref_index, reference_length):
    if reference_length:
        return (ref_index % reference_length) + 1
    return ref_index + 1


def _call_variants_from_alignment(record, reference_sequence):
    reference_length = len(reference_sequence)
    if not reference_sequence or not record.get("sequence"):
        return []

    variants = []
    read_pos = 0
    ref_pos = record["pos"]
    for length, op in record["cigar_ops"]:
        if op in "M=X":
            for _ in range(length):
                ref_index = ref_pos % reference_length
                ref_base = reference_sequence[ref_index].upper()
                read_base = record["sequence"][read_pos].upper()
                if read_base != "=" and ref_base != read_base:
                    variants.append(
                        {
                            "type": "substitution",
                            "ref_position": ref_index + 1,
                            "read_position": read_pos + 1,
                            "ref_base": ref_base,
                            "read_base": read_base,
                        }
                    )
                ref_pos += 1
                read_pos += 1
        elif op == "I":
            inserted = record["sequence"][read_pos : read_pos + length]
            variants.append(
                {
                    "type": "insertion",
                    "after_ref_position": _variant_ref_position(ref_pos - 1, reference_length),
                    "read_position": read_pos + 1,
                    "inserted_bases": inserted,
                    "length": length,
                }
            )
            read_pos += length
        elif op in "DN":
            deleted = []
            for _ in range(length):
                ref_index = ref_pos % reference_length
                deleted.append(reference_sequence[ref_index].upper())
                ref_pos += 1
            variants.append(
                {
                    "type": "deletion",
                    "ref_start": _variant_ref_position(ref_pos - length, reference_length),
                    "deleted_bases": "".join(deleted),
                    "length": length,
                }
            )
        elif op == "S":
            read_pos += length
        elif op in "HP":
            continue
    return variants


def _parse_bam_payload(bam, block_summary, reference_sequence):
    if not bam.startswith(BAM_MAGIC):
        raise ValueError("Decompressed payload is not BAM.")

    pos = 4
    header_length = struct.unpack("<i", bam[pos : pos + 4])[0]
    pos += 4
    header_text = bam[pos : pos + header_length].decode("utf-8", errors="replace")
    pos += header_length
    reference_count = struct.unpack("<i", bam[pos : pos + 4])[0]
    pos += 4

    references = []
    for ref_id in range(reference_count):
        name_length = struct.unpack("<i", bam[pos : pos + 4])[0]
        pos += 4
        name = bam[pos : pos + name_length - 1].decode("utf-8", errors="replace")
        pos += name_length
        length = struct.unpack("<i", bam[pos : pos + 4])[0]
        pos += 4
        references.append({"id": ref_id, "name": name, "length": length})

    records = []
    while pos < len(bam):
        block_size = struct.unpack("<i", bam[pos : pos + 4])[0]
        pos += 4
        record_blob = bam[pos : pos + block_size]
        pos += block_size
        if len(record_blob) < 32:
            continue

        (
            ref_id,
            ref_pos,
            read_name_length,
            mapq,
            bin_value,
            cigar_count,
            flag,
            read_length,
            next_ref_id,
            next_pos,
            template_length,
        ) = struct.unpack("<iiBBHHHiiii", record_blob[:32])
        cursor = 32
        read_name = record_blob[cursor : cursor + read_name_length - 1].decode(
            "utf-8", errors="replace"
        )
        cursor += read_name_length

        cigar_words = struct.unpack(
            "<" + "I" * cigar_count, record_blob[cursor : cursor + 4 * cigar_count]
        )
        cursor += 4 * cigar_count
        cigar_ops = [(word >> 4, CIGAR_OPS[word & 0xF]) for word in cigar_words]
        encoded_sequence_length = (read_length + 1) // 2
        sequence = _decode_bam_sequence(
            record_blob[cursor : cursor + encoded_sequence_length], read_length
        )
        cursor += encoded_sequence_length
        qualities = list(record_blob[cursor : cursor + read_length])
        cursor += read_length
        aux = _decode_bam_aux(record_blob[cursor:])

        ref = references[ref_id] if 0 <= ref_id < len(references) else None
        record = {
            "qname": read_name,
            "read_id": read_name,
            "ref_id": ref_id,
            "ref_name": ref["name"] if ref else None,
            "pos": ref_pos,
            "start": ref_pos + 1,
            "mapq": mapq,
            "bin": bin_value,
            "flag": flag,
            "read_length": read_length,
            "cigar": _cigar_to_string(cigar_ops),
            "cigar_ops": cigar_ops,
            "sequence": sequence,
            "qualities": qualities,
            "next_ref_id": next_ref_id,
            "next_pos": next_pos,
            "template_length": template_length,
            "aux": aux,
            "variants": [],
            "raw_block_type": block_summary["type"],
            "raw_block_size": block_summary["size"],
            "raw_block_offset": block_summary["offset"],
        }
        record["variants"] = _call_variants_from_alignment(record, reference_sequence)
        records.append(record)

    return {
        "format": "BAM",
        "header": header_text,
        "references": references,
        "records": records,
        "raw_block_type": block_summary["type"],
        "raw_block_size": block_summary["size"],
        "raw_block_offset": block_summary["offset"],
    }


def _gzip_payloads(block_content):
    offsets = []
    start = 0
    while True:
        offset = block_content.find(GZIP_MAGIC, start)
        if offset < 0:
            break
        offsets.append(offset)
        start = offset + 1
    if not offsets and block_content.startswith(BAM_MAGIC):
        return [(0, block_content)]

    payloads = []
    for offset in offsets:
        try:
            payloads.append((offset, gzip.decompress(block_content[offset:])))
        except Exception:
            continue
    return payloads


def _parse_comparison_block(block_content, block_summary, reference_sequence):
    comparisons = []
    for compressed_offset, payload in _gzip_payloads(block_content):
        if not payload.startswith(BAM_MAGIC):
            continue
        try:
            parsed = _parse_bam_payload(payload, block_summary, reference_sequence)
            parsed["compressed_offset"] = compressed_offset
            comparisons.append(parsed)
        except Exception as error:
            comparisons.append(
                {
                    "parse_error": str(error),
                    "compressed_offset": compressed_offset,
                    "raw_block_type": block_summary["type"],
                    "raw_block_size": block_summary["size"],
                    "raw_block_offset": block_summary["offset"],
                }
            )
    if not comparisons:
        comparisons.append(
            {
                "parse_error": "No BAM payload found in comparison block.",
                "raw_block_type": block_summary["type"],
                "raw_block_size": block_summary["size"],
                "raw_block_offset": block_summary["offset"],
            }
        )
    return comparisons


def _enrich_comparison_records(comparison_alignments, metadata_by_id):
    for comparison in comparison_alignments:
        for record in comparison.get("records", []):
            metadata = metadata_by_id.get(record.get("qname"))
            if metadata:
                record["alignable_sequence_name"] = metadata.get("name")
                record["alignable_sequence_trimmed_range"] = metadata.get("trimmed_range")
                record["alignable_sequence_use"] = metadata.get("use")


def snapgene_file_to_dict(filepath=None, fileobject=None):
    """Return a dictionary containing data from a SnapGene ``*.dna`` file."""
    close_file = False
    if filepath is not None:
        fileobject = open(filepath, "rb")
        close_file = True
    if fileobject is None:
        raise ValueError("Either filepath or fileobject must be provided.")

    try:
        if fileobject.read(1) != b"\t":
            raise ValueError("Wrong format for a SnapGene file!")

        def unpack(size, mode):
            return struct.unpack(">" + mode, fileobject.read(size))[0]

        length = unpack(4, "I")
        title = fileobject.read(8).decode("ascii")
        if length != 14 or title != "SnapGene":
            raise ValueError("Wrong format for a SnapGene file!")

        data = {
            "isDNA": unpack(2, "H"),
            "exportVersion": unpack(2, "H"),
            "importVersion": unpack(2, "H"),
            "features": [],
            "primers": [],
            "alignable_sequences": [],
            "comparison_alignments": [],
            "raw_blocks_summary": [],
        }
        alignable_sequence_blocks = []
        alignable_metadata = {}
        history_blocks = []

        while True:
            block_offset = fileobject.tell()
            next_byte = fileobject.read(1)
            if next_byte == b"":
                break

            block_size = unpack(4, "I")
            block_type = ord(next_byte)
            block_content = fileobject.read(block_size)
            block_summary = {
                "type": block_type,
                "size": block_size,
                "offset": block_offset,
            }
            data["raw_blocks_summary"].append(block_summary)

            if block_type == 0:
                props = struct.unpack(">b", block_content[:1])[0]
                data["dna"] = {
                    "topology": "circular" if props & 0x01 else "linear",
                    "strandedness": "double" if props & 0x02 > 0 else "single",
                    "damMethylated": props & 0x04 > 0,
                    "dcmMethylated": props & 0x08 > 0,
                    "ecoKIMethylated": props & 0x10 > 0,
                    "length": block_size - 1,
                }
                data["seq"] = block_content[1:].decode("ascii")
            elif block_type == 5:
                data["primers"].extend(
                    _parse_primers_block(block_content, block_summary)
                )
            elif block_type == 6:
                block_text = block_content.decode("utf-8")
                note_data = parse_dict(xmltodict.parse(block_text))
                data["notes"] = note_data["Notes"]
            elif block_type in (7, 11):
                history_blocks.append((block_summary, block_content))
            elif block_type == 10:
                _parse_features_block(block_content, data)
            elif block_type == 16:
                alignable_sequence_blocks.append(
                    _parse_alignable_sequence_block(block_content, block_summary)
                )
            elif block_type == 17:
                alignable_metadata.update(
                    _parse_alignable_metadata_block(block_content, block_summary)
                )
            elif block_type == 27:
                data["comparison_alignments"].extend(
                    _parse_comparison_block(
                        block_content,
                        block_summary,
                        data.get("seq", ""),
                    )
                )

        data["alignable_sequences"] = _merge_alignable_sequences(
            alignable_sequence_blocks, alignable_metadata
        )
        _enrich_comparison_records(data["comparison_alignments"], alignable_metadata)
        data["history"] = parse_history(
            history_blocks, data.get("seq", ""), data.get("dna", {}).get("topology")
        )
        return data
    finally:
        if close_file:
            fileobject.close()


def snapgene_file_to_seqrecord(filepath=None, fileobject=None):
    """Return a BioPython SeqRecord from a SnapGene ``*.dna`` file."""
    data = snapgene_file_to_dict(filepath=filepath, fileobject=fileobject)
    strand_dict = {"+": 1, "-": -1, ".": 0, "=": 0}

    if has_dna_alphabet:
        seq = Seq(data["seq"], alphabet=DNAAlphabet())
    else:
        seq = Seq(data["seq"])

    seqrecord = SeqRecord(
        seq=seq,
        features=[
            SeqFeature(
                location=FeatureLocation(
                    start=feature["start"],
                    end=feature["end"],
                    strand=strand_dict[feature["strand"]],
                ),
                type=feature["type"],
                qualifiers=feature["qualifiers"],
            )
            for feature in data["features"]
        ],
        annotations={"topology": data["dna"]["topology"], **data.get("notes", {})},
    )
    seqrecord.annotations["molecule_type"] = "DNA"
    return seqrecord


def _analyse_gs(dic, *args, **kwargs):
    default = kwargs.get("default")
    for arg in args:
        if isinstance(dic, dict) and arg in dic:
            dic = dic[arg]
        else:
            return default
    return dic


def snapgene_file_to_gbk(read_file_object, write_file_object):
    """Convert a SnapGene file object to a GenBank file object."""
    data = snapgene_file_to_dict(fileobject=read_file_object)
    wfo = write_file_object
    wfo.write(
        (
            "LOCUS       Exported              {0:>6} bp ds-DNA     "
            "{1:>8} SYN 15-APR-2012\n"
        ).format(len(data["seq"]), data["dna"]["topology"])
    )
    definition = _analyse_gs(data, "notes", "Description", default=".").replace(
        "\n", "\n            "
    )
    wfo.write("DEFINITION  {}\n".format(definition))
    wfo.write("ACCESSION   .\n")
    wfo.write("VERSION     .\n")
    wfo.write(
        "KEYWORDS    {}\n".format(
            _analyse_gs(data, "notes", "CustomMapLabel", default=".")
        )
    )
    wfo.write("SOURCE      .\n")
    wfo.write("  ORGANISM  .\n")

    references = _analyse_gs(data, "notes", "References")
    reference_count = 0
    if references:
        for key in references:
            reference_count += 1
            ref = references[key]
            wfo.write(
                "REFERENCE   {}  (bases 1 to {} )\n".format(
                    reference_count, _analyse_gs(data, "dna", "length")
                )
            )
            for key2 in ref:
                gb_key = key2.replace("@", "").upper()
                wfo.write("  {}   {}\n".format(gb_key, ref[key2]))

    reference_count += 1
    wfo.write(
        "REFERENCE   {}  (bases 1 to {} )\n".format(
            reference_count, _analyse_gs(data, "dna", "length")
        )
    )
    wfo.write("  AUTHORS   SnapGeneDeepReader\n")
    wfo.write("  TITLE     Direct Submission\n")
    wfo.write("  JOURNAL   Exported from SnapGene File Reader\n")
    wfo.write("            https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader\n")
    wfo.write(
        "COMMENT     {}\n".format(
            _analyse_gs(data, "notes", "Comments", default=".")
            .replace("\n", "\n            ")
            .replace("\\", "")
        )
    )
    wfo.write("FEATURES             Location/Qualifiers\n")

    for feature in _analyse_gs(data, "features", default=[]):
        strand = _analyse_gs(feature, "strand", default="")
        segments = [
            x
            for x in _analyse_gs(feature, "segments", default=[])
            if x["@type"] == "standard"
        ]
        if len(segments) > 1:
            line = "join("
            for segment in segments:
                line += _analyse_gs(segment, "@range").replace("-", "..") + ","
            line = line[:-1] + ")"
        else:
            line = "{}..{}".format(
                _analyse_gs(feature, "start", default=" "),
                _analyse_gs(feature, "end", default=" "),
            )

        if strand == "-":
            wfo.write(
                "     {} complement({})\n".format(
                    _analyse_gs(feature, "type", default=" ").ljust(15), line
                )
            )
        else:
            wfo.write(
                "     {} {}\n".format(
                    _analyse_gs(feature, "type", default=" ").ljust(15), line
                )
            )

        for qualifier_name, qualifier_value in _analyse_gs(
            feature, "qualifiers", default={}
        ).items():
            if isinstance(qualifier_value, list):
                values = qualifier_value
            else:
                values = [qualifier_value]
            for value in values:
                wfo.write(
                    "                     /{}=\"{}\"\n".format(
                        qualifier_name, str(value).replace('"', "'")
                    )
                )

    wfo.write("ORIGIN\n")
    sequence = data["seq"].lower()
    for line_start in range(0, len(sequence), 60):
        line_sequence = sequence[line_start : line_start + 60]
        groups = [
            line_sequence[group_start : group_start + 10]
            for group_start in range(0, len(line_sequence), 10)
        ]
        wfo.write("{:>9} {}\n".format(line_start + 1, " ".join(groups)))
    wfo.write("//\n")
