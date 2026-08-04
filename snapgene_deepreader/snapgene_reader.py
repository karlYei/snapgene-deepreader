"""Enhanced SnapGene reader.

This module keeps the public API of ``snapgene-reader`` and adds support for
primer and alignable sequence blocks observed in SnapGene ``*.dna`` files.
"""

import re
import struct
import zlib

import html2text
import xmltodict
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

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
            "raw_blocks_summary": [],
        }
        alignable_sequence_blocks = []
        alignable_metadata = {}

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

        data["alignable_sequences"] = _merge_alignable_sequences(
            alignable_sequence_blocks, alignable_metadata
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
