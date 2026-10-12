"""Read the saved HistoryTree and validated history recovery packet formats."""

import lzma
import struct
from xml.etree import ElementTree as ET

MAX_DECOMPRESSED_SIZE = 64 * 1024 * 1024
MAX_TOTAL_DECOMPRESSED_SIZE = 128 * 1024 * 1024
MAX_XML_DEPTH = 128
MAX_XML_ELEMENTS = 100000
XZ_MAGIC = b"\xfd7zXZ\x00"
COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def empty_history():
    return {"trees": [], "nodes": [], "snapshots": [], "parse_errors": []}


def _source(summary):
    return {
        "raw_block_type": summary["type"],
        "raw_block_size": summary["size"],
        "raw_block_offset": summary["offset"],
    }


def _decompress(payload, budget):
    if payload.startswith(XZ_MAGIC):
        decoder = lzma.LZMADecompressor(memlimit=MAX_DECOMPRESSED_SIZE)
        payload = decoder.decompress(payload, max_length=MAX_DECOMPRESSED_SIZE + 1)
        if len(payload) > MAX_DECOMPRESSED_SIZE:
            raise ValueError("History decompression size limit exceeded")
        if not decoder.eof or decoder.unused_data:
            raise ValueError("Incomplete XZ stream or trailing history data")
    if len(payload) > MAX_DECOMPRESSED_SIZE or len(payload) > budget[0]:
        raise ValueError("History decompression size limit exceeded")
    budget[0] -= len(payload)
    return payload


def _xml(payload):
    # ElementTree expands internal entities; history has no need for DTDs.
    text = payload.decode("utf-8-sig")
    if "<!DOCTYPE" in text.upper() or "<!ENTITY" in text.upper():
        raise ValueError("DTD/entity declarations are not supported in history XML")
    root = ET.fromstring(text)
    stack = [(root, 0)]
    count = 0
    while stack:
        element, depth = stack.pop()
        count += 1
        if depth > MAX_XML_DEPTH or count > MAX_XML_ELEMENTS:
            raise ValueError("History XML depth/element limit exceeded")
        stack.extend((child, depth + 1) for child in element)
    return root


def _element(element):
    result = {
        "tag": element.tag,
        "attributes": dict(element.attrib),
        "children": [_element(child) for child in element],
    }
    if element.text and element.text.strip():
        result["text"] = element.text
    return result


def _packets(payload):
    offset = 0
    while offset < len(payload):
        if len(payload) - offset < 5:
            raise ValueError("Truncated history packet header")
        size = struct.unpack_from(">I", payload, offset + 1)[0]
        end = offset + 5 + size
        if end > len(payload):
            raise ValueError("Truncated history packet payload")
        yield {"type": payload[offset], "size": size, "offset": offset}, payload[
            offset + 5 : end
        ]
        offset = end


def _compact_sequence(payload, budget):
    """Decode the observed run-1 (two-bit G/A/T/C) sequence encoding."""
    if len(payload) < 13:
        raise ValueError("Truncated compact history sequence")
    length, flags, run_count, lowercase_count = struct.unpack_from(">IBII", payload)
    if length > MAX_DECOMPRESSED_SIZE or run_count > len(payload) // 5:
        raise ValueError("Invalid compact sequence length/run count")
    if length > budget[0]:
        raise ValueError("History sequence size limit exceeded")
    budget[0] -= length
    offset = 13
    parts = []
    decoded_length = 0
    for _ in range(run_count):
        if offset + 5 > len(payload):
            raise ValueError("Truncated compact sequence run")
        encoding, count = struct.unpack_from(">BI", payload, offset)
        offset += 5
        if encoding != 1:
            raise ValueError("Unsupported compact sequence run encoding: %s" % encoding)
        size = (count + 3) // 4
        decoded_length += count
        if decoded_length > length or offset + size > len(payload):
            raise ValueError("Invalid compact sequence run length")
        packed = payload[offset : offset + size]
        full_bytes = count // 4
        bases = "".join(
            "GATC"[(byte >> shift) & 3]
            for byte in packed[:full_bytes]
            for shift in (6, 4, 2, 0)
        )
        # The final partial group is right-aligned, unlike a full group.
        if count % 4:
            bases += "".join(
                "GATC"[(packed[-1] >> shift) & 3]
                for shift in range(2 * (count % 4 - 1), -1, -2)
            )
        parts.append(bases)
        offset += size
    if decoded_length != length or offset + lowercase_count * 8 != len(payload):
        raise ValueError("Compact sequence length/trailing data mismatch")
    sequence = list("".join(parts))
    for _ in range(lowercase_count):
        start, end = struct.unpack_from(">II", payload, offset)
        offset += 8
        if not 0 <= start <= end < length:
            raise ValueError("Invalid compact sequence lowercase range")
        sequence[start : end + 1] = [base.lower() for base in sequence[start : end + 1]]
    return "".join(sequence), flags


def _tree(payload, summary, tree_index):
    root = _xml(payload)
    if root.tag != "HistoryTree":
        raise ValueError("Expected HistoryTree XML")
    tree = dict(
        _source(summary),
        attributes=dict(root.attrib),
        root_keys=[],
        raw_xml=payload.decode("utf-8-sig"),
    )
    nodes = []

    def visit(element):
        attributes = dict(element.attrib)
        key = "%s:%s" % (tree_index, len(nodes))
        length = attributes.get("seqLen")
        node = {
            "key": key,
            "id": attributes.get("ID"),
            "name": attributes.get("name"),
            "operation": attributes.get("operation"),
            "sequence_length": None,
            "topology": {"1": "circular", "0": "linear"}.get(
                attributes.get("circular")
            ),
            "strandedness": attributes.get("strandedness"),
            "attributes": attributes,
            "input_keys": [],
            "input_summaries": [],
            "properties": [],
            "recovery_status": "unavailable",
            "recovery_reason": "No validated recovery data",
        }
        if length is not None:
            try:
                node["sequence_length"] = int(length)
                if node["sequence_length"] < 0:
                    raise ValueError("Negative history sequence length")
            except ValueError as exc:
                node["parse_error"] = str(exc)
        nodes.append(node)
        for child in element:
            if child.tag == "Node":
                node["input_keys"].append(visit(child))
            else:
                node["properties"].append(_element(child))
                if child.tag == "InputSummary":
                    node["input_summaries"].append(_element(child))
        return key

    for child in root:
        if child.tag == "Node":
            tree["root_keys"].append(visit(child))
    return tree, nodes


def _snapshot(payload, summary, budget):
    entry = dict(
        _source(summary), id=None, packets=[], external_sequences={}, properties=[]
    )
    try:
        if len(payload) < 4:
            raise ValueError("Truncated history node ID")
        entry["id"] = str(struct.unpack_from(">I", payload)[0])
        for packet, content in _packets(payload[4:]):
            entry["packets"].append(packet)
            try:
                if packet["type"] == 29:
                    root = _xml(_decompress(content, budget))
                    if root.tag != "Manipulation":
                        raise ValueError("Expected Manipulation XML")
                    entry["manipulation"] = _element(root)
                    for direction in ("Undo", "Redo"):
                        entry[direction.lower()] = [
                            _element(a) for a in root.findall(direction + "/Action")
                        ]
                elif packet["type"] == 30:
                    properties = _decompress(content, budget)
                    for prop, value in _packets(properties):
                        item = dict(prop)
                        entry["properties"].append(item)
                        if value.lstrip().startswith(b"<"):
                            try:
                                item["xml"] = _element(_xml(value))
                            except (ValueError, ET.ParseError) as exc:
                                item["parse_error"] = str(exc)
                elif packet["type"] == 1:
                    if "sequence" in entry:
                        raise ValueError("Duplicate compact history sequence")
                    entry["sequence"], entry["sequence_flags"] = _compact_sequence(
                        content, budget
                    )
                elif packet["type"] == 31:
                    pos = 0
                    while pos < len(content):
                        if pos + 9 > len(content) or content[pos + 4] != 1:
                            raise ValueError(
                                "Unsupported external history sequence packet"
                            )
                        external_id, size = (
                            struct.unpack_from(">I", content, pos)[0],
                            struct.unpack_from(">I", content, pos + 5)[0],
                        )
                        end = pos + 9 + size
                        if end > len(content):
                            raise ValueError("Truncated external history sequence")
                        if str(external_id) in entry["external_sequences"]:
                            raise ValueError("Duplicate external history sequence ID")
                        sequence, _ = _compact_sequence(content[pos + 9 : end], budget)
                        entry["external_sequences"][str(external_id)] = sequence
                        pos = end
            except (ValueError, ET.ParseError, lzma.LZMAError, struct.error) as exc:
                packet["parse_error"] = str(exc)
    except (ValueError, struct.error) as exc:
        entry["parse_error"] = str(exc)
    return entry


def _apply_actions(sequence, actions, external):
    if any(base not in "ACGTRYMKSWBDHVNacgtrymkswbdhvn" for base in sequence):
        raise ValueError("Unsupported history sequence alphabet")
    for action in actions:
        attrs = action["attributes"]
        operation = attrs.get("type")
        if operation == "REMOVE":
            if set(attrs) != {"type", "range"} or action["children"]:
                raise ValueError("Unsupported history removal attributes")
            start, end = [int(value) for value in attrs["range"].split("-")]
            if not 0 <= start <= end < len(sequence):
                raise ValueError("Unsupported/out-of-bounds history removal range")
            sequence = sequence[:start] + sequence[end + 1 :]
        elif operation == "INSERT":
            if set(attrs) != {"type", "position"}:
                raise ValueError("Unsupported history insertion attributes")
            position = int(attrs["position"])
            residues = [
                child for child in action["children"] if child["tag"] == "Residues"
            ]
            if len(residues) != 1 or not 0 <= position <= len(sequence):
                raise ValueError("Invalid history insertion")
            attrs = residues[0]["attributes"]
            if attrs.get("type") == "GENERIC":
                if set(attrs) != {"type", "residues"}:
                    raise ValueError("Unsupported generic history residue attributes")
                bases = attrs["residues"]
                if any(base not in "ACGTRYMKSWBDHVNacgtrymkswbdhvn" for base in bases):
                    raise ValueError("Unsupported history insertion alphabet")
            elif attrs.get("type") == "EXTERNAL":
                if set(attrs) != {"type", "ID", "length"}:
                    raise ValueError("Unsupported external history residue attributes")
                bases = external[attrs["ID"]]
                if len(bases) != int(attrs["length"]):
                    raise ValueError("External history sequence length mismatch")
            else:
                raise ValueError("Unsupported history residue type")
            sequence = sequence[:position] + bases + sequence[position:]
        elif operation == "FLIP":
            if set(attrs) != {"type"} or action["children"]:
                raise ValueError("Unsupported partial history flip")
            sequence = sequence.translate(COMPLEMENT)[::-1]
        else:
            raise ValueError("Unsupported history action: %s" % operation)
        if len(sequence) > MAX_DECOMPRESSED_SIZE:
            raise ValueError("History sequence size limit exceeded")
    return sequence


def _recover(history, current_sequence, topology):
    by_key = {node["key"]: node for node in history["nodes"]}
    by_id = {}
    for node in history["nodes"]:
        by_id.setdefault(node["id"], []).append(node)
    snapshots = {}
    for index, snapshot in enumerate(history["snapshots"]):
        if snapshot["id"] is not None:
            snapshots.setdefault(snapshot["id"], []).append((index, snapshot))
    sequence_budget = [MAX_TOTAL_DECOMPRESSED_SIZE]

    def assign(node, sequence, source):
        if node["sequence_length"] is None or len(sequence) != node["sequence_length"]:
            raise ValueError("Recovered history sequence length mismatch")
        additional = len(sequence) - len(node.get("sequence", ""))
        if additional > sequence_budget[0]:
            raise ValueError("Recovered history sequence size limit exceeded")
        sequence_budget[0] -= additional
        node.update(
            sequence=sequence,
            recovery_status="current" if source == "document" else "recovered",
            sequence_source=source,
        )
        node.pop("recovery_reason", None)

    for node in history["nodes"]:
        if node["attributes"].get("type") != "DNA":
            node["recovery_reason"] = "Unsupported history molecule type"
            continue
        matches = snapshots.get(node["id"], [])
        if node["id"] is None or len(by_id[node["id"]]) != 1 or len(matches) > 1:
            node["recovery_reason"] = "Missing or ambiguous history node/snapshot ID"
            continue
        if matches:
            index, snapshot = matches[0]
            node["snapshot_index"] = index
            errors = [
                p["parse_error"] for p in snapshot["packets"] if "parse_error" in p
            ]
            if snapshot.get("parse_error"):
                errors.append(snapshot["parse_error"])
            if errors:
                node["recovery_reason"] = "; ".join(errors)
            if "sequence" in snapshot:
                try:
                    assign(node, snapshot["sequence"], "snapshot")
                except ValueError as exc:
                    node["recovery_reason"] = str(exc)

    roots = [key for tree in history["trees"] for key in tree["root_keys"]]
    if len(roots) == 1:
        root = by_key[roots[0]]
        if (
            root["attributes"].get("type") == "DNA"
            and root["sequence_length"] == len(current_sequence)
            and root["topology"] == topology
        ):
            try:
                if "sequence" in root and root["sequence"] != current_sequence:
                    raise ValueError("History root snapshot/document sequence mismatch")
                assign(root, current_sequence, "document")
            except ValueError as exc:
                root["recovery_warning"] = str(exc)

    # Each block 11 ID refers to an input node, whose Undo reverses its consumer's operation.
    queue = [node for node in history["nodes"] if "sequence" in node]
    for consumer in queue:
        if len(consumer["input_keys"]) != 1:
            continue
        node = by_key[consumer["input_keys"][0]]
        if "snapshot_index" not in node:
            continue
        snapshot = history["snapshots"][node["snapshot_index"]]
        undo = snapshot.get("undo", [])
        redo = snapshot.get("redo", [])
        if not undo and len(redo) == 1 and redo[0]["attributes"] == {"type": "FLIP"}:
            undo = redo
        if not undo:
            continue
        try:
            candidate = _apply_actions(
                consumer["sequence"], undo, snapshot["external_sequences"]
            )
            if (
                redo
                and _apply_actions(candidate, redo, snapshot["external_sequences"])
                != consumer["sequence"]
            ):
                raise ValueError("History Undo/Redo round trip mismatch")
            if "sequence" in node:
                if candidate != node["sequence"]:
                    raise ValueError("History snapshot/Undo sequence mismatch")
                continue
            assign(node, candidate, "undo")
            queue.append(node)
        except (ValueError, KeyError) as exc:
            node["recovery_warning" if "sequence" in node else "recovery_reason"] = str(
                exc
            )


def parse_history(blocks, current_sequence="", topology=None, current_document=None):
    """Return history trees, occurrence-keyed nodes, recovery data and errors.

    Unknown encodings remain inspectable and never replace the document sequence.
    """
    history = empty_history()
    budget = [MAX_TOTAL_DECOMPRESSED_SIZE]
    for summary, payload in blocks:
        try:
            if summary["type"] == 7:
                tree, nodes = _tree(
                    _decompress(payload, budget), summary, len(history["trees"])
                )
                history["trees"].append(tree)
                history["nodes"].extend(nodes)
            elif summary["type"] == 11:
                snapshot = _snapshot(payload, summary, budget)
                history["snapshots"].append(snapshot)
                if "parse_error" in snapshot:
                    history["parse_errors"].append(
                        dict(_source(summary), parse_error=snapshot["parse_error"])
                    )
                for packet in snapshot["packets"]:
                    if "parse_error" in packet:
                        history["parse_errors"].append(
                            dict(
                                _source(summary),
                                packet_type=packet["type"],
                                packet_offset=packet["offset"],
                                parse_error=packet["parse_error"],
                            )
                        )
        except (ValueError, ET.ParseError, lzma.LZMAError, struct.error) as exc:
            history["parse_errors"].append(dict(_source(summary), parse_error=str(exc)))
    _recover(history, current_sequence, topology)
    from .history_details import enrich_history

    enrich_history(history, current_document)
    return history
