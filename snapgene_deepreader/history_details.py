"""Expose saved node annotations, step oligos and branch-local visibility."""

from copy import deepcopy
from xml.etree import ElementTree as ET


def _xml_bytes(item):
    def build(record):
        element = ET.Element(record["tag"], record["attributes"])
        element.text = record.get("text")
        for child in record["children"]:
            element.append(build(child))
        return element

    return ET.tostring(build(item), encoding="utf-8")


def _annotations(properties, summary):
    # Imported at call time: the public reader also imports the history parser.
    from .snapgene_reader import (
        _parse_features_block,
        _parse_primers_block,
        parse_dict,
    )
    import xmltodict

    result = {"annotation_errors": []}
    for item in properties:
        tag = item["tag"]
        if tag not in ("Features", "Primers", "Notes"):
            continue
        try:
            payload = _xml_bytes(item)
            if tag == "Features":
                # A bad annotation is local to the feature, not the DNA node.
                features = result.setdefault("features", [])
                for feature in item["children"]:
                    if feature["tag"] != "Feature":
                        continue
                    container = {
                        "tag": "Features",
                        "attributes": {},
                        "children": [feature],
                    }
                    try:
                        _parse_features_block(_xml_bytes(container), result)
                    except (ValueError, KeyError, TypeError, IndexError) as exc:
                        features.append(
                            {
                                "name": feature["attributes"].get("name"),
                                "raw_xml": feature,
                                "parse_error": str(exc),
                            }
                        )
                        result["annotation_errors"].append(
                            {
                                "tag": tag,
                                "name": feature["attributes"].get("name"),
                                "parse_error": str(exc),
                            }
                        )
            elif tag == "Primers":
                result["primers"] = _parse_primers_block(payload, summary)
                result["primer_settings"] = [
                    child for child in item["children"] if child["tag"] != "Primer"
                ]
                result["annotation_errors"].extend(
                    {"tag": tag, "parse_error": p["parse_error"]}
                    for p in result["primers"]
                    if "parse_error" in p
                )
            else:
                result["notes"] = parse_dict(xmltodict.parse(payload)).get("Notes")
        except (ValueError, KeyError, TypeError, ET.ParseError) as exc:
            result["annotation_errors"].append({"tag": tag, "parse_error": str(exc)})
    return result


def _summary(source):
    return {
        "type": source["raw_block_type"],
        "size": source["raw_block_size"],
        "offset": source["raw_block_offset"],
    }


def _visibility(history):
    by_key = {node["key"]: node for node in history["nodes"]}

    def visit(key, hidden_by, uncertain_by):
        node = by_key[key]
        flag = node["attributes"].get("hideHistory", "0").lower()
        collapsed = {"0": False, "1": True, "false": False, "true": True}.get(flag)
        node.update(
            hide_history=collapsed,
            hidden_by=list(hidden_by),
            visibility_unknown_by=list(uncertain_by),
            is_visible=False if hidden_by else (None if uncertain_by else True),
        )
        if collapsed is None:
            node["visibility_error"] = "Unrecognized hideHistory value: " + flag
        for child in node["input_keys"]:
            visit(
                child,
                hidden_by + ([key] if collapsed is True else []),
                uncertain_by + ([key] if collapsed is None else []),
            )

    for tree_index, tree in enumerate(history["trees"]):
        for key in tree["root_keys"]:
            visit(key, [], [])
        members = [
            n for n in history["nodes"] if n["key"].split(":")[0] == str(tree_index)
        ]
        tree["visible_node_keys"] = [
            n["key"] for n in members if n["is_visible"] is True
        ]
        tree["hidden_node_keys"] = [
            n["key"] for n in members if n["is_visible"] is False
        ]
        tree["unknown_visibility_node_keys"] = [
            n["key"] for n in members if n["is_visible"] is None
        ]


def _used_primers(node, by_key):
    primers = []
    for item in node["properties"]:
        if item["tag"] != "Oligo":
            continue
        attrs = item["attributes"]
        primer = {
            "name": attrs.get("name"),
            "sequence": attrs.get("sequence"),
            "description": attrs.get("description", ""),
            "attributes": dict(attrs),
            "properties": item["children"],
            "source": "history_tree_oligo",
            "binding_sites": [],
        }
        for key in node["input_keys"]:
            template = by_key[key]
            # Tree primer sites are the saved operation-specific selection;
            # snapshot primers can also include unrelated/low-Tm sites.
            for source, candidates in (
                ("history_tree", template["tree_primers"]),
                (template["annotation_sources"].get("primers"), template["primers"]),
            ):
                matches = [
                    p
                    for p in candidates
                    if primer["sequence"]
                    and (p.get("sequence") or "").upper() == primer["sequence"].upper()
                ]
                named = [p for p in matches if p.get("name") == primer["name"]]
                matches = named or matches
                sites = []
                for match in matches:
                    for site in match.get("binding_sites", []):
                        sites.append(
                            dict(
                                deepcopy(site),
                                template_key=key,
                                source=source,
                                matched_primer_name=match.get("name"),
                                evidence="recorded_sequence_match",
                            )
                        )
                if sites:
                    primer["binding_sites"].extend(sites)
                    break
        primer["binding_sites_status"] = (
            "recorded" if primer["binding_sites"] else "unavailable"
        )
        primers.append(primer)
    return primers


def enrich_history(history, current_document=None):
    """Enrich only annotations belonging to a saved occurrence or matched snapshot.

    Missing annotations are not inherited from neighbouring nodes: their
    coordinates would refer to a different molecule. Oligos are recorded step
    reagents, not all primers attached to the product.
    """
    current_document = current_document or {}
    for snapshot in history["snapshots"]:
        properties = [p["xml"] for p in snapshot["properties"] if "xml" in p]
        snapshot.update(_annotations(properties, _summary(snapshot)))
    by_key = {node["key"]: node for node in history["nodes"]}
    for tree_index, tree in enumerate(history["trees"]):
        for node in history["nodes"]:
            if node["key"].split(":")[0] != str(tree_index):
                continue
            recorded = _annotations(node["properties"], _summary(tree))
            node["tree_features"] = recorded.get("features", [])
            node["tree_primers"] = recorded.get("primers", [])
            node["annotation_errors"] = list(recorded["annotation_errors"])
            node["annotation_sources"] = {}
            node["parameters"] = [
                dict(p["attributes"], properties=p["children"])
                for p in node["properties"]
                if p["tag"] == "Parameter"
            ]
            node["history_colors"] = [
                p for p in node["properties"] if p["tag"] == "HistoryColors"
            ]
            snapshot = (
                history["snapshots"][node["snapshot_index"]]
                if "snapshot_index" in node
                else {}
            )
            node["annotation_errors"].extend(snapshot.get("annotation_errors", []))
            for field in ("features", "primers", "notes", "primer_settings"):
                node[field] = [] if field != "notes" else {}
                for source, data in (
                    ("history_tree", recorded),
                    ("snapshot", snapshot),
                    (
                        "document",
                        (
                            current_document
                            if node.get("sequence_source") == "document"
                            else {}
                        ),
                    ),
                ):
                    if field in data:
                        node[field] = data[field]
                        node["annotation_sources"][field] = source
    for node in history["nodes"]:
        node["consumer_keys"] = []
    for node in history["nodes"]:
        node["inputs"] = []
        summaries = node["input_summaries"]
        for index, key in enumerate(node["input_keys"]):
            template = by_key[key]
            template["consumer_keys"].append(node["key"])
            node["inputs"].append(
                {
                    "key": key,
                    "id": template["id"],
                    "name": template["name"],
                    "summary": (
                        summaries[index]
                        if len(summaries) == len(node["input_keys"])
                        else None
                    ),
                    "snapshot_index": template.get("snapshot_index"),
                }
            )
        node["used_primers"] = _used_primers(node, by_key)
    _visibility(history)
