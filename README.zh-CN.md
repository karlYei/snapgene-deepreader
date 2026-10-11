# snapgene-deepreader

[English](README.md) | [中文](README.zh-CN.md)

`snapgene-deepreader` 是一个增强版 SnapGene `.dna` 文件解析器，基于
MIT 许可证开源项目
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)
开发。它保留原项目的高级 API，同时新增对 SnapGene 文件中引物 block 和
对比序列 block、SnapGene 保存的对比结果数据及构建历史的解析。

## 项目来源与开发说明

本项目来源于 Edinburgh Genome Foundry 维护的开源项目
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)。
原项目采用 MIT License。

将其独立为 `snapgene-deepreader` 包的全部改造均由 AI coding agent 完成，
包括包名重命名、引物 block 解析、对比序列 metadata 解析、基于 zlib 的
对比序列内容提取、基于 BAM 的对比结果解析、历史解析与序列恢复、测试以及项目文档。

## 本地开发安装

```bash
python -m pip install -e .
```

## 使用方法

```python
from snapgene_deepreader import snapgene_file_to_dict

data = snapgene_file_to_dict("example.dna")
print(data["seq"])
print(data["primers"])
print(data["alignable_sequences"])
print(data["comparison_alignments"])
print(data["history"])
```

返回的字典保留 `snapgene-reader` 原有字段：

- `dna`
- `seq`
- `features`
- `notes`
- `isDNA`
- `exportVersion`
- `importVersion`

同时新增：

- `primers`：来自 SnapGene block `5` 的引物信息。
- `alignable_sequences`：来自 SnapGene blocks `16` 和 `17` 的对比序列信息。
- `comparison_alignments`：来自 SnapGene block `27` 的对比结果信息；当该 block
  包含 BGZF/GZIP 压缩的 BAM alignment 时，会解析 BAM record、CIGAR、read
  sequence、aux tags，并根据 SnapGene 参考序列还原 `variants` 列表。
- `raw_blocks_summary`：每个 SnapGene block 的类型、大小和文件偏移，用于调试未知 block。

## 读取构建历史

0.3.0 版本新增 `data["history"]`，读取 block `7` 中保存的 HistoryTree
和 block `11` 中的节点恢复数据。即使文件没有历史，也会返回结构一致的对象，
以下四个字段为空列表：

- `trees`：历史树属性、根节点 `root_keys`、原始 XML 和源 block 位置。
- `nodes`：每次节点出现对应一个对象，包含 `key`、原始 `id`、`name`、
  `operation`、`sequence_length`、`topology`、`strandedness`、`input_keys`、
  `input_summaries`，以及保留原始 XML 内容的 `attributes` 和 `properties`。
- `snapshots`：节点 ID、Undo/Redo 操作、已解码的外部片段、当时保存的 XML
  属性（包括存在时的 features）和内部 packet 摘要。
- `parse_errors`：解析失败的 block/packet 位置和错误信息。
  snapshot 和 packet 对象还会在本地记录 `parse_error`。

`input_keys` 表示某个产物的输入来源。多输入和重复 ID 均予以保留；遍历时使用
每次出现唯一的 `key`，不要用 ID 或列表下标推断关系。节点列表按树遍历顺序排列，
不代表时间顺序。文件保存的日期、作者等属性会保留，没有保存的不会推断补齐。

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

只有能得到序列时才返回 `sequence`。`recovery_status` 为 `current`（当前文件序列）、
`recovered`（已恢复）或 `unavailable`（无法可靠恢复）；`sequence_source` 标明来自
`document`、`snapshot` 或 `undo`。不可恢复节点带有 `recovery_reason`，
校验结果存在冲突时记录 `recovery_warning`。存在 `snapshot_index` 时，
可以用它定位 `snapshots` 中的恢复数据。

目前支持已验证的二位压缩 DNA run 编码 `1`、大小写区间和外部片段，
以及保存的 INSERT/REMOVE 操作（包括先删后插的替换）和整条序列的 FLIP。
恢复长度需与历史节点一致；有 Redo 时还会用它校验 Undo 的结果。
历史操作坐标**从 0 开始，删除区间包含两端**，输出保留文件原值。

此版本覆盖已验证的格式，不代表支持全部 SnapGene 版本。
其他压缩编码、跨环状序列原点的删除、歧义 ID 和未知操作会保留诊断信息，
而不返回猜测序列。多输入构建关系会保留，但不会自动反推其所有输入序列。
历史 XML 属性与恢复序列分别保存，不承诺生成完整的历史 SeqRecord。
解析设有解压大小、XML 深度和恢复序列总大小限制，不新增运行时依赖。

原 SeqRecord 和 GenBank API 仍导出当前文件内容。
这里读取的是文件实际保存的构建历史，不是每次编辑或保存的完整操作日志。

## 兼容性

本项目提供与 `snapgene-reader` 相同调用方式的函数：

- `snapgene_file_to_dict(filepath=None, fileobject=None)`
- `snapgene_file_to_seqrecord(filepath=None, fileobject=None)`
- `snapgene_file_to_gbk(read_file_object, write_file_object)`

未知或暂不支持的 rich block 不会导致整个 `.dna` 文件读取失败。解析失败信息会记录在对应结果对象的 `parse_error` 字段中。
