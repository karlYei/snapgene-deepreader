# snapgene-deepreader

[English](README.md) | [中文](README.zh-CN.md)

`snapgene-deepreader` 是一个增强版 SnapGene `.dna` 文件解析器，基于
MIT 许可证开源项目
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)
开发。它保留原项目的高级 API，同时新增对 SnapGene 文件中引物 block 和
对比序列 block 的解析。

## 项目来源与开发说明

本项目来源于 Edinburgh Genome Foundry 维护的开源项目
[SnapGeneReader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)。
原项目采用 MIT License。

将其独立为 `snapgene-deepreader` 包的全部改造均由 AI coding agent 完成，
包括包名重命名、引物 block 解析、对比序列 metadata 解析、基于 zlib 的
对比序列内容提取、测试以及项目文档。

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
- `raw_blocks_summary`：每个 SnapGene block 的类型、大小和文件偏移，用于调试未知 block。

## 兼容性

本项目提供与 `snapgene-reader` 相同调用方式的函数：

- `snapgene_file_to_dict(filepath=None, fileobject=None)`
- `snapgene_file_to_seqrecord(filepath=None, fileobject=None)`
- `snapgene_file_to_gbk(read_file_object, write_file_object)`

未知或暂不支持的 rich block 不会导致整个 `.dna` 文件读取失败。解析失败信息会记录在对应结果对象的 `parse_error` 字段中。
