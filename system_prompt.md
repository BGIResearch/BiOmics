# BiOmics Agent 系统 Prompt

## 一、项目概述

**BiOmics Agent** 是基于 LangGraph 构建的多智能体生物信息学分析系统，结合 BRICK 知识图谱，提供自动化单细胞/组学数据分析能力。

### 技术栈
- **前端**: NiceGUI (Python Web UI)
- **Agent 框架**: LangGraph
- **沙箱环境**: Docker + Jupyter Kernel Gateway (端口 8889)
- **知识图谱**: BRICK (Neo4j + Elasticsearch)
- **LLM**: 多种大语言模型

### 项目路径
```
/home/liyuntian/Biomics_agent/
```

---

## 二、Git 仓库配置

### 主项目仓库
| 远程名称 | URL | 说明 |
|----------|-----|------|
| **origin** | `git@github.com:lyt2003-yt/Biomics_agent.git` | 个人仓库 |
| **biomics** | `https://github.com/BGIResearch/BiOmics.git` | 组织仓库 |

### BRICK 子仓库
- 路径: `/home/liyuntian/Biomics_agent/BRICK/`
- 远程: `git@github.com:caolei2-BGI/BRICK.git`
- 分支: `Biomics_Agent` (用于 BiOmics Agent 定制版)

---

## 三、已完成的重要工作

### 3.1 基因标识符标准化
从知识图谱导出 NCBI Gene 节点，建立基因映射表：

| 文件 | 说明 |
|------|------|
| `ncbi_genes.csv` | 原始查询结果（有关系的 NCBI Gene 节点） |
| `ncbi_genes_filtered.csv` | 过滤 LOC 开头的基因 |
| `gene_mapping_standard.csv` | 标准化三列：ensembl_id, gene_symbol, ncbi_id |
| `gene_mapping_filled.csv` | 补充 ENSEMBL ID 后的最终版 |

### 3.2 ENSEMBL ID 补充方案
1. **mygene 库**: 通过 NCBI ID + Gene Symbol 在线查询
2. **gene2ensembl 本地映射**: NCBI 官方映射表（1800万+ 条）

### 3.3 补充结果
- 总记录: 73,816
- 有 ENSEMBL ID: 67,587 (91.6%)
- Unknown: 6,229 (8.4%) - 大部分是 QTL/tRNA/RIKEN克隆等无映射类型

---

## 四、基因标识符知识

### 常见类型
| 类型 | 示例 | 说明 |
|------|------|------|
| Gene Symbol | TP53, GAPDH | 标准基因符号 |
| Ensembl ID | ENSG00000141510 | Ensembl 数据库 ID |
| NCBI/Entrez ID | 7157 | NCBI Gene 数据库 |
| RefSeq ID | NM_000546 | NCBI 转录本 |

### 需过滤的预测性基因
| 前缀/后缀 | 物种 | 说明 |
|-----------|------|------|
| LOC | 多物种 | NCBI 临时预测基因 |
| Gm | 小鼠 | MGI 预测基因模型 |
| *Rik | 小鼠 | RIKEN cDNA 克隆 |
| C + 纯数字 | 多物种 | 未命名预测基因 |

### 正式命名的非编码 RNA（应保留）
- LINC (lncRNA)
- MIR (microRNA)
- SNORD (snoRNA)

---

## 五、关键文件位置

### 核心代码
```
graph/
├── builder.py          # LangGraph 工作流构建
├── state.py            # BrickState 状态定义
├── node2.py            # Agent 节点实现
├── llm.py              # LLM 配置
└── configuration.py    # 运行时配置
```

### 工具函数
```
tools/
├── read_table.py       # 表格文件读取
├── summarize_data.py   # 数据摘要
├── valid_h5ad.py       # h5ad 验证
└── ...

utils/
├── sandbox_client.py   # 沙箱客户端
├── gene_id_converter.py # 基因 ID 转换
└── ...
```

### 基因映射脚本
```
get_standard_gene.py        # 从 CSV 提取标准化映射表
fill_ensembl_id.py          # 用 mygene 补充 ENSEMBL ID
fill_from_gene2ensembl.py   # 用本地映射表补充
```

---

## 六、开发规范

### 代码规范
1. 节点间使用 `relations` 键名传递数据
2. 数据处理顺序：先获取 data_info
3. 图谱查询前需进行命名标准化
4. 跨模块传递路径需转为绝对路径

### 常见问题处理
1. `search_replace` 失败 → 先 `read_file` 获取最新内容
2. DataFrame 转字符串 → 用 `to_string()` 避免截断
3. AnnData var_names → 避免使用 `fillna`
4. h5ad 文件 → 需原地修改并 write_h5ad 持久化

---

## 七、Docker 沙箱启动

```bash
docker run -it --rm \
  --network=host \
  -v /home/liyuntian/Biomics_agent/data:/workspace/data \
  biomics_agent:v11 \
  jupyter kernelgateway \
    --ip=0.0.0.0 \
    --port=8889 \
    --KernelGatewayApp.allow_origin='*' \
    --JupyterWebsocketPersonality.list_kernels=True
```

---

## 八、Git 常用命令

```bash
# 推送到个人仓库
git add .
git commit -m "提交信息"
git push origin master

# 推送到组织仓库（需要权限）
git push biomics master

# BRICK 子仓库操作
cd BRICK
git checkout Biomics_Agent
git add .
git commit -m "提交信息"
git push origin Biomics_Agent
```

---

## 九、待处理/可扩展工作

1. 6,229 条 Unknown 基因大部分是 QTL/tRNA 等，无法补充 ENSEMBL ID
2. 可考虑将 gene_mapping_filled.csv 导入 Elasticsearch 用于快速查询
3. 基因 ID 转换功能可集成到 data_analyzer 流程中

---

*最后更新: 2026-01-20*
