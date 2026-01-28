# BiOmics Agent 项目完整介绍

## 一、项目概述

**BiOmics Agent** 是一个基于 **LangGraph** 构建的多智能体生物信息学分析系统，结合 **BRICK 知识图谱**，为用户提供自动化的单细胞/组学数据分析能力。系统通过 NiceGUI 提供 Web 界面，代码在 Docker 沙箱（Jupyter Kernel Gateway）中隔离执行。

### 核心能力
- **自动化分析**: 用户输入自然语言问题，系统自动规划分析步骤、生成代码、执行并解释结果
- **知识图谱增强**: 通过 BRICK 生物医学知识图谱提供实体查询、关系分析、富集解释
- **沙箱隔离**: 所有代码在 Docker 容器中执行，确保安全性
- **Notebook 导出**: 分析过程可导出为可复现的 Jupyter Notebook

## 二、技术架构

### 2.1 核心技术栈
- **前端框架**: NiceGUI (Python Web UI)
- **Agent 框架**: LangGraph (状态图驱动的多 Agent 工作流)
- **沙箱环境**: Docker + Jupyter Kernel Gateway (端口 8888)
- **知识图谱**: BRICK (Neo4j + Elasticsearch 混合检索)
- **向量检索**: TF-IDF + 余弦相似度 (Notebook 匹配)
- **LLM**: 支持多种大语言模型调用

### 2.2 Agent 工作流节点
```
START → notebook_searcher → supervisor
          ↓
supervisor → env_checker / analyze_planner / general_responder
          ↓                    ↓
env_checker              (名词查询类问题直接走 analyze_planner)
          ↓
data_analyzer → analyze_planner
          ↓
analyze_planner → planner → plan_executor
          ↓
plan_executor → coder → code_runner
          ↓ (成功)           ↓ (失败)
plan_executor          code_debugger → code_runner
          ↓
responder / general_responder → END
```

### 2.3 关键文件结构
```
Biomics_agent/
├── app_nicegui.py          # Web UI 入口 (NiceGUI)
├── run.py                  # Agent 运行入口
├── graph/
│   ├── builder.py          # LangGraph 工作流构建
│   ├── state.py            # BrickState 状态定义 (50+ 字段)
│   ├── node2.py            # 各 Agent 节点实现
│   ├── notebook_finder.py  # 基于 RAG 的 Notebook 匹配器 (中英文支持)
│   ├── llm.py              # LLM 配置
│   └── configuration.py    # 运行时配置
├── notebooks/              # 分析模板 Notebook (12个)
├── prompt/                 # 各 Agent 的 Prompt 模板
│   ├── supervisor.md       # 路由决策 Prompt
│   ├── coder.md            # 代码生成 Prompt
│   ├── code_debugger.md    # 代码调试 Prompt
│   └── ...                 # 其他 Agent Prompt
├── utils/
│   ├── sandbox_manager.py  # Docker 沙箱管理
│   ├── sandbox_client.py   # WebSocket 沙箱客户端
│   ├── create_notebook.py  # 生成 ipynb 文件工具
│   └── ...
├── tools/                  # 工具函数
│   ├── read_table.py       # 表格文件转二维列表
│   ├── summarize_data.py   # 数据摘要
│   ├── valid_h5ad.py       # h5ad 文件验证
│   └── ...
├── BRICK/                  # BRICK 知识图谱工具包
│   ├── interpretation.py   # LLM 解释与实体抽取
│   ├── SearchClient.py     # Elasticsearch 混合检索客户端
│   ├── querygraph.py       # 图谱查询 (qr)
│   ├── rankgraph.py        # 结果排序 (rk)
│   ├── embedgraph.py       # 图表示学习 (emb)
│   ├── entity_index/       # 实体索引模块
│   └── entity_vocab/       # 实体词表数据 (~72MB)
└── data/uploaded_file/     # 用户上传文件存储 (按时间戳子目录组织)
```

## 三、BrickState 核心状态字段

| 字段 | 类型 | 说明 |
|------|------|------|
| `question` | str | 用户输入的问题 |
| `data_path` | str | 上传数据路径 |
| `data_info` | dict | 数据概览信息 |
| `current_plan` | dict/list | 当前执行计划 |
| `current_step` | int | 当前步骤编号 |
| `code` | str | 当前生成的代码 |
| `complete_output` | dict | 代码执行结果 |
| `notebook_cells` | list | 用于生成 ipynb 的单元格记录 |
| `next` | str | 下一个 Agent 节点 |
| `status` | str | 任务状态 (NOT_FINISHED/FINISHED) |
| `sandbox_id` | str | Docker 沙箱 ID |

## 四、支持的分析类型 (Notebooks)

### 4.1 组学数据分析
| Notebook | 功能 | 适用场景 |
|----------|------|----------|
| Celltype_annotation | 细胞类型注释 | 单细胞数据自动注释 |
| Celltype_refinement | 细胞类型精化 | 细胞亚群细分 |
| Differential_Expression_Gene | 差异表达基因 | 组间差异分析 |
| Drug_Discovery | 药物发现 | 疾病相关药物预测 |
| Enrichment | 基因富集分析 | GO/Pathway 富集 |
| GeneRegulatoryNetwork | 基因调控网络 | TF-Regulon 分析 |
| GWAS_Causal_SNPs | GWAS 因果 SNP | 疾病关联位点分析 |
| GWAS_Predict_Phenotype | GWAS 表型预测 | 多基因风险评分 |
| Proteome | 蛋白质组学 | 蛋白表达分析 |
| Trajectory | 细胞轨迹推断 | 发育轨迹/分化分析 |

### 4.2 知识图谱查询 (名词查询)
| Notebook | 功能 | 典型问题 |
|----------|------|----------|
| Terminology_query | 术语定义查询 | "TP53是什么?" "二型糖尿病是什么?" |
| Relation_query | 实体关系查询 | "TP53与癌症有什么关系?" |
| Most_related_query | 最相关实体 | "与阿尔茨海默病最相关的基因是什么?" |

**名词类型**: 基因、疾病、物种、细胞、化学物质、蛋白质、通路等生物医学实体

## 五、UI 功能 (app_nicegui.py)

### 5.1 界面组成
- **顶栏**: Logo + 标题 + 重置按钮
- **主区域**: 左侧 Copilot 面板 (数据信息) + 右侧 Chat 面板 (对话交互)
- **底栏**: 文件上传按钮 + 输入框
- **浮动按钮栏**: 8 个演示按钮 (页面加载 3 秒后渐现动画)

### 5.2 演示按钮配置
```python
DEMO_BUTTON_CONFIG = {
    'Demonstrate cell annotation': ('/path/adata_new1.h5ad', 'Perform cell type annotation...'),
    'Demonstrate drug discovery': ('/path/Neutrophil_adata_sub.h5ad', 'Predict therapeutic drugs for COVID-19...'),
    # ... 共 8 个预设任务
}
```

### 5.3 交互逻辑
- **用户输入问题 → 隐藏浮动按钮栏 → 启动 Agent 工作流**
- **点击重置 → 清空聊天/关闭沙箱/重新显示按钮栏**
- **Human-in-the-loop**: `env_checker` 和 `data_analyzer` 支持中断等待用户确认

## 六、Notebook 生成功能 (create_notebook.py)

将会话中的代码执行记录导出为 `.ipynb` 文件：
```python
cells_data = [
    {
        "mdata": "Step1: Data Loading",     # → Markdown 标题
        "thought": "Coder 的思考过程",       # → Markdown 内容
        "code": "import scanpy as sc\n...", # → Code Cell
        "outputs": {"stdout": [...], "images": [...]}  # → Cell Outputs
    },
    ...
]
create_notebook(cells_data, "output.ipynb")
```

## 七、Notebook 匹配器 (notebook_finder.py)

基于 TF-IDF + 余弦相似度的 RAG 系统，支持中英文双语匹配：
- 每个 Notebook 预定义 `keywords`, `functions`, `synonyms`, `scenarios`, `common_questions`, `technical_terms`
- 使用 jieba 中文分词
- 缓存向量化结果到 `notebook_cache/`

## 八、环境配置要点

1. **环境变量**: `graph/brick_test_config.env` 配置 `PROJECT_ROOT`
2. **Docker 沙箱**: 端口 8888，挂载 `data/` 目录
3. **LangGraph 递归限制**: `recursion_limit: 200`
4. **文件上传限制**: 100MB

## 九、已完成的开发工作

1. ✅ 完整的多 Agent 工作流 (15+ 节点)
2. ✅ NiceGUI Web 界面 (响应式布局)
3. ✅ 8 个演示按钮 + 点击自动执行任务
4. ✅ 浮动按钮栏动画 (渐现效果)
5. ✅ Notebook 匹配器 (中英文双语)
6. ✅ 代码执行结果导出为 ipynb
7. ✅ Human-in-the-loop 交互 (env_checker/data_analyzer)
8. ✅ 沙箱隔离执行 + 自动清理
