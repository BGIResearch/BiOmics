# -*- coding: utf-8 -*-
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
import copy
import pandas as pd

# ===============================
# 封装为类：BRICKSearchClient
# 外部可直接构造 payload 并进行模糊/混合检索
# ===============================
class BRICKSearchClient:
    def __init__(self,
                 project_root: Path | None = None,
                 env_path: Path | None = Path(r'EMBEDDING_API_KEY'),
                 gene_mapping_file: Path | None = None):
        # 1) 将仓库根目录加入 sys.path，方便直接 import entity_index 包
        if project_root is None:
            # 默认指向 SearchClient.py 所在的目录（即 BRICK 目录）
            project_root = Path(__file__).parent.resolve()
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))
        
        # 保存 project_root 用于后续加载基因映射表
        self.project_root = project_root

        # 2) 读取 .env（若存在）
        # if env_path and Path(env_path).exists():
        #     for raw_line in Path(env_path).read_text(encoding="utf-8").splitlines():
        #         line = raw_line.strip()
        #         if not line or line.startswith("#") or "=" not in line:
        #             continue
        #         key, value = line.split("=", 1)
        #         os.environ.setdefault(key, value)
        # 2) 写入环境变量（直接内嵌 .env 内容）
        os.environ.update({
            "EMBEDDING_API_KEY": "sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",
            "EMBEDDING_BASE_URL": "http://10.224.28.80:3000/v1",
            "EMBEDDING_MODEL": "text-embedding-v3",

            "ES_CONFIG": json.dumps({
                "hosts": ["http://es-cn-uax3fxy8w000e2tuf.elasticsearch.aliyuncs.com:9200"],
                "username": "elastic",
                "password": "3fg51u8ZyRYWUxPXwVpl"
            }),
            "ES_INDEX_PREFIX": "brick_index_*",
            "ES_VERIFY_CERTS": "true",

            "HYBRID_ES_INDEX": "brick_index_*",  # ES 会自动拼接 _string / _vector
            "HYBRID_TOP_K": "10",
            "HYBRID_DISABLE_VECTOR": "false",  # 明确启用向量检索
            "HYBRID_ALPHA": "0.6",  # 字符串权重
            "HYBRID_BETA": "0.4",   # 向量权重
            "HYBRID_TYPE_MIX": json.dumps({
                "Gene|Protein": 0.25,
                "Disease|Phenotype": 0.25,
                "Process|Function|Pathway|Cell_Component": 0.20,
                "Chemical": 0.10,
                "Species": 0.10,
                "Cell|Tissue": 0.10,
                "Mutation": 0.0
            }),
            "HYBRID_EMBEDDING_ENDPOINT": "http://10.224.28.80:3000/v1",
            "HYBRID_EMBEDDING_MODEL": "text-embedding-v3",
        })

        # 3) Notebook 中演示时若未显式设置类型权重，提供一个默认 JSON 值
        os.environ.setdefault("HYBRID_TYPE_MIX", json.dumps({
            "Gene|Protein": 0.25,
            "Disease|Phenotype": 0.25,
            "Process|Function|Pathway|Cell_Component": 0.20,
            "Chemical": 0.10,
            "Species": 0.10,
            "Cell|Tissue": 0.10,
            "Mutation": 0.0,
        }))

        # 4) 确保向量检索启用（避免被其他配置覆盖）
        os.environ["HYBRID_DISABLE_VECTOR"] = "false"
        os.environ["HYBRID_ES_INDEX"] = "brick_index_*"  # *_string / *_vector

        # 5) 依赖导入（放在 sys.path 注入和环境准备之后）
        from entity_index.search.settings import get_search_config
        from entity_index.search.hybrid_searcher import HybridEntitySearcher
        from entity_index.search.schema import HYBRID_TYPE_KEYS

        # 6) 构建配置与 ES 客户端
        self.search_config = get_search_config()
        self.es_client = self.search_config.es.create_client()
        self.hybrid_searcher = HybridEntitySearcher(self.es_client, self.search_config)
        self.TYPE_KEYS = tuple(HYBRID_TYPE_KEYS)
        
        # 7) 加载基因 ID 映射表（新增）
        if gene_mapping_file is None:
            # 映射表位于 BRICK 目录下
            gene_mapping_file = self.project_root / "gene_mapping.csv"
        self._load_gene_mapping(gene_mapping_file)

    # ---- 工具方法 ----
    def ensure_type_keys(self, payload: dict) -> dict:
        for k in self.TYPE_KEYS:
            payload.setdefault(k, [])
        return payload
    
    def _load_gene_mapping(self, mapping_file: Path):
        """加载基因 ID 映射表并构建索引"""
        try:
            if not Path(mapping_file).exists():
                print(f"[GeneMapping] Warning: Mapping file not found: {mapping_file}")
                self.ensembl_to_symbol = {}
                self.ncbi_to_symbol = {}
                return
            
            df = pd.read_csv(mapping_file, dtype=str)
            df = df.replace("Unknown", pd.NA)
            
            # 构建快速查询字典：Ensembl ID -> Gene Symbol
            valid_ensembl = df[df["ensembl_id"].notna() & (df["ensembl_id"] != "Unknown")]
            self.ensembl_to_symbol = dict(zip(
                valid_ensembl["ensembl_id"], 
                valid_ensembl["gene_symbol"]
            ))
            
            # 构建快速查询字典：NCBI ID -> Gene Symbol
            valid_ncbi = df[df["ncbi_id"].notna()]
            self.ncbi_to_symbol = dict(zip(
                valid_ncbi["ncbi_id"], 
                valid_ncbi["gene_symbol"]
            ))
    
        except Exception as e:
            print(f"[GeneMapping] Warning: Failed to load gene mapping: {e}")
            self.ensembl_to_symbol = {}
            self.ncbi_to_symbol = {}
    
    def normalize_gene_query(self, query: str) -> str:
        """标准化基因查询词：自动将 Ensembl/NCBI ID 转为 gene symbol
        
        Args:
            query: 原始查询词，可能是 gene symbol / ensembl id / ncbi id
        
        Returns:
            标准化后的 gene symbol（如果能转换），否则返回原值
        """
        query = query.strip()
        
        # 1. 判断是否是 Ensembl ID（以 ENS 开头）
        if query.startswith("ENS"):
            symbol = self.ensembl_to_symbol.get(query)
            if symbol:
                return symbol
        
        # 2. 判断是否是 NCBI ID
        if query.startswith("NCBI:") or query.isdigit():
            ncbi_key = query if query.startswith("NCBI:") else f"NCBI:{query}"
            symbol = self.ncbi_to_symbol.get(ncbi_key)
            if symbol:
                return symbol
        
        # 3. 无法转换，返回原值
        return query
    
    def normalize_gene_queries(self, queries: list[str]) -> list[str]:
        """批量标准化基因查询词"""
        return [self.normalize_gene_query(q) for q in queries]


    def search_hybrid(self,
                    payload: dict,
                    *,
                    top_k: int | None = None,
                    type_mix_override: dict | None = None,
                    return_diagnostics: bool | None = None,
                    debug: bool | None = None,
                    auto_normalize_gene: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, object]:
        """
        多类别 Hybrid 检索逻辑：
        - 若 payload 指定了若干类别（相应 type_key 列表非空），仅检索这些类别；
        - 否则从 options.query_text 读取查询词，对所有类别逐一检索；
        - 每个类别独立调用 hybrid_searcher.search()；若异常则回退到字符串通道；
        - 合并所有类别的诊断结果，按 final_score 全局排序，(entity_id, type_key) 去重取最高分；
        - 返回：standardized_df（按分数对齐后的 top_k）、diagnostics_df（全局 top_k）、summary 对象。
        
        Args:
            auto_normalize_gene: 是否自动将 Gene|Protein 类别的 Ensembl/NCBI ID 转为 gene symbol（默认 True）
        """
        from entity_index.search import adapters, string_client

        # ---- 统一参数与开关 ----
        payload_h = self.ensure_type_keys(copy.deepcopy(payload))
        
        # ---- 新增：对 Gene|Protein 类别自动标准化 ----
        if auto_normalize_gene and "Gene|Protein" in payload_h and payload_h["Gene|Protein"]:
            payload_h["Gene|Protein"] = self.normalize_gene_queries(payload_h["Gene|Protein"])
        opts = payload_h.setdefault("options", {})
        if top_k is not None:
            opts["top_k"] = int(top_k)
        if return_diagnostics is not None:
            opts["return_diagnostics"] = bool(return_diagnostics)
        if debug is not None:
            opts["debug"] = bool(debug)
        if type_mix_override:
            opts["type_mix_override"] = dict(type_mix_override)

        top_k_eff = int(opts.get("top_k", self.search_config.top_k))

        # ---- 判定是否“指定了类别” ----
        specified = {tk: payload_h.get(tk, []) for tk in self.TYPE_KEYS if payload_h.get(tk)}
        tasks: list[tuple[str, list[str]]] = []
        if specified:
            for tk, qlist in specified.items():
                if qlist:
                    tasks.append((tk, list(qlist)))
        else:
            qtext = opts.get("query_text")
            if not qtext:
                raise ValueError(
                    "未指定类别且缺少 options.query_text："
                    "请在 payload['options']['query_text'] 填写查询词，或在某些类型键下给出查询列表。"
                )
            for tk in self.TYPE_KEYS:
                tasks.append((tk, [str(qtext)]))

        # ---- 逐类别执行检索并合并 ----
        diag_rows = []
        std_rows = []
        per_type_summaries = []

        for tk, queries in tasks:
            per_payload = {
                "query_id": payload_h.get("query_id", "hybrid-multi"),
                "options": opts,     # 共享 options（包含 top_k、type_mix_override 等）
                tk: list(queries),   # 仅该类别有查询词，其它类别缺省/为空
            }

            try:
                # per_norm = adapters.normalize_payload(per_payload, self.search_config)
                # resp = self.hybrid_searcher.search(per_norm)
                # 让 HybridEntitySearcher 自己做 normalize（这里不要传 NormalizedQuery）
                resp = self.hybrid_searcher.search(per_payload)
                # standardized
                for type_key, entities in (resp.standardized or {}).items():
                    for entity_name in entities:
                        std_rows.append({"type_key": type_key, "entity_name": entity_name})

                # diagnostics
                for item in (resp.diagnostics or []):
                    diag_rows.append({
                        "type_key": item.type_key,
                        "entity_id": item.entity_id,
                        "primary_name": item.primary_name,
                        "final_score": item.final_score,
                        "string_score": getattr(item.channel_scores, "string_score", 0.0),
                        "vector_score": getattr(item.channel_scores, "vector_score", 0.0),
                        "matched_alias": item.matched_alias,
                    })

                per_type_summaries.append({"type_key": tk, "mode": "hybrid_ok"})

            except Exception as e:
                import traceback
                print(traceback.format_exc())
                print(e)
                # ---- 单类别回退为字符串通道 ----
                print('回退为字符串通道')
                per_norm = adapters.normalize_payload(per_payload, self.search_config)
                string_hits = string_client.search_string_channel(
                    self.es_client,
                    per_norm,
                    self.search_config.string_index_name,
                )
                for h in string_hits:
                    # standardized 用 primary_name 作为实体名
                    std_rows.append({"type_key": h.type_key, "entity_name": h.primary_name})
                    # diagnostics：final=string, vector=0
                    diag_rows.append({
                        "type_key": h.type_key,
                        "entity_id": h.entity_id,
                        "primary_name": h.primary_name,
                        "final_score": getattr(h.scores, "string_score", 0.0),
                        "string_score": getattr(h.scores, "string_score", 0.0),
                        "vector_score": 0.0,
                        "matched_alias": h.matched_alias,
                    })
                per_type_summaries.append({"type_key": tk, "mode": "string_only_fallback"})

        # ---- 汇总 DataFrame：按分数全局排序并去重 ----
        if not diag_rows:
            # 没有命中
            empty_std = pd.DataFrame(columns=["type_key", "entity_name"])
            empty_diag = pd.DataFrame(columns=[
                "type_key", "entity_id", "primary_name", "final_score",
                "string_score", "vector_score", "matched_alias"
            ])
            return empty_std, empty_diag, {"mode": "multi-type", "types": [t for t, _ in tasks]}

        diagnostics_df = pd.DataFrame(diag_rows)
        diagnostics_df = diagnostics_df.sort_values("final_score", ascending=False)
        # 同一 (entity_id, type_key) 保留最高分
        diagnostics_df = diagnostics_df.drop_duplicates(subset=["entity_id", "type_key"], keep="first")
        diagnostics_df = diagnostics_df.head(top_k_eff)

        # standardized 与诊断对齐（按 primary_name/type_key 关联拿到 top_k）
        standardized_df = pd.DataFrame(std_rows).drop_duplicates()
        standardized_df = standardized_df.merge(
            diagnostics_df[["type_key", "primary_name", "final_score"]],
            left_on=["type_key", "entity_name"],
            right_on=["type_key", "primary_name"],
            how="inner"
        ).rename(columns={"final_score": "score"}).drop(columns=["primary_name"])
        standardized_df = standardized_df.sort_values("score", ascending=False).head(top_k_eff)

        summary = {"mode": "multi-type", "types": [t for t, _ in tasks], "per_type": per_type_summaries}
        return standardized_df, diagnostics_df, summary


    def build_untyped_payload(self, query_text: str, top_k: int = 10):
        """
        构造“未指定类别”的通用 payload：
        - 所有类型键置空 []
        - 查询词放到 options.query_text
        - search_fuzzy 若发现所有类型都为空，就会对 self.TYPE_KEYS 逐类检索并合并
        """
        payload = {
            "query_id": "demo-generic-001",
            "options": {
                "top_k": int(top_k),
                "return_diagnostics": True,
                "query_text": str(query_text),
            },
        }
        # 显式置空所有类型键，表示“未指定类别”
        for type_key in self.TYPE_KEYS:
            payload[type_key] = []
        return payload


    def search_fuzzy(self, payload: dict, auto_normalize_gene: bool = True) -> pd.DataFrame:
        """
        字符串通道检索（可多类别）：
        - 若 payload 中存在"非空"的类型键，则仅检索这些已指定的类别（每类独立搜）
        - 若全部类型键为空，则读取 options.query_text，对 self.TYPE_KEYS 全部类别逐一检索
        - 合并所有类别的命中，按 entity_id+type_key 去重保留最高 string_score，再全局排序取 top_k
            
        Args:
            auto_normalize_gene: 是否自动将 Gene|Protein 类别的 Ensembl/NCBI ID 转为 gene symbol（默认 True）
        """
        from entity_index.search import adapters, string_client
    
        # 让 payload 至少包含所有类型键（空列表）
        payload = self.ensure_type_keys(dict(payload))
            
        # ---- 新增：对 Gene|Protein 类别自动标准化 ----
        if auto_normalize_gene and "Gene|Protein" in payload and payload["Gene|Protein"]:
            payload["Gene|Protein"] = self.normalize_gene_queries(payload["Gene|Protein"])
        opts = payload.setdefault("options", {})
        top_k = int(opts.get("top_k", self.search_config.top_k))

        # === 关键修复：用“原始 payload”来判断哪些类别被指定了 ===
        # 指定类别：该 type_key 的列表非空
        specified = {tk: payload.get(tk, []) for tk in self.TYPE_KEYS if payload.get(tk)}

        # 准备任务列表
        tasks: list[tuple[str, list[str]]] = []
        if specified:
            # 只查这些被显式指定的类别
            for tk, qlist in specified.items():
                if qlist:
                    tasks.append((tk, list(qlist)))
        else:
            # 未指定类别：从 options.query_text 读取查询词，并对所有类别逐一检索
            qtext = opts.get("query_text")
            if not qtext:
                raise ValueError(
                    "未指定类别且缺少 options.query_text："
                    "请在 payload['options']['query_text'] 填写查询词，或在某些类型键下给出查询列表。"
                )
            for tk in self.TYPE_KEYS:
                tasks.append((tk, [str(qtext)]))

        # 逐类别检索并收集结果
        rows = []
        for tk, queries in tasks:
            per_payload = {
                "query_id": payload.get("query_id", "string-multi"),
                "options": opts,
                tk: list(queries),
            }
            per_normalized = adapters.normalize_payload(per_payload, self.search_config)
            string_hits = string_client.search_string_channel(
                self.es_client,
                per_normalized,
                self.search_config.string_index_name,
            )
            for hit in string_hits:
                rows.append({
                    "entity_id": hit.entity_id,
                    "primary_name": hit.primary_name,
                    "type_key": hit.type_key,
                    "node_type": hit.node_type,
                    "string_score": getattr(hit.scores, "string_score", 0.0),
                    "matched_alias": hit.matched_alias,
                })

        if not rows:
            return pd.DataFrame(columns=[
                "entity_id", "primary_name", "type_key", "node_type", "string_score", "matched_alias"
            ])

        df = pd.DataFrame(rows)
        df = df.sort_values("string_score", ascending=False)
        df = df.drop_duplicates(subset=["entity_id", "type_key"], keep="first")
        return df.head(top_k)


if __name__ == "__main__":
    # ===============================
    # 配置查询词（全局参数，只改这里）
    # ===============================
    query_text = "isl1"  # 示例，可改成任何英文词，如 "lung cancer"
    # query_cells =[
    #     "excitatory (EX)", 
    #     "inhibitory (INH)", 
    #     "astrocytic (ASC)", 
    #     "microglial (MG)", 
    #     "oligodendrocytic (ODC)", 
    #     "oligodendrocyte precursor (OPC)"
    # ]

    # ===============================
    # 初始化客户端与环境信息
    # ===============================
    client = BRICKSearchClient()
    print(f"字符串索引: {client.search_config.string_index_name}")
    print(f"向量索引: {client.search_config.vector_index_name or '已禁用'}\n")

    # ===============================
    # A) build_untyped_payload：未指定类别 → 全类别模糊检索
    # ===============================
    # payload_fuzzy_auto = client.build_untyped_payload(query_text=query_text, top_k=1)
    # print(f">>> [Fuzzy] 未指定类别：全类别检索（query_text='{query_text}'）")
    # df_all = client.search_fuzzy(payload_fuzzy_auto)
    # entity_id = df_all.head(1)['entity_id'].iloc[0]
    # print(entity_id)
    # print(df_all.head(1)['entity_id'])
    # print(f"命中类别数：{df_all['type_key'].nunique()}\n")

    # ===============================
    # B) build_untyped_payload + 指定类别 → 只查该类
    # ===============================
    # payload_fuzzy_gene = client.build_untyped_payload(query_text=query_text, top_k=10)
    # payload_fuzzy_gene["Gene|Protein"] = [query_text]
    # print(f">>> [Fuzzy] 指定类别：仅 Gene|Protein（query_text='{query_text}'）")
    # df_gene = client.search_fuzzy(payload_fuzzy_gene)
    # print(df_gene.head(10))
    # print(f"Gene|Protein 命中数量：{len(df_gene)}\n")

    # ===============================
    # C) 未指定类别 → 全类别 Hybrid（字符串+向量，带回退）
    # ===============================
    # payload_hybrid_auto = {
    #     "query_id": "demo-auto-hybrid",
    #     "options": {"top_k": 10, "return_diagnostics": True, "query_text": query_text},
    # }
    # print(f">>> [Hybrid] 未指定类别：全类别混合检索（query_text='{query_text}'）")
    # std_all, diag_all, resp_all = client.search_hybrid(payload_hybrid_auto)
    # print("standardized_df:")
    # print(std_all)
    # print("\ndiagnostics_df:")
    # print(diag_all.head(10))
    # print("\nsummary/resp:")
    # print(resp_all, "\n")

    # ===============================
    # D) 指定单一类别 Hybrid：仅 Gene|Protein
    # ===============================
    payload_hybrid_gene = {
        "query_id": "demo-gene-hybrid",
        "options": {"top_k": 10, "return_diagnostics": True},
        "Gene|Protein": [query_text],
    }
    print(f">>> [Hybrid] 指定单一类别：仅 Gene|Protein（query_text='{query_text}'）")
    std_g, diag_g, resp_g = client.search_hybrid(
        payload_hybrid_gene,
        type_mix_override={"Gene|Protein": 1.0},
    )
    print(diag_g["primary_name"].to_list())

    # print("standardized_df:")
    # print("type(std_g):", type(std_g))
    # print(std_g)exit
    # print("\ndiagnostics_df:")
    # print("type(diag_g):", type(diag_g))
    # print(diag_g.head(10))
    # print("\nsummary/resp:")
    # print("type(resp_g):", type(resp_g))
    # print(resp_g, "\n")

    # ===============================
    # E) 指定多类别 Hybrid：Gene + Disease
    # ===============================
    # payload_hybrid_mix = {
    #     "query_id": "demo-mix-hybrid",
    #     "options": {"top_k": 10, "return_diagnostics": True},
    #     "Gene|Protein": [query_text],
    #     "Disease|Phenotype": [query_text],
    # }
    # print(f">>> [Hybrid] 多类别：Gene 与 Disease（query_text='{query_text}'）")
    # std_m, diag_m, resp_m = client.search_hybrid(
    #     payload_hybrid_mix,
    #     type_mix_override={"Gene|Protein": 0.6, "Disease|Phenotype": 0.4},
    # )
    # print("standardized_df:")
    # print(std_m)
    # print("\ndiagnostics_df:")
    # print(diag_m.head(10))
    # print("\nsummary/resp:")
    # print(resp_m)
    # for cell in query_cells:
    #     payload_hybrid_gene = {
    #         "query_id": "demo-gene-hybrid",
    #         "options": {"top_k": 10, "return_diagnostics": True},
    #         "Cell|Tissue": [cell],
    #     }
    #     print(f">>> [Hybrid] 指定单一类别：仅 Cell|Tissue（query_text='{cell}'）")
    #     std_g, diag_g, resp_g = client.search_hybrid(
    #         payload_hybrid_gene,
    #         type_mix_override={"Cell|Tissue": 1.0},
    #     )
    #     print("standardized_df:")
    #     print(std_g)
    #     print("\ndiagnostics_df:")
    #     print(diag_g.head(10))
    #     print("\nsummary/resp:")
    #     print(resp_g, "\n")