"""
功能：串联混合检索流程的各个步骤，对外提供统一的查询入口，并输出标准化结果与诊断信息。
构建逻辑：预期实现为类或函数，内部依次调用 `adapters`、`string_client`、`vector_client`、`fusion`、`schema` 中定义的数据模型，以及可选的 LLM 校验模块。
数据流：上游接收外部传入的查询 JSON；下游输出 `schema` 描述的标准化结果（standardized、diagnostics、logs），供调用方消费。
调用链：外部调用 `HybridEntitySearcher`（占位）后，会依次触发输入适配、字符串检索、向量检索、分数融合、标准化输出和日志记录。
"""

from __future__ import annotations

import logging
from typing import Dict, Mapping, Optional

from elasticsearch import Elasticsearch  # type: ignore

from . import adapters, fusion, logging_utils, settings, string_client, vector_client
from .schema import DiagnosticEntry, HybridHit, HybridResponse, NormalizedQuery
from .settings import HybridSearchConfig

logger = logging.getLogger(__name__)


class HybridEntitySearcher:
    """
    功能：混合检索的对外入口类，负责生命周期内的 ES 客户端、配置管理以及检索流程编排。
    构建逻辑：在初始化阶段加载配置与 ES 客户端，`search` 方法串联 `adapters`、通道检索与融合流程。
    数据流：
        - 输入：原始 payload（JSON -> dict）。
        - 输出：`HybridResponse`，包含标准化结果、诊断信息与日志。
    调用链：由上层业务或服务通过实例化该类并调用 `search` 方法获取检索结果。
    """

    def __init__(self, es_client: Elasticsearch, config: HybridSearchConfig | None = None) -> None:
        """
        功能：初始化混合检索器，保存 ES 客户端与配置。
        构建逻辑：若未显式提供配置，则调用 `settings.get_search_config()` 获取默认配置。
        出参入参：入参为 ES 客户端和可选配置；无返回值。
        数据流：为后续检索流程提供运行所需的信息。
        模块内调用链：`search`、`_run_channels` 等方法依赖成员属性。
        模块外调用链：上层代码实例化此类。
        """
        self._es = es_client
        self._config = config or settings.get_search_config()

    def _run_channels(self, normalized: NormalizedQuery) -> tuple[list[HybridHit], list[HybridHit]]:
        """
        功能：执行字符串与向量两个检索通道，返回各自命中列表。
        构建逻辑：调用 `string_client.search_string_channel` 与 `vector_client.search_vector_channel`，并根据配置决定是否启用向量检索。
        出参入参：入参为 `NormalizedQuery`；出参为 (string_hits, vector_hits)。
        数据流：结果将传递给 `fusion.apply_fusion` 做融合。
        模块内调用链：供 `search` 调用。
        模块外调用链：无。
        """
        string_hits = string_client.search_string_channel(
            self._es,
            normalized_query=normalized,
            index_name=self._config.string_index_name,
        )
        vector_hits: list[HybridHit] = []
        if self._config.vector_index_name:
            vector_hits = vector_client.search_vector_channel(
                self._es,
                normalized_query=normalized,
                config=self._config,
            )
        return string_hits, vector_hits

    def _build_response(
        self,
        normalized: NormalizedQuery,
        fused_hits: list[HybridHit],
    ) -> HybridResponse:
        """
        功能：将融合后的命中转换为最终响应，包括标准化结果、诊断信息与日志。
        构建逻辑：
            1. 根据类型键整理标准化列表；
            2. 若选项开启诊断/调试，生成 `DiagnosticEntry` 列表；
            3. 记录必要的日志信息并封装为 `HybridResponse`。
        出参入参：入参为标准化查询与融合命中列表；出参为 `HybridResponse`。
        数据流：作为最终结果返回给调用方。
        模块内调用链：`search` 调用本方法。
        模块外调用链：无。
        """
        standardized: Dict[str, list[str]] = {}
        for hit in fused_hits:
            standardized.setdefault(hit.type_key, []).append(hit.entity_id)

        diagnostics: Optional[list[DiagnosticEntry]] = None
        if normalized.options.return_diagnostics or normalized.options.debug:
            diagnostics = [
                DiagnosticEntry(
                    type_key=hit.type_key,
                    entity_id=hit.entity_id,
                    primary_name=hit.primary_name,
                    node_type=hit.node_type,
                    final_score=hit.final_score,
                    channel_scores=hit.scores,
                    matched_alias=hit.matched_alias,
                    extra=hit.metadata,
                )
                for hit in fused_hits
            ]

        logs: Optional[Dict[str, object]] = None
        if normalized.options.debug:
            logs = {
                "string_hits": sum(1 for hit in fused_hits if hit.scores.string_score > 0),
                "vector_hits": sum(1 for hit in fused_hits if hit.scores.vector_score > 0),
                "top_k": normalized.options.top_k,
                "type_mix_override": normalized.options.type_mix_override,
            }

        return HybridResponse(
            query_id=normalized.query_id,
            standardized=standardized,
            diagnostics=diagnostics,
            logs=logs,
        )

    def search(self, payload: Mapping[str, object]) -> HybridResponse:
        """
        功能：执行一次混合检索流程，外部仅需提供原始 payload。
        构建逻辑：
            1. 调用 `adapters.normalize_payload` 标准化输入；
            2. 调用 `_run_channels` 获取两通道命中；
            3. 调用 `fusion.apply_fusion` 合并；
            4. 调用 `_build_response` 构造返回值；
            5. 调用 `logging_utils` 记录必要日志。
        出参入参：入参为原始 payload；出参为 `HybridResponse`。
        数据流：完整执行混合检索的数据路径。
        模块内调用链：依赖 `_run_channels`、`_build_response`。
        模块外调用链：业务方调用本方法以获取检索结果。
        """
        payload_dict = dict(payload)
        normalized = adapters.normalize_payload(payload_dict, self._config)
        string_hits, vector_hits = self._run_channels(normalized)

        type_mix = dict(self._config.type_mix)
        if normalized.options.type_mix_override:
            type_mix.update(normalized.options.type_mix_override)

        fused_hits = fusion.apply_fusion(
            normalized_query=normalized,
            string_hits=string_hits,
            vector_hits=vector_hits,
            type_mix=type_mix,
            alpha=self._config.alpha,
            beta=self._config.beta,
        )

        response = self._build_response(normalized, fused_hits)

        if response.diagnostics and logging_utils.should_trigger_llm(
            response.diagnostics, threshold=0.1
        ):
            response = logging_utils.run_llm_verification(payload_dict, response)

        logging_utils.record_query_log(
            payload_dict,
            response,
            extra={
                "string_hits": len(string_hits),
                "vector_hits": len(vector_hits),
            },
        )
        return response


__all__ = ["HybridEntitySearcher"]
