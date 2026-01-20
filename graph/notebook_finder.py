#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import json
import re
import pickle
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import jieba
import jieba.analyse


class NotebookFinder:
    """基于RAG的Notebook查找器"""
    
    def __init__(self, cache_dir: str = "./notebook_cache"):
        """
        初始化NotebookFinder
        
        Args:
            cache_dir: 缓存目录路径
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        # 初始化中文分词
        jieba.initialize()
        
        # 预定义的notebook功能描述映射 - 增强版，支持多元化问题解析
        self.notebook_descriptions = {
            "Celltype_annotation.ipynb": {
                "title": "细胞类型注释",
                "description": "使用BRICK进行细胞类型注释，基于聚类结果的差异表达基因定位对应细胞类型",
                "keywords": ["细胞注释", "细胞类型", "聚类", "差异表达基因", "单细胞", "scRNA-seq", "cell annotation", "cell type", "clustering", "differential expression", "single cell"],
                "functions": ["细胞类型识别", "细胞注释", "细胞分类", "单细胞分析", "cell type identification", "cell annotation", "cell classification", "single cell analysis"],
                "synonyms": ["细胞标注", "细胞识别", "细胞分型", "细胞类型鉴定", "细胞类型标记", "cell labeling", "cell identification", "cell typing"],
                "scenarios": ["单细胞测序数据分析", "细胞群体识别", "细胞类型分类", "免疫细胞分析", "single cell sequencing", "cell population identification", "immune cell analysis"],
                "common_questions": [
                    "如何进行细胞注释", "怎样识别细胞类型", "细胞分类怎么做", 
                    "单细胞数据如何注释", "BRICK细胞注释", "细胞类型识别方法",
                    "我想给细胞打标签", "如何确定细胞身份", "细胞群体分析",
                    "how to annotate cells", "cell type identification", "cell annotation method",
                    "single cell annotation", "identify cell types", "cell labeling"
                ],
                "technical_terms": ["marker基因", "细胞标记物", "细胞表型", "细胞身份", "细胞群体", "marker gene", "cell marker", "cell phenotype", "cell identity"]
            },
            "Celltype_refinement.ipynb": {
                "title": "细胞类型精化",
                "description": "对已注释的细胞类型进行进一步精化和优化，识别细胞亚型",
                "keywords": ["细胞精化", "细胞类型优化", "注释优化", "细胞亚型", "细分", "精细化", "cell refinement", "cell type optimization", "annotation refinement", "cell subtype"],
                "functions": ["细胞类型精化", "细胞亚型分析", "注释优化", "cell type refinement", "cell subtype analysis", "annotation optimization"],
                "synonyms": ["细胞类型细分", "细胞亚群分析", "细胞类型优化", "注释精细化", "cell type subdivision", "cell subpopulation analysis"],
                "scenarios": ["细胞亚型识别", "细胞类型细分", "注释结果优化", "细胞群体精细分析", "cell subtype identification", "annotation optimization"],
                "common_questions": [
                    "如何细分细胞类型", "细胞亚型怎么分析", "注释结果如何优化",
                    "细胞类型太粗糙怎么办", "如何识别细胞亚群", "细胞分类不够精细",
                    "想要更详细的细胞分类", "细胞类型需要进一步划分",
                    "refine cell types", "cell subtype analysis", "optimize annotation",
                    "cell classification refinement", "identify cell subtypes"
                ],
                "technical_terms": ["细胞亚型", "细胞亚群", "细胞状态", "细胞异质性", "cell subtype", "cell subpopulation", "cell state", "cell heterogeneity"]
            },
            "Differential_Expression_Gene.ipynb": {
                "title": "差异表达基因分析",
                "description": "进行差异表达基因分析，识别不同条件、组别或细胞类型间的基因表达差异",
                "keywords": ["差异表达", "基因分析", "DEG", "基因表达", "统计分析", "差异基因", "表达差异", "differential expression", "gene analysis", "gene expression", "statistical analysis"],
                "functions": ["差异基因分析", "基因表达分析", "统计检验", "differential gene analysis", "gene expression analysis", "statistical testing"],
                "synonyms": ["差异基因分析", "表达差异分析", "基因差异表达", "DEG分析", "DEG analysis", "expression difference analysis"],
                "scenarios": ["疾病vs正常对比", "处理vs对照分析", "不同细胞类型比较", "时间序列分析", "disease vs normal", "treatment vs control", "cell type comparison"],
                "common_questions": [
                    "如何找差异表达基因", "怎样分析基因表达差异", "DEG分析怎么做",
                    "哪些基因表达有差异", "如何比较基因表达", "差异基因筛选",
                    "基因表达变化分析", "找显著差异的基因", "基因上调下调分析",
                    "find DEGs", "differential gene expression", "how to do DEG analysis",
                    "gene expression comparison", "identify differentially expressed genes"
                ],
                "technical_terms": ["fold change", "p值", "FDR", "上调基因", "下调基因", "显著性检验", "p-value", "upregulated", "downregulated", "significance test"]
            },
            "Drug_Discovery.ipynb": {
                "title": "药物发现",
                "description": "基于BRICK进行药物发现分析，探索潜在的治疗药物和药物靶点",
                "keywords": ["药物发现", "药物筛选", "治疗", "药物分析", "药物靶点", "治疗靶点", "drug discovery", "drug screening", "therapy", "drug target", "therapeutic target"],
                "functions": ["药物发现", "药物筛选", "治疗分析", "药物靶点", "drug discovery", "drug screening", "therapeutic analysis"],
                "synonyms": ["药物开发", "药物研发", "治疗药物发现", "药物候选物筛选", "drug development", "therapeutic drug discovery", "drug candidate screening"],
                "scenarios": ["疾病治疗药物筛选", "药物重定位", "新药发现", "药物机制研究", "disease treatment", "drug repurposing", "new drug discovery"],
                "common_questions": [
                    "如何发现新药", "怎样筛选药物", "药物发现流程", "治疗药物分析",
                    "找潜在的治疗药物", "药物靶点分析", "药物机制研究",
                    "什么药物可以治疗", "药物重定位分析", "寻找治疗方案",
                    "discover new drugs", "drug screening", "therapeutic drug analysis",
                    "find potential drugs", "drug target analysis", "drug repositioning"
                ],
                "technical_terms": ["药物靶点", "药物机制", "药物相互作用", "治疗效果", "药物重定位", "drug target", "drug mechanism", "drug interaction", "therapeutic effect"]
            },
            "Enrichment.ipynb": {
                "title": "基因富集分析",
                "description": "基于BRICK进行基因集富集分析，解析基因功能、生物学通路和分子机制",
                "keywords": ["基因富集", "功能分析", "通路分析", "GO分析", "KEGG", "功能注释", "通路富集", "gene enrichment", "functional analysis", "pathway analysis", "GO analysis", "functional annotation"],
                "functions": ["基因富集分析", "功能注释", "通路分析", "GO分析", "gene enrichment analysis", "functional annotation", "pathway analysis"],
                "synonyms": ["功能富集", "通路富集", "基因集分析", "功能注释分析", "functional enrichment", "pathway enrichment", "gene set analysis"],
                "scenarios": ["基因功能解释", "生物学通路分析", "分子机制研究", "功能分类", "gene function interpretation", "biological pathway analysis", "molecular mechanism"],
                "common_questions": [
                    "基因富集分析怎么做", "如何分析基因功能", "通路分析方法",
                    "基因的生物学功能", "GO富集分析", "KEGG通路分析",
                    "基因集功能注释", "生物学过程分析", "分子功能研究",
                    "gene enrichment analysis", "GO enrichment", "KEGG pathway",
                    "functional annotation", "pathway enrichment"
                ],
                "technical_terms": ["Gene Ontology", "生物学过程", "分子功能", "细胞组分", "信号通路", "biological process", "molecular function", "cellular component", "signaling pathway"]
            },
            "GeneRegulatoryNetwork.ipynb": {
                "title": "基因调控网络",
                "description": "构建和分析基因调控网络，解释细胞轨迹推断后的调控关系和转录调控机制",
                "keywords": ["基因调控网络", "调控关系", "网络分析", "转录因子", "靶基因", "调控机制", "GRN", "gene regulatory network", "regulatory relationship", "network analysis", "transcription factor", "target gene"],
                "functions": ["基因调控网络", "调控分析", "网络构建", "转录调控", "gene regulatory network", "regulatory analysis", "network construction"],
                "synonyms": ["转录调控网络", "基因网络", "调控网络分析", "转录网络", "GRN分析", "transcriptional regulatory network", "gene network", "GRN analysis"],
                "scenarios": ["转录调控机制研究", "基因网络构建", "调控关系分析", "网络拓扑分析", "transcriptional regulation", "network construction", "network topology"],
                "common_questions": [
                    "如何构建基因调控网络", "基因调控关系分析", "转录因子靶基因",
                    "基因网络分析", "调控机制研究", "转录调控网络",
                    "基因间的调控关系", "网络拓扑分析", "调控网络构建", "GRN分析",
                    "gene regulatory network", "transcription factor analysis", "regulatory relationship",
                    "network construction", "GRN analysis"
                ],
                "technical_terms": ["转录因子", "靶基因", "调控元件", "网络拓扑", "调控强度", "regulon", "transcription factor", "target gene", "regulatory element", "network topology"]
            },
            "GWAS_Causal_SNPs.ipynb": {
                "title": "GWAS因果SNP分析",
                "description": "基于BRICK进行全基因组关联研究(GWAS)的因果SNP识别和分析，找出与疾病或性状相关的关键遗传变异位点",
                "keywords": ["GWAS", "SNP", "因果变异", "遗传变异", "全基因组关联", "疾病关联", "遗传标记", "causal variant", "genetic variation", "genome-wide association", "disease association"],
                "functions": ["GWAS分析", "因果SNP识别", "遗传变异分析", "疾病关联分析", "GWAS analysis", "causal SNP identification", "genetic variant analysis"],
                "synonyms": ["全基因组关联分析", "单核苷酸多态性分析", "致病变异识别", "遗传位点分析", "genome-wide association study", "single nucleotide polymorphism", "pathogenic variant"],
                "scenarios": ["疾病易感基因识别", "遗传风险评估", "致病位点筛选", "表型关联研究", "disease susceptibility", "genetic risk assessment", "phenotype association"],
                "common_questions": [
                    "如何识别因果SNP", "GWAS分析怎么做", "如何找致病位点",
                    "遗传变异分析方法", "SNP与疾病的关联", "如何筛选关键位点",
                    "致病基因定位", "遗传标记分析", "疾病易感位点识别",
                    "identify causal SNPs", "GWAS analysis", "find pathogenic variants",
                    "genetic variation analysis", "SNP disease association"
                ],
                "technical_terms": ["单核苷酸多态性", "连锁不平衡", "显著性阈值", "效应量", "遗传度", "精细定位", "single nucleotide polymorphism", "linkage disequilibrium", "significance threshold", "effect size", "heritability"]
            },
            "GWAS_Predict_Phenotype.ipynb": {
                "title": "GWAS表型预测",
                "description": "利用GWAS数据和BRICK知识图谱进行表型预测，基于遗传变异信息预测个体性状和疾病风险",
                "keywords": ["GWAS", "表型预测", "遗传风险", "性状预测", "多基因评分", "风险预测", "基因型", "phenotype prediction", "genetic risk", "trait prediction", "polygenic score", "risk prediction", "genotype"],
                "functions": ["表型预测", "风险评分", "遗传风险评估", "性状预测模型", "phenotype prediction", "risk scoring", "genetic risk assessment", "trait prediction model"],
                "synonyms": ["性状预测", "疾病风险预测", "遗传评分", "多基因风险评分", "PRS分析", "trait prediction", "disease risk prediction", "genetic scoring", "PRS analysis"],
                "scenarios": ["疾病风险预测", "个体化医疗", "精准医学", "遗传咨询", "健康风险评估", "disease risk prediction", "personalized medicine", "precision medicine", "genetic counseling"],
                "common_questions": [
                    "如何预测表型", "疾病风险怎么计算", "遗传风险评估方法",
                    "多基因评分怎么做", "如何预测性状", "基因型如何预测表型",
                    "个体疾病风险预测", "遗传信息如何用于预测", "风险评分计算",
                    "phenotype prediction", "disease risk calculation", "genetic risk assessment",
                    "polygenic score", "predict traits", "PRS analysis"
                ],
                "technical_terms": ["多基因风险评分", "PRS", "遗传风险", "预测模型", "基因型-表型关联", "遗传贡献", "polygenic risk score", "genetic risk", "prediction model", "genotype-phenotype association"]
            },
            "Proteome.ipynb": {
                "title": "蛋白质组学分析",
                "description": "分析蛋白质组数据，研究疾病条件下的蛋白质相互作用网络变化和蛋白质功能",
                "keywords": ["蛋白质组", "蛋白质相互作用", "疾病机制", "蛋白质网络", "蛋白质功能", "proteome", "protein interaction", "disease mechanism", "protein network", "protein function"],
                "functions": ["蛋白质组分析", "蛋白质相互作用", "蛋白质网络分析", "proteomics analysis", "protein interaction", "protein network analysis"],
                "synonyms": ["蛋白组学", "蛋白质分析", "蛋白质研究", "蛋白质功能分析", "proteomics", "protein analysis", "protein research"],
                "scenarios": ["蛋白质功能研究", "疾病蛋白质分析", "蛋白质网络构建", "蛋白质相互作用研究", "protein function study", "disease protein analysis", "protein network construction"],
                "common_questions": [
                    "蛋白质组学分析", "蛋白质相互作用分析", "蛋白质网络研究",
                    "蛋白质功能分析", "疾病相关蛋白质", "蛋白质数据分析",
                    "蛋白质组数据处理", "蛋白质网络构建",
                    "proteomics analysis", "protein interaction analysis", "protein network",
                    "protein function analysis", "disease protein"
                ],
                "technical_terms": ["蛋白质相互作用", "蛋白质复合物", "蛋白质域", "蛋白质修饰", "protein interaction", "protein complex", "protein domain", "protein modification"]
            },
            "Spatial_Transcriptomics.ipynb": {
                "title": "空间转录组学",
                "description": "空间转录组学数据分析，研究基因在组织空间中的表达模式和空间异质性",
                "keywords": ["空间转录组", "空间分析", "组织结构", "空间表达", "空间异质性", "组织切片", "spatial transcriptomics", "spatial analysis", "tissue structure", "spatial expression", "spatial heterogeneity", "tissue section"],
                "functions": ["空间转录组分析", "空间基因表达", "组织空间分析", "spatial transcriptomics analysis", "spatial gene expression", "tissue spatial analysis"],
                "synonyms": ["空间组学", "空间基因组学", "组织空间分析", "空间表达分析", "spatial omics", "spatial genomics", "tissue spatial analysis", "spatial expression analysis"],
                "scenarios": ["组织结构分析", "空间表达模式", "细胞空间分布", "组织功能区域", "tissue structure analysis", "spatial expression pattern", "cell spatial distribution", "tissue functional region"],
                "common_questions": [
                    "空间转录组分析", "空间基因表达", "组织空间结构",
                    "空间异质性分析", "细胞空间分布", "组织切片分析",
                    "空间表达模式", "组织功能区域分析",
                    "spatial transcriptomics analysis", "spatial gene expression", "tissue spatial structure",
                    "spatial heterogeneity analysis", "cell spatial distribution", "tissue section analysis",
                    "spatial expression pattern", "tissue functional region analysis"
                ],
                "technical_terms": ["空间坐标", "组织切片", "空间分辨率", "空间聚类", "spatial coordinates", "tissue section", "spatial resolution", "spatial clustering"]
            },
            "Trajectory.ipynb": {
                "title": "细胞轨迹推断",
                "description": "使用BRICK解释细胞发育轨迹，分析细胞分化过程和发育路径",
                "keywords": ["细胞轨迹", "发育轨迹", "细胞分化", "伪时间", "轨迹分析", "发育路径", "cell trajectory", "developmental trajectory", "cell differentiation", "pseudotime", "trajectory analysis", "developmental pathway"],
                "functions": ["细胞轨迹分析", "发育轨迹", "细胞分化分析", "伪时间分析", "cell trajectory analysis", "developmental trajectory", "cell differentiation analysis", "pseudotime analysis"],
                "synonyms": ["细胞发育分析", "分化轨迹", "细胞命运", "发育过程分析", "cell development analysis", "differentiation trajectory", "cell fate", "developmental process analysis"],
                "scenarios": ["细胞分化研究", "发育生物学", "细胞命运决定", "干细胞分化", "cell differentiation study", "developmental biology", "cell fate determination", "stem cell differentiation"],
                "common_questions": [
                    "细胞轨迹分析", "细胞分化过程", "发育轨迹推断",
                    "伪时间分析", "细胞命运分析", "分化路径研究",
                    "细胞发育过程", "轨迹推断方法", "细胞状态转换",
                    "cell trajectory analysis", "cell differentiation process", "trajectory inference",
                    "pseudotime analysis", "cell fate analysis", "differentiation pathway",
                    "cell development process", "trajectory inference method", "cell state transition"
                ],
                "technical_terms": ["伪时间", "分化轨迹", "细胞状态", "发育路径", "细胞命运", "pseudotime", "differentiation trajectory", "cell state", "developmental pathway", "cell fate"]
            },
            "Terminology_query.ipynb": {
                "title": "名词术语查询",
                "description": "使用BRICK查询生物医学术语的定义和相关信息，回答'XXX是什么'类型的问题",
                "keywords": ["名词查询", "术语定义", "实体查询", "基因查询", "疾病查询", "细胞查询", "知识查询", "terminology query", "entity definition", "gene query", "disease query", "cell query", "knowledge query"],
                "functions": ["术语查询", "实体定义查询", "名词解释", "知识图谱查询", "邻居节点查询", "terminology lookup", "entity definition query", "term explanation", "knowledge graph query", "neighbor query"],
                "synonyms": ["定义查询", "概念查询", "实体信息查询", "术语解释", "名词定义", "definition query", "concept query", "entity information query", "term explanation", "term definition"],
                "scenarios": ["基因功能查询", "疾病信息查询", "细胞类型查询", "蛋白质信息查询", "组织信息查询", "gene function query", "disease information query", "cell type query", "protein information query", "tissue information query"],
                "common_questions": [
                    "什么是", "XXX是什么", "XXX的定义", "XXX的功能",
                    "查询XXX", "XXX的信息", "告诉我XXX", "解释XXX",
                    "XXX有哪些相关的", "与XXX相关的节点", "XXX的邻居节点",
                    "what is", "definition of", "what does XXX mean", "tell me about XXX",
                    "information about XXX", "explain XXX", "describe XXX",
                    "what is related to XXX", "neighbors of XXX"
                ],
                "technical_terms": ["实体节点", "邻居节点", "节点属性", "知识图谱", "语义查询", "entity node", "neighbor node", "node attribute", "knowledge graph", "semantic query"]
            },
            "Relation_query.ipynb": {
                "title": "实体关系查询",
                "description": "使用BRICK查询两个或多个生物医学实体之间的关系，回答'XXX与YYY有什么关系'类型的问题",
                "keywords": ["关系查询", "实体关系", "关联分析", "路径查询", "关系类型", "知识图谱关系", "relation query", "entity relationship", "association analysis", "path query", "relationship type", "knowledge graph relation"],
                "functions": ["关系查询", "实体关联查询", "路径分析", "关系类型识别", "多实体关系查询", "relation lookup", "entity association query", "path analysis", "relationship type identification", "multi-entity relation query"],
                "synonyms": ["关联查询", "实体间关系", "关系分析", "连接查询", "路径发现", "association query", "inter-entity relationship", "relationship analysis", "connection query", "path discovery"],
                "scenarios": ["基因-疾病关联", "基因-蛋白关系", "疾病-药物关系", "细胞-基因关系", "多实体关系分析", "gene-disease association", "gene-protein relationship", "disease-drug relationship", "cell-gene relationship", "multi-entity relationship analysis"],
                "common_questions": [
                    "XXX与YYY有什么关系", "XXX和YYY之间的关系", "XXX如何影响YYY",
                    "XXX是否关联YYY", "XXX与YYY的关联", "XXX和YYY的联系",
                    "XXX、YYY和ZZZ之间的关系", "查询XXX与YYY的关系",
                    "what is the relationship between XXX and YYY", "how XXX relates to YYY",
                    "relationship between XXX and YYY", "association between XXX and YYY",
                    "how does XXX affect YYY", "connection between XXX and YYY",
                    "relationship among XXX, YYY and ZZZ"
                ],
                "technical_terms": ["关系类型", "路径长度", "关系置信度", "信息来源", "关系方向", "多跳关系", "relation type", "path length", "relation confidence", "information source", "relation direction", "multi-hop relation"]
            },
            "Most_related_query.ipynb": {
                "title": "最相关实体查询",
                "description": "使用BRICK查询与指定实体最相关的特定类型实体，回答'与XXX最相关的YYY是什么'类型的问题",
                "keywords": ["最相关查询", "相关性排序", "实体排名", "邻居查询", "匹配计数", "关联度", "most related query", "relevance ranking", "entity ranking", "neighbor query", "match count", "association degree"],
                "functions": ["相关性查询", "实体排序", "最相关实体识别", "关联度计算", "邻居节点排名", "relevance query", "entity ranking", "most related entity identification", "association degree calculation", "neighbor node ranking"],
                "synonyms": ["最相关实体", "关联度最高", "最强关联", "最密切相关", "top关联", "most relevant entity", "highest association", "strongest connection", "most closely related", "top association"],
                "scenarios": ["查找相关疾病", "查找相关基因", "查找相关细胞", "查找相关通路", "查找相关药物", "find related diseases", "find related genes", "find related cells", "find related pathways", "find related drugs"],
                "common_questions": [
                    "与XXX最相关的YYY是什么", "XXX最相关的YYY", "哪个YYY与XXX最相关",
                    "XXX关联最强的YYY", "XXX最相关的疾病", "XXX最相关的基因",
                    "与XXX关系最密切的YYY", "XXX最可能相关的YYY",
                    "what YYY is most related to XXX", "most related YYY to XXX",
                    "which YYY is most associated with XXX", "top YYY related to XXX",
                    "most relevant YYY for XXX", "strongest YYY associated with XXX",
                    "which disease is most related to XXX", "which gene is most related to XXX"
                ],
                "technical_terms": ["匹配计数", "关联强度", "邻居节点", "目标实体类型", "源实体", "排序权重", "match count", "association strength", "neighbor node", "target entity type", "source entity", "ranking weight"]
            }
        }
        
        # 初始化向量化器
        self.vectorizer = None
        self.document_vectors = None
        self.notebook_files = []
        
    def _extract_notebook_content(self, notebook_path: str) -> Dict[str, str]:
        """
        提取notebook文件的内容
        
        Args:
            notebook_path: notebook文件路径
            
        Returns:
            包含提取内容的字典
        """
        try:
            with open(notebook_path, 'r', encoding='utf-8') as f:
                notebook_data = json.load(f)
            
            # 提取markdown和代码内容
            content_parts = []
            title = ""
            
            for cell in notebook_data.get('cells', []):
                if cell.get('cell_type') == 'markdown':
                    source = ''.join(cell.get('source', []))
                    # 提取标题
                    if source.startswith('#') and not title:
                        title = re.sub(r'^#+\s*', '', source.split('\n')[0])
                    content_parts.append(source)
                elif cell.get('cell_type') == 'code':
                    # 只提取注释部分
                    source = ''.join(cell.get('source', []))
                    comments = re.findall(r'#.*', source)
                    content_parts.extend(comments)
            
            content = '\n'.join(content_parts)
            
            return {
                'title': title,
                'content': content,
                'path': notebook_path
            }
            
        except Exception as e:
            print(f"Error extracting content from {notebook_path}: {e}")
            return {'title': '', 'content': '', 'path': notebook_path}
    
    def _preprocess_text(self, text: str) -> str:
        """
        预处理文本，包括中文分词和清理
        
        Args:
            text: 原始文本
            
        Returns:
            预处理后的文本
        """
        # 移除markdown标记
        text = re.sub(r'[#*`]', '', text)
        # 移除URL
        text = re.sub(r'http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\\(\\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+', '', text)
        # 移除特殊字符，保留中英文和数字
        text = re.sub(r'[^\u4e00-\u9fa5a-zA-Z0-9\s]', ' ', text)
        
        # 中文分词
        words = jieba.cut(text)
        return ' '.join(words)
    
    def _build_document_corpus(self, target_dir: str) -> List[str]:
        """
        构建文档语料库
        
        Args:
            target_dir: notebook目录路径
            
        Returns:
            文档语料库列表
        """
        corpus = []
        self.notebook_files = []
        
        for filename in os.listdir(target_dir):
            if filename.endswith('.ipynb'):
                notebook_path = os.path.join(target_dir, filename)
                
                # 获取预定义描述
                if filename in self.notebook_descriptions:
                    desc = self.notebook_descriptions[filename]
                    # 构建更丰富的文档文本，包含所有描述信息
                    text_parts = [
                        desc['title'],
                        desc['description'],
                        ' '.join(desc['keywords']),
                        ' '.join(desc['functions']),
                        ' '.join(desc.get('synonyms', [])),
                        ' '.join(desc.get('scenarios', [])),
                        ' '.join(desc.get('common_questions', [])),
                        ' '.join(desc.get('technical_terms', []))
                    ]
                    document_text = ' '.join(text_parts)
                else:
                    # 如果没有预定义描述，则提取notebook内容
                    extracted = self._extract_notebook_content(notebook_path)
                    document_text = f"{extracted['title']} {extracted['content']}"
                
                # 预处理文本
                processed_text = self._preprocess_text(document_text)
                corpus.append(processed_text)
                self.notebook_files.append(notebook_path)
        
        return corpus
    
    def _load_or_build_vectors(self, target_dir: str) -> None:
        """
        加载或构建文档向量
        
        Args:
            target_dir: notebook目录路径
        """
        cache_file = self.cache_dir / f"vectors_{hash(target_dir)}.pkl"
        
        if cache_file.exists():
            # 加载缓存的向量
            try:
                with open(cache_file, 'rb') as f:
                    cache_data = pickle.load(f)
                    self.vectorizer = cache_data['vectorizer']
                    self.document_vectors = cache_data['document_vectors']
                    self.notebook_files = cache_data['notebook_files']
                    return
            except Exception as e:
                print(f"Error loading cache: {e}")
        
        # 构建新的向量
        corpus = self._build_document_corpus(target_dir)
        
        # 使用TF-IDF向量化
        self.vectorizer = TfidfVectorizer(
            max_features=1000,
            stop_words=None,  # 不使用英文停用词，因为我们有中文
            ngram_range=(1, 2),
            min_df=1,
            max_df=0.95
        )
        
        self.document_vectors = self.vectorizer.fit_transform(corpus)
        
        # 保存到缓存
        try:
            cache_data = {
                'vectorizer': self.vectorizer,
                'document_vectors': self.document_vectors,
                'notebook_files': self.notebook_files
            }
            with open(cache_file, 'wb') as f:
                pickle.dump(cache_data, f)
        except Exception as e:
            print(f"Error saving cache: {e}")
    
    def find_notebook(self, question: str, target_dir: str, top_k: int = 3) -> List[Tuple[str, float]]:
        """
        根据问题查找最匹配的notebook文件
        
        Args:
            question: 用户问题
            target_dir: notebook目录路径
            top_k: 返回前k个最匹配的结果
            
        Returns:
            包含(文件路径, 相似度分数)的元组列表
        """
        # 加载或构建向量
        self._load_or_build_vectors(target_dir)
        
        if self.vectorizer is None or self.document_vectors is None:
            return []
        
        # 预处理问题
        processed_question = self._preprocess_text(question)
        
        # 向量化问题
        question_vector = self.vectorizer.transform([processed_question])
        
        # 计算相似度
        similarities = cosine_similarity(question_vector, self.document_vectors).flatten()
        
        # 获取top_k结果
        top_indices = np.argsort(similarities)[::-1][:top_k]
        
        results = []
        for idx in top_indices:
            if similarities[idx] > 0:  # 只返回有相似度的结果
                results.append((self.notebook_files[idx], float(similarities[idx])))
        
        return results


def find_notebook(ques: str, target_dir: str) -> Optional[str]:
    """
    便捷函数：根据问题查找最匹配的notebook文件
    
    Args:
        ques: 用户问题
        target_dir: notebook目录路径
        
    Returns:
        最匹配的notebook文件路径，如果没有找到则返回None
    """
    finder = NotebookFinder()
    results = finder.find_notebook(ques, target_dir, top_k=1)
    
    if results:
        return results[0][0]
    return None


def find_notebooks_with_scores(ques: str, target_dir: str, top_k: int = 3) -> List[Tuple[str, float]]:
    """
    便捷函数：根据问题查找多个匹配的notebook文件及其相似度分数
    
    Args:
        ques: 用户问题
        target_dir: notebook目录路径
        top_k: 返回前k个最匹配的结果
        
    Returns:
        包含(文件路径, 相似度分数)的元组列表
    """
    finder = NotebookFinder()
    return finder.find_notebook(ques, target_dir, top_k=top_k)


if __name__ == "__main__":
    # 测试示例
    test_questions = [
        # 原有测试问题
        "我想使用BRICK进行细胞注释",
        "如何进行药物发现分析",
        "怎样分析基因调控网络",
        "如何开始使用BRICK",
        "进行差异表达基因分析",
        "细胞轨迹分析怎么做",
        
        # 新增多元化测试问题
        "我想给细胞打标签",  # 细胞注释的同义表达
        "细胞分类不够精细，需要进一步细分",  # 细胞精化
        "找显著差异的基因",  # 差异表达基因的另一种表达
        "寻找治疗方案",  # 药物发现的另一种表达
        "NK细胞治疗药物筛选",  # NK细胞药物发现
        "基因的生物学功能分析",  # 基因富集分析
        "转录因子靶基因关系",  # 基因调控网络
        "新手教程",  # BRICK入门
        "蛋白质相互作用网络",  # 蛋白质组学
        "组织空间结构分析",  # 空间转录组学
        "细胞分化过程研究",  # 细胞轨迹
        
        # 更口语化的问题
        "我是新手，怎么开始用BRICK？",
        "细胞类型识别不准确怎么办？",
        "想找一些治疗COVID-19的药物",
        "基因表达有什么差异？",
        "细胞是怎么发育的？",
        "蛋白质之间有什么关系？"
    ]
    
    target_directory = "/home/lyt/checker_finallap/notebooks"
    
    print("=== Notebook Finder 测试 ===\n")
    
    for question in test_questions:
        print(f"问题: {question}")
        
        # 获取最佳匹配
        best_match = find_notebook(question, target_directory)
        if best_match:
            print(best_match)
        
        # 获取前3个匹配结果
        results = find_notebooks_with_scores(question, target_directory, top_k=3)
        print("前3个匹配结果:")
        for i, (path, score) in enumerate(results, 1):
            print(f"  {i}. {os.path.basename(path)} (相似度: {score:.3f})")
        
        print("-" * 50)