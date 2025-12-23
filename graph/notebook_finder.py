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
                "keywords": ["细胞注释", "细胞类型", "聚类", "差异表达基因", "单细胞", "scRNA-seq"],
                "functions": ["细胞类型识别", "细胞注释", "细胞分类", "单细胞分析"],
                "synonyms": ["细胞标注", "细胞识别", "细胞分型", "细胞类型鉴定", "细胞类型标记"],
                "scenarios": ["单细胞测序数据分析", "细胞群体识别", "细胞类型分类", "免疫细胞分析"],
                "common_questions": [
                    "如何进行细胞注释", "怎样识别细胞类型", "细胞分类怎么做", 
                    "单细胞数据如何注释", "BRICK细胞注释", "细胞类型识别方法",
                    "我想给细胞打标签", "如何确定细胞身份", "细胞群体分析"
                ],
                "technical_terms": ["marker基因", "细胞标记物", "细胞表型", "细胞身份", "细胞群体"]
            },
            "Celltype_refinement.ipynb": {
                "title": "细胞类型精化",
                "description": "对已注释的细胞类型进行进一步精化和优化，识别细胞亚型",
                "keywords": ["细胞精化", "细胞类型优化", "注释优化", "细胞亚型", "细分", "精细化"],
                "functions": ["细胞类型精化", "细胞亚型分析", "注释优化"],
                "synonyms": ["细胞类型细分", "细胞亚群分析", "细胞类型优化", "注释精细化"],
                "scenarios": ["细胞亚型识别", "细胞类型细分", "注释结果优化", "细胞群体精细分析"],
                "common_questions": [
                    "如何细分细胞类型", "细胞亚型怎么分析", "注释结果如何优化",
                    "细胞类型太粗糙怎么办", "如何识别细胞亚群", "细胞分类不够精细",
                    "想要更详细的细胞分类", "细胞类型需要进一步划分"
                ],
                "technical_terms": ["细胞亚型", "细胞亚群", "细胞状态", "细胞异质性"]
            },
            "Differential_Expression_Gene.ipynb": {
                "title": "差异表达基因分析",
                "description": "进行差异表达基因分析，识别不同条件、组别或细胞类型间的基因表达差异",
                "keywords": ["差异表达", "基因分析", "DEG", "基因表达", "统计分析", "差异基因", "表达差异"],
                "functions": ["差异基因分析", "基因表达分析", "统计检验"],
                "synonyms": ["差异基因分析", "表达差异分析", "基因差异表达", "DEG分析"],
                "scenarios": ["疾病vs正常对比", "处理vs对照分析", "不同细胞类型比较", "时间序列分析"],
                "common_questions": [
                    "如何找差异表达基因", "怎样分析基因表达差异", "DEG分析怎么做",
                    "哪些基因表达有差异", "如何比较基因表达", "差异基因筛选",
                    "基因表达变化分析", "找显著差异的基因", "基因上调下调分析"
                ],
                "technical_terms": ["fold change", "p值", "FDR", "上调基因", "下调基因", "显著性检验"]
            },
            "Drug_Discovery.ipynb": {
                "title": "药物发现",
                "description": "基于BRICK进行药物发现分析，探索潜在的治疗药物和药物靶点",
                "keywords": ["药物发现", "药物筛选", "治疗", "药物分析", "COVID-19", "药物靶点", "治疗靶点"],
                "functions": ["药物发现", "药物筛选", "治疗分析", "药物靶点"],
                "synonyms": ["药物开发", "药物研发", "治疗药物发现", "药物候选物筛选"],
                "scenarios": ["疾病治疗药物筛选", "药物重定位", "新药发现", "药物机制研究"],
                "common_questions": [
                    "如何发现新药", "怎样筛选药物", "药物发现流程", "治疗药物分析",
                    "找潜在的治疗药物", "药物靶点分析", "药物机制研究",
                    "什么药物可以治疗", "药物重定位分析", "寻找治疗方案"
                ],
                "technical_terms": ["药物靶点", "药物机制", "药物相互作用", "治疗效果", "药物重定位"]
            },
            "Drug_Discovery_NK.ipynb": {
                "title": "NK细胞药物发现",
                "description": "针对NK细胞的药物发现分析，研究NK细胞相关的治疗策略",
                "keywords": ["NK细胞", "药物发现", "免疫细胞", "药物筛选", "免疫治疗", "细胞治疗"],
                "functions": ["NK细胞分析", "免疫细胞药物发现"],
                "synonyms": ["自然杀伤细胞药物", "NK细胞治疗", "免疫细胞药物"],
                "scenarios": ["免疫治疗", "细胞治疗", "肿瘤免疫", "NK细胞功能增强"],
                "common_questions": [
                    "NK细胞药物发现", "免疫细胞治疗药物", "NK细胞相关药物",
                    "如何增强NK细胞功能", "NK细胞治疗方案", "免疫治疗药物筛选"
                ],
                "technical_terms": ["自然杀伤细胞", "细胞毒性", "免疫监视", "肿瘤杀伤"]
            },
            "Enrichment.ipynb": {
                "title": "基因富集分析",
                "description": "基于BRICK进行基因集富集分析，解析基因功能、生物学通路和分子机制",
                "keywords": ["基因富集", "功能分析", "通路分析", "GO分析", "KEGG", "功能注释", "通路富集"],
                "functions": ["基因富集分析", "功能注释", "通路分析", "GO分析"],
                "synonyms": ["功能富集", "通路富集", "基因集分析", "功能注释分析"],
                "scenarios": ["基因功能解释", "生物学通路分析", "分子机制研究", "功能分类"],
                "common_questions": [
                    "基因富集分析怎么做", "如何分析基因功能", "通路分析方法",
                    "基因的生物学功能", "GO富集分析", "KEGG通路分析",
                    "基因集功能注释", "生物学过程分析", "分子功能研究"
                ],
                "technical_terms": ["Gene Ontology", "生物学过程", "分子功能", "细胞组分", "信号通路"]
            },
            "GRN.ipynb": {
                "title": "基因调控网络",
                "description": "构建和分析基因调控网络，解释细胞轨迹推断后的调控关系和转录调控机制",
                "keywords": ["基因调控网络", "调控关系", "网络分析", "转录因子", "靶基因", "调控机制"],
                "functions": ["基因调控网络", "调控分析", "网络构建", "转录调控"],
                "synonyms": ["转录调控网络", "基因网络", "调控网络分析", "转录网络"],
                "scenarios": ["转录调控机制研究", "基因网络构建", "调控关系分析", "网络拓扑分析"],
                "common_questions": [
                    "如何构建基因调控网络", "基因调控关系分析", "转录因子靶基因",
                    "基因网络分析", "调控机制研究", "转录调控网络",
                    "基因间的调控关系", "网络拓扑分析", "调控网络构建"
                ],
                "technical_terms": ["转录因子", "靶基因", "调控元件", "网络拓扑", "调控强度"]
            },
            "Get_start.ipynb": {
                "title": "BRICK入门教程",
                "description": "BRICK的基础使用教程，展示基本功能和用法，适合初学者快速上手",
                "keywords": ["入门", "教程", "基础", "开始", "配置", "查询", "初学者", "快速上手"],
                "functions": ["入门教程", "基础教程", "BRICK使用", "配置教程"],
                "synonyms": ["新手教程", "基础指南", "快速入门", "使用指南"],
                "scenarios": ["初次使用BRICK", "学习基本操作", "环境配置", "功能演示"],
                "common_questions": [
                    "如何开始使用BRICK", "BRICK怎么用", "新手教程",
                    "BRICK入门指南", "基础操作教程", "快速上手",
                    "第一次使用BRICK", "BRICK配置方法", "基本功能介绍"
                ],
                "technical_terms": ["知识图谱", "图数据库", "查询语言", "配置文件"]
            },
            "Proteome.ipynb": {
                "title": "蛋白质组学分析",
                "description": "分析蛋白质组数据，研究疾病条件下的蛋白质相互作用网络变化和蛋白质功能",
                "keywords": ["蛋白质组", "蛋白质相互作用", "疾病机制", "蛋白质网络", "蛋白质功能"],
                "functions": ["蛋白质组分析", "蛋白质相互作用", "蛋白质网络分析"],
                "synonyms": ["蛋白组学", "蛋白质分析", "蛋白质研究", "蛋白质功能分析"],
                "scenarios": ["蛋白质功能研究", "疾病蛋白质分析", "蛋白质网络构建", "蛋白质相互作用研究"],
                "common_questions": [
                    "蛋白质组学分析", "蛋白质相互作用分析", "蛋白质网络研究",
                    "蛋白质功能分析", "疾病相关蛋白质", "蛋白质数据分析",
                    "蛋白质组数据处理", "蛋白质网络构建"
                ],
                "technical_terms": ["蛋白质相互作用", "蛋白质复合物", "蛋白质域", "蛋白质修饰"]
            },
            "Spatial_Transcriptomics.ipynb": {
                "title": "空间转录组学",
                "description": "空间转录组学数据分析，研究基因在组织空间中的表达模式和空间异质性",
                "keywords": ["空间转录组", "空间分析", "组织结构", "空间表达", "空间异质性", "组织切片"],
                "functions": ["空间转录组分析", "空间基因表达", "组织空间分析"],
                "synonyms": ["空间组学", "空间基因组学", "组织空间分析", "空间表达分析"],
                "scenarios": ["组织结构分析", "空间表达模式", "细胞空间分布", "组织功能区域"],
                "common_questions": [
                    "空间转录组分析", "空间基因表达", "组织空间结构",
                    "空间异质性分析", "细胞空间分布", "组织切片分析",
                    "空间表达模式", "组织功能区域分析"
                ],
                "technical_terms": ["空间坐标", "组织切片", "空间分辨率", "空间聚类"]
            },
            "Trajectory.ipynb": {
                "title": "细胞轨迹推断",
                "description": "使用BRICK解释细胞发育轨迹，分析细胞分化过程和发育路径",
                "keywords": ["细胞轨迹", "发育轨迹", "细胞分化", "伪时间", "轨迹分析", "发育路径"],
                "functions": ["细胞轨迹分析", "发育轨迹", "细胞分化分析", "伪时间分析"],
                "synonyms": ["细胞发育分析", "分化轨迹", "细胞命运", "发育过程分析"],
                "scenarios": ["细胞分化研究", "发育生物学", "细胞命运决定", "干细胞分化"],
                "common_questions": [
                    "细胞轨迹分析", "细胞分化过程", "发育轨迹推断",
                    "伪时间分析", "细胞命运分析", "分化路径研究",
                    "细胞发育过程", "轨迹推断方法", "细胞状态转换"
                ],
                "technical_terms": ["伪时间", "分化轨迹", "细胞状态", "发育路径", "细胞命运"]
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