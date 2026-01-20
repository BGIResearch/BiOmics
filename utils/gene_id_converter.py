"""
基因ID识别与转换工具
用于识别h5ad文件中var属性的基因ID类型，并自动将var_names转换为Gene Symbol
"""

import re
import scanpy as sc
import pandas as pd
from collections import Counter
from typing import Dict, List, Tuple, Optional


# ========== 基因ID类型识别 ==========

def identify_gene_id_type(gene_id: str) -> str:
    """
    根据基因ID的格式判断其类型
    
    Args:
        gene_id: 单个基因ID字符串
        
    Returns:
        基因ID类型名称
    """
    gene_id = str(gene_id).strip()
    
    # Ensembl Gene ID (人/小鼠/大鼠等)
    if re.match(r'^ENS[A-Z]*G\d{11}$', gene_id):
        if gene_id.startswith('ENSG'):
            return 'Ensembl_Human_Gene'
        elif gene_id.startswith('ENSMUSG'):
            return 'Ensembl_Mouse_Gene'
        elif gene_id.startswith('ENSRNOG'):
            return 'Ensembl_Rat_Gene'
        elif gene_id.startswith('ENSDARG'):
            return 'Ensembl_Zebrafish_Gene'
        return 'Ensembl_Gene'
    
    # Ensembl Transcript ID
    if re.match(r'^ENS[A-Z]*T\d{11}$', gene_id):
        return 'Ensembl_Transcript'
    
    # Ensembl Protein ID
    if re.match(r'^ENS[A-Z]*P\d{11}$', gene_id):
        return 'Ensembl_Protein'
    
    # RefSeq mRNA ID (NM_) / ncRNA (NR_) / Protein (NP_) / Predicted (XM_, XR_, XP_)
    if re.match(r'^[NX][MRP]_\d+(\.\d+)?$', gene_id):
        return 'RefSeq'
    
    # Entrez Gene ID (纯数字, 1-8位)
    if re.match(r'^\d{1,8}$', gene_id):
        return 'Entrez_Gene_ID'
    
    # UniProt ID (6字符格式)
    if re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$', gene_id):
        return 'UniProt'
    
    # HGNC ID
    if re.match(r'^HGNC:\d+$', gene_id, re.IGNORECASE):
        return 'HGNC'
    
    # MGI ID (小鼠)
    if re.match(r'^MGI:\d+$', gene_id):
        return 'MGI'
    
    # RGD ID (大鼠)
    if re.match(r'^RGD:\d+$', gene_id):
        return 'RGD'
    
    # FlyBase ID (果蝇)
    if re.match(r'^FBgn\d+$', gene_id):
        return 'FlyBase'
    
    # WormBase ID (线虫)
    if re.match(r'^WBGene\d+$', gene_id):
        return 'WormBase'
    
    # NCBI Accession (如 BC032353, AK129345)
    if re.match(r'^[A-Z]{2}\d{6}$', gene_id):
        return 'NCBI_Accession'
    
    # Gene Symbol: 短名称，通常1-15个字符，字母开头，可含数字和连字符
    # 典型格式: CD4, TP53, IL6, BRCA1, HLA-A, etc.
    if re.match(r'^[A-Za-z][A-Za-z0-9\-\.]{0,14}$', gene_id):
        # 进一步验证是否像 Gene Symbol (通常大写或首字母大写)
        if gene_id.isupper() or gene_id[0].isupper():
            return 'Gene_Symbol'
    
    return 'Unknown'


def detect_gene_id_type_from_list(gene_list: List[str], sample_size: int = 50) -> Dict:
    """
    根据一批基因ID判断整体类型
    
    Args:
        gene_list: 基因ID列表
        sample_size: 用于检测的样本数量
        
    Returns:
        包含检测结果的字典
    """
    # 取样本
    sample = gene_list[:sample_size] if len(gene_list) > sample_size else gene_list
    
    types = [identify_gene_id_type(g) for g in sample]
    type_counts = Counter(types)
    
    # 返回最常见的类型
    most_common = type_counts.most_common(1)[0]
    confidence = most_common[1] / len(sample)
    
    return {
        'detected_type': most_common[0],
        'confidence': round(confidence, 3),
        'sample_size': len(sample),
        'type_distribution': dict(type_counts)
    }


# ========== h5ad 文件 var 属性分析 ==========

def analyze_var_attributes(adata) -> Dict:
    """
    分析 AnnData 对象 var 中所有属性的基因ID类型
    
    Args:
        adata: AnnData 对象
        
    Returns:
        包含各属性分析结果的字典
    """
    results = {
        'var_names': None,
        'columns': {},
        'gene_symbol_candidates': []
    }
    
    # 分析 var_names
    var_names_list = adata.var_names.tolist()
    var_names_analysis = detect_gene_id_type_from_list(var_names_list)
    results['var_names'] = {
        'type': var_names_analysis['detected_type'],
        'confidence': var_names_analysis['confidence'],
        'is_gene_symbol': var_names_analysis['detected_type'] == 'Gene_Symbol',
        'sample': var_names_list[:5]
    }
    
    # 分析 var 中的每一列
    for col in adata.var.columns:
        col_values = adata.var[col].dropna().astype(str).tolist()
        if len(col_values) == 0:
            results['columns'][col] = {'type': 'Empty', 'confidence': 0}
            continue
            
        col_analysis = detect_gene_id_type_from_list(col_values)
        results['columns'][col] = {
            'type': col_analysis['detected_type'],
            'confidence': col_analysis['confidence'],
            'is_gene_symbol': col_analysis['detected_type'] == 'Gene_Symbol',
            'sample': col_values[:3]
        }
        
        # 记录可能是 Gene Symbol 的列
        if col_analysis['detected_type'] == 'Gene_Symbol' and col_analysis['confidence'] > 0.7:
            results['gene_symbol_candidates'].append({
                'column': col,
                'confidence': col_analysis['confidence']
            })
    
    return results


def find_best_gene_symbol_column(adata) -> Optional[str]:
    """
    找到最佳的 Gene Symbol 列名
    
    Args:
        adata: AnnData 对象
        
    Returns:
        Gene Symbol 列名，如果没有找到返回 None
    """
    # 常见的 Gene Symbol 列名
    common_symbol_names = [
        'gene_short_name', 'gene_symbol', 'gene_name', 'symbol', 
        'Gene', 'gene', 'SYMBOL', 'Symbol', 'name', 'Name',
        'external_gene_name', 'gene_symbols', 'genes'
    ]
    
    # 首先检查常见列名
    for col_name in common_symbol_names:
        if col_name in adata.var.columns:
            col_values = adata.var[col_name].dropna().astype(str).tolist()
            if len(col_values) > 0:
                analysis = detect_gene_id_type_from_list(col_values)
                if analysis['detected_type'] == 'Gene_Symbol' and analysis['confidence'] > 0.6:
                    return col_name
    
    # 如果常见列名中没有找到，遍历所有列
    analysis_results = analyze_var_attributes(adata)
    candidates = analysis_results.get('gene_symbol_candidates', [])
    
    if candidates:
        # 按置信度排序，返回最高的
        candidates.sort(key=lambda x: x['confidence'], reverse=True)
        return candidates[0]['column']
    
    return None


# ========== 主要功能：转换 var_names ==========

def convert_var_names_to_symbol(adata, inplace: bool = True, verbose: bool = True, save_path: str = None):
    """
    将 AnnData 的 var_names 转换为 Gene Symbol
    
    如果 var_names 已经是 Gene Symbol，或者找不到 Gene Symbol 列，则不做任何修改
    
    Args:
        adata: AnnData 对象或 h5ad 文件路径
        inplace: 是否原地修改（仅当 adata 是对象时有效）
        verbose: 是否打印详细信息
        save_path: 保存路径，如果为 None 且需要转换，则自动生成 _symbol_converted 后缀
        
    Returns:
        如果输入是文件路径，返回新文件路径（如果进行了转换）或原路径（如果未转换）
        如果输入是对象且 inplace=False，返回修改后的 AnnData 对象
        如果输入是对象且 inplace=True，返回 (converted: bool, save_path: str)
    """
    # 处理文件路径输入
    if isinstance(adata, str):
        file_path = adata
        if verbose:
            print(f"[INFO] 正在加载文件: {file_path}")
        adata = sc.read_h5ad(file_path)
        is_file_input = True
        if save_path is None:
            save_path = file_path.replace('.h5ad', '_symbol_converted.h5ad')
    else:
        is_file_input = False
        if not inplace:
            adata = adata.copy()
    
    # 1. 分析当前 var_names 的类型
    var_names_list = adata.var_names.tolist()
    current_type = detect_gene_id_type_from_list(var_names_list)
    
    if verbose:
        print(f"[INFO] 当前 var_names 类型: {current_type['detected_type']} (置信度: {current_type['confidence']:.1%})")
        print(f"[INFO] 示例: {var_names_list[:5]}")
    
    # 2. 如果已经是 Gene Symbol，不需要转换
    if current_type['detected_type'] == 'Gene_Symbol' and current_type['confidence'] > 0.7:
        if verbose:
            print("[INFO] var_names 已经是 Gene Symbol，无需转换")
        if is_file_input:
            return file_path
        else:
            return (False, None) if inplace else adata
    
    # 3. 查找 Gene Symbol 列
    symbol_col = find_best_gene_symbol_column(adata)
    
    if symbol_col is None:
        if verbose:
            print("[WARN] 未找到 Gene Symbol 列，无法进行转换")
            print(f"[INFO] var 中可用的列: {list(adata.var.columns)}")
        if is_file_input:
            return file_path
        else:
            return (False, None) if inplace else adata
    
    if verbose:
        print(f"[INFO] 找到 Gene Symbol 列: '{symbol_col}'")
        sample_symbols = adata.var[symbol_col].dropna().astype(str).tolist()[:5]
        print(f"[INFO] 示例值: {sample_symbols}")
    
    # 4. 执行转换
    # 保存原始 var_names
    adata.var['original_var_names'] = adata.var_names.tolist()
    
    # 获取新的 var_names，处理空值和空字符串
    original_names = adata.var_names.tolist()
    symbol_values = adata.var[symbol_col].tolist()
    
    new_var_names = []
    for orig, sym in zip(original_names, symbol_values):
        # 如果 symbol 有效（非空、非NaN、非空字符串），使用 symbol
        if pd.notna(sym) and str(sym).strip() != '':
            new_var_names.append(str(sym))
        else:
            new_var_names.append(str(orig))
    
    # 设置新的 var_names
    adata.var_names = pd.Index(new_var_names)
    
    # 处理重复名称
    adata.var_names_make_unique()
    
    # 在 var_names 去重后的最终名称基础上，构建原始ID到新ID的映射
    final_var_names = adata.var_names.tolist()
    id_mapping = {orig: new for orig, new in zip(original_names, final_var_names)}
    
    # 同步更新 raw.var_names（如果存在）
    if getattr(adata, "raw", None) is not None:
        try:
            # 获取 raw 的 var DataFrame 和表达矩阵
            raw_var = adata.raw.var.copy()
            raw_X = adata.raw.X
            raw_obs_names = adata.raw.obs_names
            
            # 用映射表把 raw 中的 var_names 替换成对应的 symbol
            raw_var_names = raw_var.index.tolist()
            new_raw_var_names = [id_mapping.get(name, name) for name in raw_var_names]
            raw_var.index = pd.Index(new_raw_var_names)
            
            # 创建临时 AnnData 对象并重新赋值给 raw
            from anndata import AnnData
            temp_adata = AnnData(X=raw_X, obs=pd.DataFrame(index=raw_obs_names), var=raw_var)
            adata.raw = temp_adata
            
            if verbose:
                print("[INFO] 已同步更新 adata.raw.var_names 为 Gene Symbol")
        except Exception as e:
            if verbose:
                print(f"[WARN] 更新 adata.raw.var_names 失败: {e}")
    
    if verbose:
        print(f"[SUCCESS] 已将 var_names 从 {current_type['detected_type']} 转换为 Gene Symbol")
        print(f"[INFO] 原始 ID 已保存在 adata.var['original_var_names']")
    
    # 5. 保存文件（如果需要）
    if is_file_input or (save_path and inplace):
        adata.write_h5ad(save_path)
        if verbose:
            print(f"[SUCCESS] 已保存转换后的文件: {save_path}")
        return save_path
    
    if inplace:
        return (True, save_path)
    else:
        return adata


def analyze_and_report(h5ad_path: str) -> Dict:
    """
    分析 h5ad 文件并生成完整报告
    
    Args:
        h5ad_path: h5ad 文件路径
        
    Returns:
        分析报告字典
    """
    print(f"[INFO] 正在加载文件: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"[INFO] 数据维度: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"[INFO] var 属性列: {list(adata.var.columns)}")
    print()
    
    # 分析所有属性
    analysis = analyze_var_attributes(adata)
    
    # 打印报告
    print("=" * 60)
    print("基因ID类型分析报告")
    print("=" * 60)
    
    print(f"\n【var_names 分析】")
    vn = analysis['var_names']
    print(f"  类型: {vn['type']}")
    print(f"  置信度: {vn['confidence']:.1%}")
    print(f"  是否为 Gene Symbol: {'是' if vn['is_gene_symbol'] else '否'}")
    print(f"  示例: {vn['sample']}")
    
    print(f"\n【var 列属性分析】")
    for col, info in analysis['columns'].items():
        symbol_tag = " ★ Gene Symbol" if info['is_gene_symbol'] else ""
        print(f"  {col}:")
        print(f"    类型: {info['type']}{symbol_tag}")
        print(f"    置信度: {info['confidence']:.1%}")
        print(f"    示例: {info.get('sample', [])[:3]}")
    
    print(f"\n【Gene Symbol 候选列】")
    if analysis['gene_symbol_candidates']:
        for candidate in analysis['gene_symbol_candidates']:
            print(f"  - {candidate['column']} (置信度: {candidate['confidence']:.1%})")
    else:
        print("  未找到 Gene Symbol 候选列")
    
    print("=" * 60)
    
    return analysis


# ========== CLI 入口 ==========

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("用法: python gene_id_converter.py <h5ad_file> [--convert]")
        print("  --convert: 自动转换 var_names 为 Gene Symbol 并保存")
        sys.exit(1)
    
    h5ad_path = sys.argv[1]
    do_convert = '--convert' in sys.argv
    
    # 分析文件
    analyze_and_report(h5ad_path)
    
    # 如果指定了 --convert，执行转换
    if do_convert:
        print("\n[INFO] 正在执行转换...")
        adata = sc.read_h5ad(h5ad_path)
        converted = convert_var_names_to_symbol(adata, inplace=True, verbose=True)
        
        if converted:
            output_path = h5ad_path.replace('.h5ad', '_symbol_converted.h5ad')
            adata.write_h5ad(output_path)
            print(f"[SUCCESS] 已保存转换后的文件: {output_path}")
        else:
            print("[INFO] 未进行转换")
