# 原生python
from collections.abc import Sequence
from typing import List, Optional, Union

# 第三方包
import pandas as pd
import networkx as nx

# 目录下包
from ._settings import get_driver

# ===============================
# 混合检索：名称标准化模块
# ===============================

# 图谱类型 -> 混合检索类型 映射表
_TYPE_TO_HYBRID = {
    "Gene": "Gene|Protein",
    "Protein": "Gene|Protein",
    "Gene|Protein": "Gene|Protein",
    "Cell": "Cell|Tissue",
    "Tissue": "Cell|Tissue",
    "Cell|Tissue": "Cell|Tissue",
    "Species": "Species",
    "Disease": "Disease|Phenotype",
    "Phenotype": "Disease|Phenotype",
    "Disease|Phenotype": "Disease|Phenotype",
    "Chemical": "Chemical",
    "Process": "Process|Function|Pathway|Cell_Component",
    "Function": "Process|Function|Pathway|Cell_Component",
    "Pathway": "Process|Function|Pathway|Cell_Component",
    "Cell_Component": "Process|Function|Pathway|Cell_Component",
    "Mutation": "Mutation",
}

_search_client = None


def _get_search_client():
    """懒加载 SearchClient 单例"""
    global _search_client
    if _search_client is None:
        from .SearchClient import BRICKSearchClient
        _search_client = BRICKSearchClient()
    return _search_client


def standardize_names(
    names: Union[str, List[str]],
    entity_type: Optional[str] = None,
    min_score: float = 0.0,
) -> List[str]:
    """
    使用混合检索将用户输入的名字标准化为图谱中的 primary_name。
    
    Args:
        names: 用户输入的名字（可以是中文，如"白细胞"）
        entity_type: 图谱实体类型，如 'Cell', 'Gene', 'Species' 等
        min_score: 分数阈值，低于此值则保留原名（默认 0.0 不过滤）
    
    Returns:
        标准化后的名字列表，查不到则保留原名
    
    Example:
        >>> standardize_names(["白细胞", "红细胞"], entity_type="Cell")
        ["WBC", "RBC"]  # 假设图谱中的标准名
    """
    # 标准化为 list
    if isinstance(names, str):
        names = [names]
    else:
        names = list(names)
    
    if not names:
        return names
    
    # 类型不在映射表里，直接原样返回
    if entity_type is None:
        return names
    
    hybrid_type = _TYPE_TO_HYBRID.get(entity_type)
    if hybrid_type is None:
        return names
    
    try:
        client = _get_search_client()
    except Exception:
        # SearchClient 初始化失败，不影响原有功能
        return names
    
    # ========== 新增：对 Gene 类型先做 Ensembl/NCBI ID 预转换 ==========
    if hybrid_type == "Gene|Protein":
        names = [client.normalize_gene_query(name) for name in names]
    # ================================================================
    
    standardized = []
    for name in names:
        payload = {
            "query_id": f"normalize-{entity_type}",
            "options": {"top_k": 1, "return_diagnostics": True},
            hybrid_type: [name],
        }
        try:
            std_df, diag_df, resp = client.search_hybrid(
                payload,
                type_mix_override={hybrid_type: 1.0},
                auto_normalize_gene=False,  # 已经预处理过，不需要重复转换
            )
            # 没命中，保留原名
            if diag_df is None or len(diag_df) == 0:
                standardized.append(name)
                continue
            
            best = diag_df.iloc[0]
            score = float(best.get("final_score", 0.0))
            if score < min_score:
                # 分数太低，认为不可靠，保留原名
                standardized.append(name)
            else:
                standardized.append(best["primary_name"])
        except Exception:
            # 搜索出问题时，不影响图谱查询功能
            standardized.append(name)
    
    return standardized


def _standardize_entity_set(
    entity_set: List[str],
    entity_type: Optional[str],
    query_attribution: str,
) -> List[str]:
    """
    内部函数：只在 query_attribution == 'name' 时做标准化
    """
    if query_attribution == "name" and entity_type is not None:
        return standardize_names(entity_set, entity_type=entity_type)
    return entity_set


def _flatten_dict(nested_dict, parent_key='', sep='.'):
    """ 
    Recursively flattens a nested dictionary (including lists) into a single-layer dictionary, 
    where keys represent the nested paths.

    Args:
        nested_dict (dict): The nested dictionary, which may contain lists or other nested dictionaries.
        parent_key (str, optional): The parent key for the current recursion level, defaults to an empty string.
        sep (str, optional): The separator used to connect the nested keys, defaults to '.'.

    Returns:
        dict: A flattened dictionary with keys as the path of the nested structure and values as the corresponding elements.
    """
    items = []
    for key, value in nested_dict.items():
        new_key = f"{parent_key}{sep}{key}" if parent_key else key
        if isinstance(value, dict):
            items.extend(_flatten_dict(value, new_key, sep).items())
        elif isinstance(value, list):
            for i, sub_value in enumerate(value):
                if isinstance(sub_value, (dict, list)):
                    sub_items = _flatten_dict(sub_value, f"{new_key}{sep}{i}", sep)
                    items.extend(sub_items.items())
                else:
                    items.append((f"{new_key}{sep}{i}", sub_value))
        else:
            items.append((new_key, value))
    return dict(items)

""" def _neo4j2networkx(neo4j_g):
    edges = []
    name2id = {}
    for n in neo4j_g.relationships:
        x, y = n.nodes
        name2id[x['name']] = x['id']
        name2id[y['name']] = y['id']
        x = x['name']
        y = y['name']
        if x is not None and y is not None:
            r = n.type
            rela_conf = n['relation_confidence']
            #if rela_conf is not None and rela_conf<0:
            #    c = 'red'
            #else:
            #    c = 'green'
                
            edges.append([x,y, {'type':r, 'relation_confidence':rela_conf, 'color':'black'}])
    nxg = nx.MultiDiGraph()
    nxg.add_edges_from(edges)
    return nxg
 """

def _neo4j2networkx(neo4j_g):
    """ 
    Converts a Neo4j graph object to a NetworkX MultiDiGraph object.

    Args:
        neo4j_g (neo4j.Graph): A Neo4j graph object that contains nodes and relationships.

    Returns:
        networkx.MultiDiGraph: A NetworkX MultiDiGraph object with nodes and edges from the Neo4j graph.
    """
    nodes = []
    edges = []
    
    for x in neo4j_g.nodes:
        n = x._properties
        nodes.append([n['id'], n])
    
    for y in neo4j_g.relationships:
        i = y.start_node._properties['id']
        j = y.end_node._properties['id']
        edges.append([i, j, y._properties])
    
    nxg = nx.MultiDiGraph()
    nxg.add_nodes_from(nodes)
    nxg.add_edges_from(edges)
    return nxg

def _neo4j2networkx(neo4j_g):
    """
    Converts a Neo4j graph into a NetworkX MultiDiGraph.

    Args:
        neo4j_g (neo4j.Graph): A Neo4j graph object obtained from a database query, containing nodes and relationships.

    Returns:
        networkx.MultiDiGraph: A directed multigraph representation of the Neo4j graph, where:
            - Nodes are identified by their `id` and include a dictionary of Neo4j attributes.
            - Edges connect nodes via their `id` and contain a dictionary of relationship attributes.
    """    
    if not hasattr(neo4j_g, "nodes") or not hasattr(neo4j_g, "relationships"):
        raise AttributeError("Your input neo4j_g must contains `nodes` and `relationships` attributes.")

    nodes = []
    edges = []
    
    for x in neo4j_g.nodes:
        n = x._properties
        node_id = n.get("id") 
        if node_id is None:
            raise KeyError(f"Node {n} is missing key `id`")
        nodes.append([node_id, n])
    
    for y in neo4j_g.relationships:
        start_node_id = y.start_node._properties.get("id")
        end_node_id = y.end_node._properties.get("id")

        if start_node_id is None or end_node_id is None:
            raise KeyError(f"Node in relationship {y._properties} is missing key `id`.")
        
        edges.append([start_node_id, end_node_id, y._properties])
    
    nxg = nx.MultiDiGraph()
    nxg.add_nodes_from(nodes)
    nxg.add_edges_from(edges)
    return nxg

def _stardard_list(entity_set):
    """ 
    Converts the input entity set (`entity_set`) into a standardized list. 
    Supports both strings and iterable objects (e.g., lists, tuples).

    Args:
        entity_set (str or iterable): The entity set to be converted. It can be a string or any iterable object such as a list or tuple.

    Returns:
        list: A standardized list. If the input is a string, it is wrapped as a single-element list. 
              If the input is another iterable, it is converted into a list.
    """
    assert len(entity_set) > 0, 'entity_set must be iterrowable and length > 0'
    if isinstance(entity_set, str):
        entity_set = [entity_set]
    else:
        entity_set = list(entity_set)
    return entity_set    

def _stardard_type_str(types=None):
    """ 
    Generates a standardized type string based on the given input (`types`). 
    Supports a single string or a list of strings.

    Args:
        types (str or list, optional): The type(s) to process. It can be a single string or a list of strings. 
                                       If `None`, an empty string is returned.

    Returns:
        str: A standardized type string:
            - If the input is a string, returns `:string`.
            - If the input is a list of strings, returns `:type1|type2|...`.
            - If `types` is `None`, returns an empty string.
    """
    if types is None:
        return ""
    elif isinstance(types, str) :
        return f":{types}"
    elif isinstance(types, list):
        return f":{'|'.join(types)}"
    else:
        raise ValueError('type must be a string or list.')
        

def query_cypher(cypher, parameters=None, return_type='dataframe'):
    """
    Executes a Cypher query and returns the query result. The output format is determined by the `return_type` parameter, 
    which supports returning a Pandas DataFrame, a NetworkX MultiDiGraph, or a list.

    Args:
        cypher (str): The Cypher query string to be executed.
        parameters (dict or None, optional): Query parameters (default is None). If needed, parameters can be passed to the query.
        return_type (str): Specifies the format of the returned object. Available options:
            - 'dataframe': Returns the result as a Pandas DataFrame.
            - 'graph': Returns the result as a NetworkX MultiDiGraph.
            - 'list': Returns the result as a list of records.

    Returns:
        - If `return_type` is 'dataframe', returns a Pandas DataFrame containing the query results.
        - If `return_type` is 'graph', returns a NetworkX MultiDiGraph containing nodes and edges.
        - If `return_type` is 'list', returns a list of records.
    """
    driver = get_driver()

    if return_type not in ['dataframe', 'graph', 'list']:
        raise ValueError("return_type must be chosen from ['dataframe', 'graph']")

    if return_type == 'dataframe':
        with driver.session() as session:
            #print(driver)
            result = session.run(cypher, parameters)
            original_result = [record for record in result]
        
        result = [record.data() for record in original_result]
        result = [_flatten_dict(x) for x in result]
        result = pd.DataFrame(result)

        if 'path.1' in result:
            edges_properties = []
            for record in original_result:
                tmp_dict = {}
                for i, j in enumerate([i._properties for i in record['path'].relationships]):
                    for x, y in j.items():
                        tmp_dict[f"path.{2*i+1}.{x}"] = y
                edges_properties.append(tmp_dict)
            edges_properties = pd.DataFrame(edges_properties)
    
            result = pd.concat([result, edges_properties], axis = 1)
            result = result[sorted(result.columns)]

    elif return_type =='graph':
        with driver.session() as session:
            #print(driver)
            result = session.run(cypher, parameters)
            result = result.graph()
        assert len(result.nodes) > 0, "The cypher should return nodes, edges or paths.\
         i.e. 'MATCH path=(n)-[r]-(m) RETURN path LIMIT 10'"
        result = _neo4j2networkx(result)
    
    elif return_type == 'list':
        with driver.session() as session:
            #print(driver)
            result = session.run(cypher, parameters)
            result = [record for record in result]

    return result


def query_node(source_entity_set, source_entity_type=None, query_attribution=None, limit_number=None):
    """
    Queries node information from the graph database based on specified conditions.

    Args:
        source_entity_set (list or set): The set of entities to query, such as gene names or cell names. 
            If a set is provided, it will be converted into a list internally.
        source_entity_type (str, optional): The entity type, e.g., 'Cell' or 'Gene'. Processing rules:
            - If 'Cell', the first letter of each entity in the set is capitalized.
            - If 'Gene', ensure proper capitalization:
                - For human genes, names should be in all uppercase (e.g., 'FBXO34').
                - For mouse genes, names should have only the first letter capitalized (e.g., 'Fbxo34').
        query_attribution (str, optional): The attribute used for matching nodes. Must be either 'name' or 'id'.
        limit_number (int, optional): The maximum number of results to return. If None, there is no limit; 
            if an integer is provided, the number of returned results is restricted accordingly.

    Returns:
        pandas.DataFrame: A DataFrame containing information about nodes that meet the specified conditions.
    """
    source_entity_type_str = _stardard_type_str(source_entity_type)
 
    if isinstance(source_entity_set,set):
        source_entity_set = list(source_entity_set)
    source_entity_set = _stardard_list(source_entity_set)
    
    # 名称标准化：将用户输入转为图谱标准名
    source_entity_set = _standardize_entity_set(source_entity_set, source_entity_type, query_attribution)

    if query_attribution not in ["name","id"]:
        raise ValueError("query_attribution must be 'name' or 'id'.")
    
    if limit_number is None:
        limit_str = ""
    elif isinstance(limit_number, int):
        limit_str = f"LIMIT {limit_number}"
    else:
        raise ValueError('limit_number must be an int.')
    
    cypher = f"""UNWIND $source_entity_set AS entity MATCH (n{source_entity_type_str}) WHERE entity IN COALESCE([n.{query_attribution}] + SPLIT(n.synonym, '|'), []) RETURN n {limit_str}"""

    driver = get_driver()
    parameters = {'source_entity_set':source_entity_set}

    with driver.session() as session:
        result = session.run(cypher, parameters)
        original_result = [record for record in result]
        result = [record.data() for record in original_result]
        result = [_flatten_dict(x) for x in result]
        result = pd.DataFrame(result)
        result = result[sorted(result.columns)]

    return result

def query_path(source_entity_set, target_entity_set=None, relation=None, source_entity_type=None, target_entity_type=None,  
               limit_number=None, multi_hop=1, query_attribution=None, relation_info_source=None, directed=False, 
               link=2, return_type='dataframe'):
    """
    Queries path information from a Neo4j graph database, supporting both single-hop and multi-hop queries 
    with multiple return formats.

    Args:
        source_entity_set (list or set): The set of source entities, such as gene names or cell names. 
            Must be an iterable and cannot be empty.
        target_entity_set (list or set, optional): The set of target entities. Defaults to None, meaning 
            all paths related to the source entities will be queried.
        relation (list of str, optional): A list of relationship types. The length of this list must match 
            the `link` parameter. If None, all relationships are queried.
        source_entity_type (str, optional): The type of the source entity (e.g., 'Cell' or 'Gene'), used to 
            restrict the search scope.
        target_entity_type (str or list of str, optional): The type of the target entity. Can be a single 
            string or a list (matching the length of `link`), used to restrict the endpoint types of paths.
        limit_number (int, optional): The maximum number of paths to return. Defaults to None (no limit).
        multi_hop (int or tuple(int, int), optional): The allowed number of hops. Can be an integer (maximum 
            hops) or a tuple (min, max) to define a range.
        query_attribution (str, optional): The attribute field for matching entities. Must be either 'name' 
            or 'id'.
        relation_info_source (str or list of str, optional): The source of relationship information, which 
            can be a string or a list used to filter specific relationship sources.
        directed (bool, default=False): Whether to query directed paths. Defaults to False (undirected paths).
        link (int, default=2): The number of steps in the queried path, i.e., the allowed number of edges.
        return_type (str, default='dataframe'): The format of the returned data. Available options:
            - 'dataframe': Returns the data as a Pandas DataFrame.
            - 'graph': Returns the data as a NetworkX MultiDiGraph.
            - 'list': Returns the raw query results as a list.

    Returns:
        pandas.DataFrame or networkx.MultiDiGraph or list:
            - If `return_type='dataframe'`, returns a Pandas DataFrame containing the queried path information.
            - If `return_type='graph'`, returns a NetworkX MultiDiGraph representing the graph structure.
            - If `return_type='list'`, returns the raw query result as a list.
    """
    source_entity_set = _stardard_list(source_entity_set)
    
    # 名称标准化：将用户输入转为图谱标准名
    source_entity_set = _standardize_entity_set(source_entity_set, source_entity_type, query_attribution)

    if target_entity_set is None:
        target_entity_set = ""
    else:
        target_entity_set = _stardard_list(target_entity_set)
        # 对 target 也做标准化（取 target_entity_type 的第一个类型）
        _target_type = target_entity_type[0] if isinstance(target_entity_type, list) else target_entity_type
        target_entity_set = _standardize_entity_set(target_entity_set, _target_type, query_attribution)

    if relation is None:
        relation_lst = ["" for _ in range(link)]
    elif isinstance(relation, list):
        if len(relation) != link:
            raise ValueError(f"The length of relation type: {len(relation)} does not match the path link length: {link}.")
        relation_lst = relation
    else:
        raise ValueError("relation must be a list or None.")

    source_entity_type_str = _stardard_type_str(source_entity_type)

    if target_entity_type is None:
        target_entity_type_str = ["" for _ in range(link)]
    elif isinstance(target_entity_type, list):
        if len(target_entity_type) != link:
            raise ValueError(f"The length of target entity type: {len(target_entity_type)} does not match the path link length: {link}.")
        target_entity_type_str = target_entity_type
    elif isinstance(target_entity_type, str):
        target_entity_type_str = [target_entity_type for _ in range(link)]
    else:
        raise ValueError("target_entity_type must be a string or list.")

    if limit_number is None:
        limit_str = ""
    elif isinstance(limit_number, int):
        limit_str = f"LIMIT {limit_number}"
    else:
        raise ValueError('limit_number must be an int.')

    if isinstance(multi_hop, int):
        multi_hop_str = f"*..{multi_hop}"
    elif isinstance(multi_hop, Sequence):
        a, b = multi_hop
        multi_hop_str = f"*{a}..{b}"
    else:
        raise ValueError('multi_hop must be an int or a tuple of int with length 2.')

    if relation_info_source is None:
        relation_info_source_str = ''
    elif isinstance(relation_info_source, str):
        relation_info_source_str = f'AND apoc.text.join(r.info_source, ", ") CONTAINS "{relation_info_source}"'
    elif isinstance(relation_info_source, list):
        relation_info_source_str = f'AND ANY(substring IN {relation_info_source} WHERE apoc.text.join(r.info_source, ", ") CONTAINS substring)'
    else:
        raise ValueError('relation_info_source must be a string or a list of strings.')

    directed_str = '>' if directed else ''

    inner_text = ""

    if relation_lst and relation_lst[0]:
        begin_text = f"MATCH path=(n{source_entity_type_str})-[r0:{relation_lst[0]}{multi_hop_str}]-{directed_str}"
        for i in range(link - 1):
            inner_text += f"(m{i}:{target_entity_type_str[i]})-[r{i+1}:{relation_lst[i+1]}{multi_hop_str}]-{directed_str}"
        inner_text += f"(m{link-1}:{target_entity_type_str[-1]})"
    else:
        begin_text = f"MATCH path=(n{source_entity_type_str})-[r0{multi_hop_str}]-{directed_str}"
        for i in range(link - 1):
            inner_text += f"(m{i}:{target_entity_type_str[i]})-[r{i+1}{multi_hop_str}]-{directed_str}"
        inner_text += f"(m{link-1}:{target_entity_type_str[-1]})"

    condition_text = f" WHERE n.{query_attribution} IN $source_entity_set AND m{link-1}.{query_attribution} IN $target_entity_set {relation_info_source_str}"

    return_text = f" RETURN path {limit_str}"

    cypher = begin_text + inner_text + condition_text + return_text
    parameters = {'source_entity_set': source_entity_set, 'target_entity_set': target_entity_set}

    return query_cypher(cypher, parameters, return_type)


def query_shortest_path(source_entity_name, target_entity_name, relation=None, multi_hop=1, 
                        directed=False, source_entity_type=None, target_entity_type=None, relation_info_source=None, 
                        query_attribution='name', limit_number=None, return_type='dataframe'):
    """
    Queries the shortest path between two entities in a Neo4j graph database.

    Args:
        source_entity_name (str): The name of the source entity.
        target_entity_name (str or list of str): The name of the target entity.
        relation (str or list of str, optional): The relationship type(s) to filter paths. If None, all relationships are considered.
        multi_hop (int or tuple(int, int), optional): The allowed number of hops in the path. Can be:
            - An integer specifying the maximum number of hops.
            - A tuple (min_hops, max_hops) defining a range.
        directed (bool, default=False): Whether to query directed paths. Defaults to False (undirected paths).
        source_entity_type (str, optional): The type of the source entity, such as 'Cell' or 'Gene'.
        target_entity_type (str, optional): The type of the target entity. Defaults to the same type as `source_entity_type` if not specified.
        relation_info_source (str or list of str, optional): The source(s) of relationship information for filtering specific relationship sources.
        query_attribution (str, default='name'): The attribute field for matching entities. Must be either 'name' or 'id'.
        limit_number (int, optional): The maximum number of paths to return. If None, there is no limit.
        return_type (str, default='dataframe'): The format of the returned data. Available options:
            - 'dataframe': Returns a Pandas DataFrame.
            - 'graph': Returns a NetworkX MultiDiGraph.
            - 'list': Returns the raw query results as a list.

    Returns:
        pandas.DataFrame or networkx.MultiDiGraph or list:
            - If `return_type='dataframe'`, returns a Pandas DataFrame containing the queried shortest paths.
            - If `return_type='graph'`, returns a NetworkX MultiDiGraph representing the shortest path(s).
            - If `return_type='list'`, returns the raw query result as a list.

    Raises:
        ValueError: If `multi_hop` is neither an integer nor a tuple of two integers.
        ValueError: If `relation_info_source` is neither a string nor a list of strings.
        ValueError: If `limit_number` is not an integer.
    """
    source_entity_type_str = _stardard_type_str(source_entity_type)

    if target_entity_type is None:
        target_entity_type_str = source_entity_type_str
    else:
        target_entity_type_str = _stardard_type_str(target_entity_type)
    
    # 名称标准化：将用户输入转为图谱标准名
    if query_attribution == "name":
        if source_entity_name:
            result = standardize_names([source_entity_name], entity_type=source_entity_type)
            source_entity_name = result[0] if result else source_entity_name
        
        _target_type = target_entity_type or source_entity_type
        if isinstance(target_entity_name, list):
            target_entity_name = standardize_names(target_entity_name, entity_type=_target_type)
        elif isinstance(target_entity_name, str) and target_entity_name:
            result = standardize_names([target_entity_name], entity_type=_target_type)
            target_entity_name = result[0] if result else target_entity_name

    relation_str = _stardard_type_str(relation)
    
    if isinstance(multi_hop, int):
        multi_hop_str = f"*..{multi_hop}"
    elif isinstance(multi_hop, Sequence):
        a, b = multi_hop
        multi_hop_str = f"*{a}..{b}"
    else:
        raise ValueError('multi_hop must be a int or a tuple of int with length of 2.')

    if directed:
        directed_str='>'
    else:
        directed_str=''

    if relation_info_source is None:
        relation_info_source_str = ''
    elif isinstance(relation_info_source, str):
        relation_info_source_str = f'AND apoc.text.join(r.info_source, ", ") CONTAINS "{relation_info_source}"'
    elif isinstance(relation_info_source, list):
        relation_info_source_str = f'AND ANY(substring IN {relation_info_source} WHERE apoc.text.join(r.info_source, ", ") CONTAINS substring)'
    else:
        raise ValueError('relation_info_source must be a string or a list of string.')
    
    if limit_number is None:
        limit_str = ""
    elif isinstance(limit_number, int):
        limit_str = f"LIMIT {limit_number}"
    else:
        raise ValueError('limit_number must be an int.')

    if type(target_entity_name) == list:
        cypher = f"MATCH path=shortestpath((n{source_entity_type_str})-[r{relation_str}{multi_hop_str}]-{directed_str}(m{target_entity_type_str}))\
               WHERE n.{query_attribution} = $source_entity_name AND m.{query_attribution} IN $target_entity_name {relation_info_source_str}\
               RETURN path {limit_str}"
    if type(target_entity_name) == str:
        cypher = f"MATCH path=shortestpath((n{source_entity_type_str})-[r{relation_str}{multi_hop_str}]-{directed_str}(m{target_entity_type_str}))\
               WHERE n.{query_attribution} = $source_entity_name AND m.{query_attribution} = $target_entity_name {relation_info_source_str}\
               RETURN path {limit_str}"
    
    parameters = {'source_entity_name':source_entity_name, 'target_entity_name':target_entity_name}

    return query_cypher(cypher,parameters,return_type)

def query_relation(source_entity_set, target_entity_set=None, relation=None, multi_hop=1, directed=False,
                   source_entity_type=None, target_entity_type=None, relation_info_source=None, query_attribution='name',
                   limit_number=None, return_type='dataframe',):
    """
    Queries relationships between a set of entities, including direct and multi-hop connections.

    Args:
        source_entity_set (str or list): The starting entity or a list of entities to query.
            - If a single entity is provided as a string, it will be converted into a list internally.
            - This parameter is required.
        target_entity_set (str, list, or None, optional): The target entity or a list of target entities.
            - If None (default), the target entity set is considered the same as `source_entity_set`.
            - If a single entity is provided, it will be converted into a list internally.
        relation (str or list, optional): Specifies the relationship type(s) to filter, such as 'knows' or 'marker_of'.
            - If None (default), all relationship types are considered.
            - The function standardizes the format internally.
        multi_hop (int or tuple(int, int), optional): Specifies the allowed number of hops in the relationship path.
            - An integer N (default: 1) allows paths with up to N hops (expressed as `*..N` in Cypher).
            - A tuple (a, b) sets a range with at least `a` hops and at most `b` hops (`*a..b`).
        directed (bool, default=False): Whether the queried paths are directed.
            - False (default): Undirected paths.
            - True: Directed paths (n → m).
        source_entity_type (str, optional): Specifies the type of the source entity, such as 'Cell' or 'Person'.
            - If None, no type restriction is applied.
        target_entity_type (str, optional): Specifies the type of the target entity, such as 'Gene' or 'Company'.
            - If None, it defaults to the same type as `source_entity_type`.
        relation_info_source (str or list, optional): Filters relationships by their source information.
            - A string filters relationships that contain the specified source.
            - A list filters relationships that contain at least one of the specified sources.
            - If None (default), no filtering is applied.
        query_attribution (str, default='name'): Specifies the attribute field to match entities (e.g., 'name' or 'id').
        limit_number (int, optional): Limits the number of returned paths.
            - If None (default), there is no limit.
            - If an integer N is provided, at most N paths will be returned.
        return_type (str, default='dataframe'): Specifies the format of the returned data.
            - 'dataframe' (default): Returns a Pandas DataFrame.
            - 'graph': Returns a NetworkX MultiDiGraph.

    Returns:
        pandas.DataFrame or networkx.MultiDiGraph:
            - If `return_type='dataframe'`, returns a Pandas DataFrame containing the queried relationships.
            - If `return_type='graph'`, returns a NetworkX MultiDiGraph representing the relationships.

    Raises:
        ValueError: If `multi_hop` is neither an integer nor a tuple of two integers.
        ValueError: If `relation_info_source` is neither a string nor a list of strings.
        ValueError: If `limit_number` is not an integer.
    """
    source_entity_set = _stardard_list(source_entity_set)
    
    # 名称标准化：将用户输入转为图谱标准名
    source_entity_set = _standardize_entity_set(source_entity_set, source_entity_type, query_attribution)
    
    if target_entity_set is None :
        target_entity_set = source_entity_set
    else:
        target_entity_set = _stardard_list(target_entity_set)
        # 对 target 也做标准化
        _target_type = target_entity_type or source_entity_type
        target_entity_set = _standardize_entity_set(target_entity_set, _target_type, query_attribution)
    
    # relation-str
    relation_str = _stardard_type_str(relation)
    
    # multi-hop-str
    if isinstance(multi_hop, int):
        multi_hop_str = f"*..{multi_hop}"
    elif isinstance(multi_hop, Sequence):
        a, b = multi_hop
        multi_hop_str = f"*{a}..{b}"
    else:
        raise ValueError('multi_hop must be a int or a tuple of int with length of 2.')
    
    if directed:
        directed_str='>'
    else:
        directed_str=''
        
    source_entity_type_str = _stardard_type_str(source_entity_type)
    if target_entity_type is None:
        target_entity_type_str =source_entity_type_str
    else:
        target_entity_type_str = _stardard_type_str(target_entity_type)
    
    if relation_info_source is None:
        relation_info_source_str = ''
    elif isinstance(relation_info_source, str):
        relation_info_source_str = f'AND apoc.text.join(r.info_source, ", ") CONTAINS "{relation_info_source}"'
    elif isinstance(relation_info_source, list):
        relation_info_source_str = f'AND ANY(substring IN {relation_info_source} WHERE apoc.text.join(r.info_source, ", ") CONTAINS substring)'
    else:
        raise ValueError('relation_info_source must be a string or a list of string.')
    
    if limit_number is None:
        limit_str = ""
    elif isinstance(limit_number, int):
        limit_str = f"LIMIT {limit_number}"
    else:
        raise ValueError('limit_number must be an int.')
    
    #cypher_templete = "MATCH path=(n:Cell)-[r:marker_of*1..2]-(m:Gene) WHERE n.name IN $source_entity_set AND m.name IN $target_entity_set  RETURN path LIMIT 10"
    cypher_templete = f"MATCH path=(n{source_entity_type_str})-[r{relation_str}{multi_hop_str}]-{directed_str}(m{target_entity_type_str}) \
                        WHERE n.{query_attribution} IN $source_entity_set \
                        AND m.{query_attribution} IN $target_entity_set  \
                        {relation_info_source_str} \
                        RETURN path {limit_str}"

    parameters = {'source_entity_set':source_entity_set, 'target_entity_set':target_entity_set}
    return query_cypher(cypher_templete, parameters=parameters, return_type=return_type)


def query_neighbor(source_entity_set, relation=None, multi_hop=1, directed=False, source_entity_type=None,
                   target_entity_type=None, relation_info_source=None, query_attribution='name',
                   limit_number=None, return_type='dataframe',):
    """
    Queries the neighboring nodes of a set of source entities along with their relationships.

    This function retrieves neighboring nodes of specified source entities in a knowledge graph and returns 
    relationship information. Users can specify relationship types, hop counts, direction, and other query conditions.

    Args:
        source_entity_set (str or list): The source entity or a list of entities for querying neighboring nodes.
            - If a single entity is provided, it will be converted into a list internally.
            - This parameter is required.
        relation (str or list, optional): The relationship type(s) to query, such as 'knows' or 'marker_of'.
            - If None, all relationship types are considered.
            - The function standardizes the format internally.
        multi_hop (int or tuple(int, int), optional): The number of hops allowed in the relationship path.
            - An integer N (default: 1) allows paths with exactly N hops.
            - A tuple (a, b) allows paths with at least `a` and at most `b` hops.
        directed (bool, default=False): Specifies whether the queried paths should be directed.
            - False (default): Retrieves both incoming and outgoing relationships (undirected).
            - True: Retrieves only outgoing relationships (source entity → target entity).
        source_entity_type (str, optional): The type of the source entity, such as 'Cell' or 'Person'.
            - If None, no type restriction is applied.
        target_entity_type (str, optional): The type of the target entity, such as 'Gene' or 'Company'.
            - If None, no type restriction is applied.
        relation_info_source (str or list, optional): Filters relationships based on their `info_source` attribute.
            - A string filters relationships containing the specified source.
            - A list filters relationships containing any of the specified sources.
            - If None (default), no filtering is applied.
        query_attribution (str, default='name'): The attribute used to match entities from `source_entity_set`.
        limit_number (int, optional): Limits the number of returned results.
            - If None (default), there is no limit.
            - If an integer N is provided, at most N results will be returned.
        return_type (str, default='dataframe'): Specifies the format of the returned data.
            - 'dataframe' (default): Returns a Pandas DataFrame containing the queried relationship paths.
            - 'triplet': Only applicable when `multi_hop=1`. Returns a DataFrame in triplet format (source node, relationship, target node).

    Returns:
        pandas.DataFrame or networkx.MultiDiGraph:
            - If `return_type='dataframe'`, returns a Pandas DataFrame or a NetworkX MultiDiGraph with query results.
            - If `return_type='triplet'`, returns a Pandas DataFrame containing:
              (source node ID, source node name, source node type, relationship, target node ID, target node name, target node type).

    Raises:
        ValueError: If `multi_hop` is neither an integer nor a tuple of two integers.
        ValueError: If `relation_info_source` is neither a string nor a list of strings.
        ValueError: If `limit_number` is not an integer.
    """
    source_entity_set = _stardard_list(source_entity_set)
    
    # 名称标准化：将用户输入转为图谱标准名
    source_entity_set = _standardize_entity_set(source_entity_set, source_entity_type, query_attribution)
    
    # relation-str
    relation_str = _stardard_type_str(relation)
    
    # multi-hop-str
    if isinstance(multi_hop, int):
        multi_hop_str = f"*..{multi_hop}"
    elif isinstance(multi_hop, Sequence):
        a, b = multi_hop
        multi_hop_str = f"*{a}..{b}"
    else:
        raise ValueError('multi_hop must be a int or a tuple of int with length of 2.')
    
    if directed:
        directed_str='>'
    else:
        directed_str=''
        
    source_entity_type_str = _stardard_type_str(source_entity_type)
    target_entity_type_str = _stardard_type_str(target_entity_type)

    if relation_info_source is None:
        relation_info_source_str = ''
    elif isinstance(relation_info_source, str):
        relation_info_source_str = f'AND apoc.text.join(r.info_source, ", ") CONTAINS "{relation_info_source}"'
    elif isinstance(relation_info_source, list):
        relation_info_source_str = f'AND ANY(substring IN {relation_info_source} WHERE apoc.text.join(r.info_source, ", ") CONTAINS substring)'
    else:
        raise ValueError('relation_info_source must be a string or a list of string.')
    
    if limit_number is None:
        limit_str = ""
    elif isinstance(limit_number, int):
        limit_str = f"LIMIT {limit_number}"
    else:
        raise ValueError('limit_number must be an int.')

    if return_type == 'triplet':
        assert multi_hop == 1, "multi_hop must be 1 if return_type  == 'triplet' "
        # cypher_templete = "MATCH path=(n:Cell)-[r:marker_of*1..2]-(m:Gene) WHERE n.name IN $source_entity_set AND m.name IN $target_entity_set  RETURN path LIMIT 10"
        cypher_templete = f"MATCH path=(n{source_entity_type_str})-[r{relation_str}]-{directed_str}(m{target_entity_type_str}) \
                                WHERE n.{query_attribution} IN $source_entity_set \
                                {relation_info_source_str} \
                                RETURN n.id, n.name, n.type, r.relation, m.id, m.name, m.type {limit_str}"
        parameters = {'source_entity_set': source_entity_set}
        # return cypher_templete
        triplet_df = query_cypher(cypher_templete, parameters=parameters, return_type='dataframe')
        triplet_df.columns = ['path.0.id', 'path.0.name', 'path.0.type', 'path.1',
                              'path.2.id', 'path.2.name', 'path.2.type']
        return triplet_df
    else:
        #cypher_templete = "MATCH path=(n:Cell)-[r:marker_of*1..2]-(m:Gene) WHERE n.name IN $source_entity_set AND m.name IN $target_entity_set  RETURN path LIMIT 10"
        cypher_templete = f"MATCH path=(n{source_entity_type_str})-[r{relation_str}{multi_hop_str}]-{directed_str}(m{target_entity_type_str}) \
                            WHERE n.{query_attribution} IN $source_entity_set \
                            {relation_info_source_str} \
                            RETURN path {limit_str}"

        parameters = {'source_entity_set':source_entity_set}
        #return cypher_templete
        return query_cypher(cypher_templete, parameters=parameters, return_type=return_type)