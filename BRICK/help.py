import pandas as pd

from .querygraph import query_cypher, _stardard_type_str

def check_schema():
    labels = query_cypher("CALL db.labels() YIELD label RETURN label")['label'].tolist()
    nodes = {}
    for label in labels:
        query = f"MATCH (n:`{label}`) WITH n LIMIT 1 RETURN keys(n) AS keys"
        result = query_cypher(query)
        if not result.empty:
            key_columns = [col for col in result.columns if col.startswith('keys.')]
            if key_columns:
                keys = result[key_columns].values.flatten().tolist()
                keys = [k for k in keys if pd.notnull(k)]  
                nodes[label] = keys

    rtypes = query_cypher("CALL db.relationshipTypes() YIELD relationshipType AS rtype")['rtype'].tolist()
    relations = {}
    for rtype in rtypes:
        query = f"MATCH ()-[r:`{rtype}`]-() WITH r LIMIT 1 RETURN keys(r) AS keys"
        r_result = query_cypher(query)
        if not r_result.empty:
            key_columns = [col for col in r_result.columns if col.startswith('keys.')]
            if key_columns:
                keys = r_result[key_columns].values.flatten().tolist()
                keys = [k for k in keys if pd.notnull(k)]  
                relations[rtype] = keys
    
    return nodes, relations

def check_nodes_type():
    """
    Queries and retrieves all node types (labels) from the database.

    Returns:
        list: A list containing all node types present in the database.
    """
    return list(query_cypher("CALL db.labels() YIELD label RETURN label;")['label'])

def check_relation_type(source_entity_type=None, target_entity_type=None, directed=False):
    """
    Queries the database for relationship types between entities and their occurrence counts.

    Args:
        source_entity_type (str, optional): The type of the source entity. Defaults to None, meaning no restriction on type.
        target_entity_type (str, optional): The type of the target entity. Defaults to None, meaning no restriction on type.
        directed (bool, optional): Whether to consider the directionality of relationships. Defaults to False (undirected query).

    Returns:
        DataFrame: A DataFrame containing relationship types and their occurrence counts, sorted in descending order.
    """
    #all_nodes_type = check_nodes_type()
    #if node_type in all_nodes_type:
    source_entity_type = _stardard_type_str(source_entity_type)
    target_entity_type = _stardard_type_str(target_entity_type)
    if directed:
        directed_str = '>'
    else:
        directed_str = ''
    return query_cypher(f"MATCH (n{source_entity_type})-[r]-{directed_str}(m{target_entity_type}) \
                          RETURN DISTINCT type(r) AS relationship_type, COUNT(r) AS count \
                          ORDER BY count DESC;")

# TODO: 写一个函数可以查询每个relation可能的info_source，应该统计出来作为字典存在这里
def check_relation_info_source(relation=None, source_entity_type=None, target_node_type=None, ):
    """
    Queries the database for the information sources of a specific relationship.

    Args:
        relation (str, optional): The type of relationship. Defaults to None, meaning no restriction.
        source_entity_type (str, optional): The type of the source entity. Defaults to None, meaning no restriction.
        target_node_type (str, optional): The type of the target entity. Defaults to None, meaning no restriction.

    Returns:
        list: A list of database names that serve as information sources for the specified relationship.
    """
    relation_str = _stardard_type_str(relation)
    source_entity_type = _stardard_type_str(source_entity_type)
    target_node_type = _stardard_type_str(target_node_type)

    cypher = f"MATCH ({source_entity_type})-[r{relation_str}]-({target_node_type}) \
               UNWIND r.info_source AS i \
               return DISTINCT SPLIT(i, ':')[0] AS DB"
    df = query_cypher(cypher)
    return list(df['DB'])