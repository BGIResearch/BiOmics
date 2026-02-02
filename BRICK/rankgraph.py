from collections import defaultdict

import pandas as pd
from scipy.stats import fisher_exact, hypergeom

from .querygraph import query_cypher, query_relation

def _query_df2count_df(df, group_col = 'path.2.name', 
                       retained_col = ['path.0.name', 
                                       'path.1.relation', 'path.1.info_source_length', 'path.1.relation_confidence', 
                                       'path.2.id', 'path.2.name', 'path.2.type']):
    """
    Groups the given DataFrame by a specified column and processes the retained columns, returning a new DataFrame.

    This function takes a DataFrame, groups it based on a specified column, and processes the retained columns 
    as per the requirements. The result is a new DataFrame that contains the processed information.

    Args:
        df (pandas.DataFrame): The input DataFrame to be processed.
        group_col (str, default='path.2.name'): The column name to group by.
        retained_col (list of str, optional): A list of column names to retain and process. Default includes 'path.0.name', 
                                              'path.1.relation', 'path.1.info_source_length', 
                                              'path.1.relation_confidence', 'path.2.id', 'path.2.name', and 'path.2.type'.

    Returns:
        pandas.DataFrame: A new DataFrame with the processed results after grouping by the specified column 
                          and retaining the specified columns.
    """
    group_path = '.'.join(group_col.split('.')[:2])

    retained_col = [x for x in retained_col if x in df]
    df_count = defaultdict(list)
    for x, y in df.groupby(group_col):
        #df_count[group_col].append(x)
        for i in retained_col:
            tmp = list(y[i])
            if group_path in i:
                df_count[i].append(tmp[0])
            else:
                df_count[i].append(tmp)
    
    df_count = pd.DataFrame(df_count)
    return df_count

def _prune_parent(query_df):
    """
    Prune the parent nodes of the query DataFrame.

    Args:
        query_df (pandas.DataFrame): The input DataFrame to be processed.
        
    Returns:
        pandas.DataFrame: A new DataFrame with the processed results after pruning the parent nodes.
    """
    all_possible_nodes = list(query_df['path.2.name'].unique())
    all_possible_types = "|".join(query_df['path.2.type'].unique())
    parent2children = query_relation(all_possible_nodes, target_entity_type=all_possible_types, relation='is_a', directed=True)
    parent_cells = [x for x, y in parent2children.groupby('path.2.name') if y.shape[0] == 1]
    only_child_parent = parent2children.loc[parent2children['path.2.name'].isin(parent_cells)]
    only_child_parent = dict(zip(only_child_parent['path.2.name'], only_child_parent['path.0.name']))
    query_df['path.2.name'] = [ only_child_parent[x] if x in only_child_parent else x for x in query_df['path.2.name'] ]
    return query_df

def match_count(df, target_group_df=None,  prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame 
    with statistical results.

    This function allows grouping the input DataFrame by the specified columns, processes the retained columns, 
    and computes statistical results based on those columns. It can also integrate additional statistics from 
    a pre-grouped target DataFrame when multiple sorting indicators are used.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values. 
                                                      This is used when computing additional statistics 
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can be 
                                   passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed 
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and calculated results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    df_count['path.2.match_count'] = [len(set(x)) for x in df_count['path.0.name']]
    return df_count.sort_values(by='path.2.match_count', ascending=False)

def info_source_count(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame 
    with statistical results, including the count of information sources.

    This function allows grouping the input DataFrame by specified columns, processes the retained columns, 
    and computes the total count of information sources based on the 'path.1.info_source_length' column.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values. 
                                                      This is used when computing additional statistics 
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can be 
                                   passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed 
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and calculated results, 
                           including the information source count.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    df_count['path.2.info_source_count'] = [sum(x) for x in df_count['path.1.info_source_length']]
    return df_count.sort_values(by='path.2.info_source_count', ascending=False)

def match_probability(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame 
    with statistical results, including the calculated match probability.

    This function groups the input DataFrame based on specified columns, processes the retained columns, and 
    calculates the match probability for each group, providing a new DataFrame with these results.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values. 
                                                      This is used when calculating additional statistics 
                                                      with multiple sorting indicators.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can 
                                   be passed via **kwargs.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed 
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and match probability results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df

    adj = defaultdict(dict)
    for (a, b), y in df.groupby(['path.0.name', 'path.2.name']):
        adj[a][b] = 1

    adj = pd.DataFrame(adj)
    adj.fillna(0, inplace=True)
    adj = adj / adj.sum()
    target2prob = dict(adj.mean(axis=1))
    df_count['path.2.match_probability'] = [target2prob[x] for x in df_count['path.2.name']]
    return df_count.sort_values(by='path.2.match_probability', ascending=False)

def info_source_probability(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame 
    with statistical results, including the calculated information source probability.

    This function groups the input DataFrame based on specified columns, processes the retained columns, and 
    calculates the information source probability for each group, providing a new DataFrame with these results.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values. 
                                                      This is used when calculating additional statistics 
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can 
                                   be passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed 
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and information source probability results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df

    adj = defaultdict(dict)
    for (a, b), y in df.groupby(['path.0.name', 'path.2.name']):
        adj[a][b] = y['path.1.info_source_length'].sum()

    adj = pd.DataFrame(adj)
    adj.fillna(0, inplace=True)
    adj = adj / adj.sum()
    target2prob = dict(adj.mean(axis=1))
    df_count['path.2.info_source_probability'] = [target2prob[x] for x in df_count['path.2.name']]
    return df_count.sort_values(by='path.2.info_source_probability', ascending=False)


def enrich(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, 
    and returns a new DataFrame with statistical results.

    This function enhances the input DataFrame by grouping it based on specified columns 
    and processing the retained columns. If a pre-grouped DataFrame is provided, 
    new computed indicators are added to it.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values. 
                                                      Used when additional sorting indicators need 
                                                      to be added to the existing results.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        group_col (str, optional): The column name to group by. The default is 'path.2.name', 
                                   and it can be passed via **kwargs.
        retained_col (list of str, optional): A list of column names to retain and process. 
                                              This can be passed via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistical information 
                          and computed results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    source_type = ':' + '|'.join(list(df['path.0.type'].unique()))
    relation_type = ':' + '|'.join(list(df['path.1'].unique()))
    target_type = ':' + '|'.join(list(df['path.2.type'].unique()))
    
    target_set = list(df_count['path.2.id'])
    source_set_length = df['path.0.id'].nunique()    
    targetid2bg_count = query_cypher(f"MATCH (n{target_type})-[r{relation_type}]-(m{source_type}) WHERE n.id IN $target_set return n.id, count(m)", 
                                      parameters={"target_set":target_set})
    
    targetid2bg_count = dict(zip(targetid2bg_count['n.id'], targetid2bg_count['count(m)']))

    df_count = match_count(df, target_group_df=df_count, **kwargs)
    #df_count['path.2.match_count'] = [len(x) for x in df_count['path.0.name']]
    df_count['path.2.background_count'] = [targetid2bg_count[x] for x in df_count['path.2.id']]

    test_list = [fisher_exact([[x, y], [source_set_length, 20000]], alternative='greater') \
                          for x, y in zip(df_count['path.2.match_count'], df_count['path.2.background_count'])]

    df_count[['path.2.enrich_statistic', 'path.2.enrich_pvalue']] = [[x.statistic, x.pvalue] for x in test_list]
    return df_count.sort_values(by='path.2.enrich_pvalue', ascending=True)



def rank_voting(df, target_group_df=None,
                metrics=['match_count', 'match_probability', 'info_source_count', 'info_source_probability'],
                ascending=False,
                prune_parent=False,
                **kwargs):
    """
    Ranks data using multiple metrics and applies a voting mechanism to determine the final ranking score.

    This function sorts the given DataFrame based on multiple specified metrics, assigns ranking scores, 
    and aggregates them to produce a final ranking score.

    Args:
        df (pandas.DataFrame): The input DataFrame containing data to be ranked.
        target_group_df (pandas.DataFrame, optional): A pre-processed grouped DataFrame. 
                                                      If None (default), the ranking metrics are computed.
        metrics (list of str, optional): A list of metric names used for ranking. 
                                         The default includes 'match_count', 'match_probability', 
                                         'info_source_count', and 'info_source_probability'.
        ascending (bool or list of bool, optional): Sorting order. If a single boolean is provided, 
                                                    all metrics follow the same sorting order. 
                                                    If a list is given, each metric follows its respective order.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        kwargs: Additional parameters passed to the `_query_df2count_df` function.

    Returns:
        pandas.DataFrame: A DataFrame with ranking scores, including a new column 'path.2.rank_voting' 
                          representing the final ranking score.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    for x in metrics:
        df_count = eval(x)(df, df_count)

    if isinstance(ascending, bool):
        ascending_list = [ascending] * len(metrics)
    elif isinstance(ascending, list):
        ascending_list = [bool(x) for x in ascending]
    else:
        raise ValueError('ascending must be bool type or a list of bool type')
    rank_list = []
    for met, asc in zip(metrics, ascending_list):
        if met == 'enrich':
            met = 'enrich_pvalue'
            asc = True
        rank_list.append(df_count['path.2.' + met].rank(ascending=asc))

    rank_list = pd.DataFrame(rank_list).T
    df_count['path.2.rank_voting'] = list(rank_list.mean(axis=1))
    return df_count.sort_values(by='path.2.rank_voting', ascending=True)

