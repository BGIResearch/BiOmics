import scanpy as sc
from scipy import sparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


def radial_enrich_plot(query_target_df, adata_graph=None, n_neighbors=3, top_k=7, leiden_resolution=1, palette='tab10',
              target_entity_type=None, return_cluster_result=False):
    """
    Analyzes gene enrichment within target biological entities (e.g., cells, tissues, diseases) 
    and visualizes the results using a radial plot.

    This function constructs a gene network and visualizes enrichment relationships 
    based on gene-target entity associations provided in `query_target_df`:
    - Computes gene co-occurrence relationships and constructs an `AnnData` object.
    - Performs Leiden clustering and UMAP dimensionality reduction.
    - Generates a radial plot where target entities are distributed along the circumference, 
      and genes are positioned at the center.

    Args:
        query_target_df (pd.DataFrame): A DataFrame containing gene-target entity associations. 
            Must include the `path.0.name` column, which stores lists of gene names in pathways.
        adata_graph (anndata.AnnData, optional): A precomputed `AnnData` object containing gene embeddings.
            If None (default), the graph is constructed based on co-occurrence relationships.
        n_neighbors (int, optional): Number of neighbors used to compute gene-gene connectivity, 
            affecting graph sparsity. Defaults to 3.
        top_k (int, optional): Number of top target entities selected for enrichment analysis, 
            sorted by `pvalue`. Defaults to 7.
        leiden_resolution (float, optional): Resolution parameter for Leiden clustering, controlling 
            clustering granularity. Defaults to 1.
        palette (str, optional): Color palette used for visualization. Defaults to 'tab10'.
        target_entity_type (str or list, optional): Specifies the types of target entities to include. 
            Available options include `'Cell'`, `'Tissue'`, `'Disease'`, `'Phenotype'`, `'Chemical'`, 
            `'Process'`, `'Function'`, `'Cell_Component'`, and `'Pathway'`. Defaults to all types.
        return_cluster_result (bool, optional): Whether to return Leiden clustering results. Defaults to False.

    Returns:
        matplotlib.figure.Figure: Radial enrichment analysis plot.
        anndata.AnnData (optional): Computed `AnnData` object containing UMAP and Leiden clustering results, 
            returned only if `return_cluster_result=True`.
    """
    gene_nodes_df = set()
    for x in query_target_df['path.0.name']:
        for y in x:
            gene_nodes_df.add(y)

    gene_nodes_df = pd.DataFrame({'gene': list(gene_nodes_df)})
    gene_nodes_df.set_index('gene', inplace=True)

    if target_entity_type is None:
        target_entity_type = ['Cell', 'Tissue', 'Disease', 'Phenotype', 'Chemical', 'Process', 'Function',
                              'Cell_Component', 'Pathway']
    if isinstance(target_entity_type, str):
        target_entity_type = [target_entity_type]

    if adata_graph is None:
        adata = sc.AnnData(obs=gene_nodes_df)

        tmp = np.zeros((adata.shape[0], adata.shape[0]))
        tmp = pd.DataFrame(tmp)
        tmp.columns = adata.obs_names
        tmp.index = adata.obs_names

        for x in query_target_df['path.0.name']:
            for i in range(len(x) - 1):
                for j in range(i + 1, len(x)):
                    a = x[i]
                    b = x[j]
                    tmp[a][b] += 1
                    tmp[b][a] += 1

        for x in tmp:
            top_gene = tmp.sort_values(x, ascending=False).head(n_neighbors).index
            tmp.loc[[x not in top_gene for x in tmp.index], x] = 0

        # tmp[tmp < min_edge_weight] = 0
        tmp = sparse.csr_array(tmp)

        # tmp = tmp/tmp.max()

        adata.obsp['connectivities'] = tmp
        adata.uns['neighbors'] = {'connectivities_key': 'connectivities',
                                  'distances_key': 'distances',
                                  'params': {'n_neighbors': n_neighbors,
                                             'method': 'umap',
                                             'random_state': 0,
                                             'metric': 'euclidean'}}
    else:
        gl = list(gene_nodes_df.index)
        gl_emb = adata_graph[gl].X
        adata = sc.AnnData(X=gl_emb, obs=gene_nodes_df)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, )

    sc.tl.umap(adata)
    sc.tl.leiden(adata, leiden_resolution)

    # color map
    def _rgb_to_hex(rgb):
        """
        Converts an RGB color to a hexadecimal color code.

        Args:
            rgb (tuple[int, int, int]): A tuple containing three integer values (R, G, B), 
                each ranging from 0 to 255.

        Returns:
            str: The corresponding hexadecimal color code in the format `#RRGGBB`.
        """
        r, g, b = rgb
        hex_code = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex_code

    def _hex_to_rgb(hex_code):
        """
        Converts a hexadecimal color code to a normalized RGB color value.

        Args:
            hex_code (str): The hexadecimal color code in the format `#RRGGBB` or `RRGGBB`.

        Returns:
            tuple[float, float, float]: The normalized RGB color values as floating-point numbers (R, G, B), 
                each ranging from 0 to 1.
        """
        hex_code = hex_code.lstrip('#')  # 移除十六进制代码中的 "#"
        r = int(hex_code[0:2], 16)  # 提取红色分量并转换为十进制
        g = int(hex_code[2:4], 16)  # 提取绿色分量并转换为十进制
        b = int(hex_code[4:6], 16)  # 提取蓝色分量并转换为十进制
        return (r / 255, g / 255, b / 255)

    colors = matplotlib.colormaps.get_cmap(palette)
    colors = [_rgb_to_hex([int(y * 255) for y in x]) for x in colors.colors]

    leiden2colors = {x: y for x, y in zip(adata.obs['leiden'].unique(), colors)}
    target_type2colors = {x: y for x, y in zip(query_target_df['path.2.type'].unique(), colors[::-1])}

    def _add_Arc(ax, radius, theta, length=2, color='black'):
        """ 
        Draws an arc-like line segment on the specified matplotlib axis.
        
        This function converts polar coordinates to Cartesian coordinates 
        and draws a line segment with the specified radius, angle, length, and color.
        It is commonly used to add decorative or identifying arc-like elements in visualizations.

        Args:
            ax (matplotlib.axes.Axes): The axis object on which to draw the arc-like segment.
            radius (float): The starting radius of the arc, determining its initial distance from the origin.
            theta (float): The angle (in radians) at which the arc segment is positioned.
            length (float, optional): The length of the arc segment. The actual extension length 
                is calculated as `0.01 + length * 0.3`. Defaults to 2.
            color (str, optional): The color of the arc segment. Defaults to `'black'`.

        Returns:
            None
        """
        ax.plot([radius * np.cos(theta), (radius + 0.01 + length * 0.3) * np.cos(theta)],
                [radius * np.sin(theta), (radius + 0.01 + length * 0.3) * np.sin(theta)],
                c=color,
                linewidth=5, solid_capstyle='butt'
                )

    def adjust_theta4count_df(query_target_df, min_interval_degree=3):
        """ 
        Adjusts the `theta` angles in `query_target_df` to ensure a minimum spacing between them.
        
        This function ensures that the `theta` values, after sorting, have a minimum angular separation, 
        preventing overly dense distributions. It iteratively adjusts `theta` values to maintain 
        at least `min_interval_degree` degrees (converted to radians) while preserving the overall 
        distribution trend.

        Args:
            query_target_df (pandas.DataFrame): A DataFrame containing a `theta` column representing 
                the angular positions of target entities.
            min_interval_degree (float, optional): The minimum angle interval (in degrees) between 
                consecutive `theta` values. Defaults to 3 degrees. The value is converted to radians 
                for computation.

        Returns:
            pandas.DataFrame: The adjusted `query_target_df` with an additional `theta_adj` column 
            representing the modified `theta` values.
        """
        tmp_df = query_target_df.sort_values('theta').copy()
        min_interval_rad = np.radians(min_interval_degree)

        theta_list = list(tmp_df['theta'])
        theta_list = [-np.pi] + theta_list
        delta_theta = [theta_list[i + 1] - theta_list[i] for i in range(len(theta_list) - 1)]
        sum_delta_theta = sum(delta_theta)
        sum_delta_theta = min(6.23, sum_delta_theta)
        delta_theta = np.array(delta_theta)

        count = 0
        while delta_theta.min() < min_interval_rad:
            delta_theta = np.sqrt(0.01 + delta_theta ** 2)
            delta_theta_sum_new = np.sum(delta_theta)
            delta_theta = delta_theta * sum_delta_theta / delta_theta_sum_new
            count += 1
            # print(delta_theta.min())
            if count == 5:
                break

        new_theta = []
        for i in range(len(delta_theta)):
            new_theta.append(-np.pi + delta_theta[:i + 1].sum())

        tmp_df['theta_adj'] = new_theta
        tmp_df['theta_delta'] = tmp_df['theta'] - tmp_df['theta_adj']
        tmp_df['theta_adj'] = tmp_df['theta_adj'] + tmp_df['theta_delta'].mean()
        target2theta_adj = {x: y for x, y in zip(tmp_df['path.2.name'], tmp_df['theta_adj'])}
        query_target_df['theta_adj'] = [target2theta_adj[x] for x in query_target_df['path.2.name']]
        return query_target_df

    fig, ax = plt.subplots(1, figsize=(6, 6))

    # postion of genes
    pos = adata.obsm['X_umap'].T
    pos_center = pos.mean(axis=1, keepdims=True)
    pos_ptp = np.abs(pos - pos_center).max(axis=1, keepdims=True)
    pos = (pos - pos_center) / pos_ptp
    gene2pos = {x: y for x, y in zip(adata.obs_names, pos.T)}

    ax.scatter(pos[0], pos[1], c=[leiden2colors[x] for x in adata.obs['leiden']], zorder=10)

    df_count_concat = []
    for i, x in query_target_df.groupby('path.2.type'):
        if i in target_entity_type:
            df_count_concat.append(x.sort_values('path.2.enrich_pvalue', ascending=True).head(top_k))
    query_target_df = pd.concat(df_count_concat)

    query_target_df['source_length_scaled'] = query_target_df['path.2.match_count'] / query_target_df[
        'path.2.background_count']
    # get theta
    theta_list = []
    for _, row in query_target_df.iterrows():
        tmp_center = np.array([gene2pos[x] for x in row['path.0.name']]).mean(axis=0)
        theta = np.arctan2(tmp_center[1], tmp_center[0])
        theta_list.append(theta)

    query_target_df['theta'] = theta_list

    query_target_df = adjust_theta4count_df(query_target_df)
    # return query_target_df
    base_r = 1.5
    query_target_df['x'] = np.cos(query_target_df['theta_adj']) * base_r
    query_target_df['y'] = np.sin(query_target_df['theta_adj']) * base_r

    #
    query_target_df['c'] = [target_type2colors[x] for x in query_target_df['path.2.type']]

    for _, row in query_target_df.iterrows():
        tmp_x = row['x']
        tmp_y = row['y']
        tmp_s = row['path.2.name']
        tmp_alpha = 1 - (0.8 / query_target_df['path.2.enrich_pvalue'].max()) * row['path.2.enrich_pvalue']
        tmp_theta = row['theta_adj']
        tmp_length = row['source_length_scaled'] / query_target_df['source_length_scaled'].max()
        tmp_color = row['c']

        # add target entity
        _add_Arc(ax, base_r, tmp_theta, length=tmp_length, color=tmp_color)

        for g_pos in [gene2pos[x] for x in row['path.0.name']]:
            g_x, g_y = g_pos
            ax.plot([g_x, tmp_x], [g_y, tmp_y], alpha=tmp_alpha, linewidth=.1, color='red', zorder=1)

        d = np.degrees(tmp_theta)
        if (d > 90) or (d < -90):
            d = d + 180
            ax.text(tmp_x * 1.2, tmp_y * 1.2, tmp_s, va='center', ha='right', rotation=d, rotation_mode='anchor')
        else:
            ax.text(tmp_x * 1.2, tmp_y * 1.2, tmp_s, va='center', ha='left', rotation=d, rotation_mode='anchor')

    # 添加圆到 Axes
    # circle = patches.Circle((0,0), 1.41, edgecolor='black', fill=False)
    # ax.add_patch(circle)

    for t, c in target_type2colors.items():
        ax.plot([0, 0], [0, 0], linewidth=5, color=c, label=t)

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)
    # ax.set_xlim([-1.6, 1.6])
    ax.legend(frameon=False)
    ax.set_aspect('equal')
    if return_cluster_result:
        return fig, adata
    else:
        return fig