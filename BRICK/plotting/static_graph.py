import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Ellipse
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

def compute_node_levels(G, root):
    """
    Compute the depth level of each node in a directed graph starting from a given root.

    This function performs a breadth-first traversal of the graph, assigning a depth 
    level to each node based on its distance from the root node.

    Args:
        G (networkx.DiGraph): The directed graph in which node levels are computed.
        root (hashable): The root node from which the levels are determined.

    Returns:
        dict: A dictionary mapping each node to its depth level, where the root 
              has a level of 0, and deeper nodes have increasing level values.
    """
    levels = {}
    queue = [(root, 0)]
    
    while queue:
        node, level = queue.pop(0)
        levels[node] = level
        for child in G.successors(node): 
            if child not in levels:
                queue.append((child, level + 1))
    
    return levels

def _auto_change_line(s, max_length=10):
    """
    Automatically inserts line breaks in a string if its length exceeds a specified limit.

    This function splits the input string by spaces and joins the words with newline characters
    if the total length of the string exceeds `max_length`. Otherwise, it returns the original string.

    Args:
        s (str): The input string to process.
        max_length (int, optional): The maximum allowed length before inserting line breaks. Defaults to 10.

    Returns:
        str: The processed string with line breaks if necessary.
    """
    if len(s) > max_length:
        return '\n'.join(s.split(' '))
    else:
        return s

def bezier_curve(p0, p1, p2, p3, num_points=100):
    """
    生成三阶贝塞尔曲线的点
    p0: 起点
    p1: 第一个控制点
    p2: 第二个控制点
    p3: 终点
    num_points: 生成曲线上的点数
    """
    t = np.linspace(0, 1, num_points)
    curve = (1-t)**3 * p0[:, None] + 3*(1-t)**2*t * p1[:, None] + 3*(1-t)*t**2 * p2[:, None] + t**3 * p3[:, None]
    return curve.T 

def compute_node_positions(G, root):
    """计算每个节点的位置，使得同一层级的节点在同一水平线上"""
    levels = compute_node_levels(G, root)
    max_level = max(levels.values())
    
    level_counts = {}
    for node, level in levels.items():
        if level not in level_counts:
            level_counts[level] = 0
        level_counts[level] += 1

    node_positions = {}
    for level in range(max_level + 1):
        nodes_in_level = [node for node, lvl in levels.items() if lvl == level]
        num_nodes = len(nodes_in_level)
        x_positions = np.linspace(0, num_nodes, num_nodes, endpoint=False) 
        x_positions -= np.mean(x_positions)  
        for i, node in enumerate(nodes_in_level):
            node_positions[node] = (x_positions[i], -level) 

    return node_positions, levels

def draw_tree(ax, node_positions, node_levels, final_G, node_colors, node_shapes, node_sizes, default_node_color, edge_colors, node_scale_factor):
    for node_name, (x, y) in node_positions.items():
        num_children = sum(1 for _ in final_G.successors(node_name))  
        text_y = y + 0.3 if num_children >= 3 else y - 0.15 

        color = node_colors.get(node_name, default_node_color)
        shape = node_shapes.get(node_name, 'circle')
        size = node_sizes.get(node_name, 10) / node_scale_factor

        if shape == 'circle':
            plt.scatter(x, y, s=200, color=color, zorder=2)
        elif shape == 'hexagon':
            hexagon = RegularPolygon((x, y), numVertices=6, radius=size, color=color, zorder=2) # radius=0.15写死了节点的大小
            ax.add_patch(hexagon)
        elif shape == 'ellipse':
            ellipse = Ellipse((x, y), width=0.3, height=0.15, color=color, zorder=2)
            ax.add_patch(ellipse)
        else:
            plt.scatter(x, y, s=200, color=color, zorder=2) 

        plt.text(x, text_y, _auto_change_line(node_name), fontsize=7, ha='center', va='top',zorder=4)

    for i, (parent, child) in enumerate(final_G.edges()):
        if parent in node_positions and child in node_positions:
            p0 = np.array(node_positions[parent])  
            p3 = np.array(node_positions[child]) 
            
            edge_data = final_G[parent][child]
            is_verified = any(inner_data.get('arrows', 'none') == 'to' for inner_data in edge_data.values()) 
            edge_color = edge_colors[i]  

            if node_levels[parent] == node_levels[child]: 
                mid_x = (p0[0] + p3[0]) / 2 
                mid_y = p0[1] + 0.5  
                p1 = np.array([mid_x, mid_y]) 
                p2 = np.array([mid_x, mid_y]) 

                curve = bezier_curve(p0, p1, p2, p3)
                plt.plot(curve[:, 0], curve[:, 1], color=edge_color, linewidth=2)
                
                if is_verified:
                    arrow = FancyArrowPatch(
                        (curve[-2, 0], curve[-2, 1]), 
                        (curve[-1, 0], curve[-1, 1]),  
                        arrowstyle='-|>', mutation_scale=20, color=edge_color, zorder=3
                    )
                    ax.add_patch(arrow)
            else:  
                if is_verified:
                    arrow = FancyArrowPatch(p0, p3, arrowstyle='wedge', mutation_scale=8, color=edge_color, zorder=3)
                    ax.add_patch(arrow)
                else:
                    plt.plot([p0[0], p3[0]], [p0[1], p3[1]], color=edge_color, linewidth=2)
    plt.axis("off")

def create_label(nxg, pos, default_node_color, node_shape, default_edge_color):
    node_legend_items = []
    drawn_nodes = set(pos.keys())
    for n in nxg.nodes:
        if n not in drawn_nodes:
            continue

        color = nxg.nodes[n].get('color', default_node_color)
        label = nxg.nodes[n].get('label', n)
        shape = nxg.nodes[n].get(node_shape, 'o') 
        
        marker = 'h' if shape == 'hexagon' else 'o'
        
        node_legend_items.append(
            plt.plot([], [], color=color, marker=marker, markersize=10, linestyle='None', label=label)[0]
        )
    edge_legend_items = []
    edge_styles = {}

    if isinstance(nxg, (nx.MultiGraph, nx.MultiDiGraph)):
        for u, v, k in nxg.edges(keys=True):
            edge_data = nxg[u][v][k]
            edge_color = edge_data.get('color', default_edge_color)
            edge_arrows = edge_data.get('arrows', 'none')
            key = (edge_color, edge_arrows)
            if key not in edge_styles:
                edge_styles[key] = f"{'Verified' if edge_arrows == 'to' else 'Unverified'} Edge"
    else:
        for u, v in nxg.edges:
            edge_data = nxg[u][v]
            edge_color = edge_data.get('color', default_edge_color)
            edge_arrows = edge_data.get('arrows', 'none')
            key = (edge_color, edge_arrows)
            if key not in edge_styles:
                edge_styles[key] = f"{'Verified' if edge_arrows == 'to' else 'Unverified'} Edge"

    for (color, arrows), label in edge_styles.items():
        line_style = '-'  
        edge_item = Line2D(
            [0], [1],
            color=color,
            lw=2, 
            linestyle=line_style,
            label=label
        )
        edge_legend_items.append(edge_item)

        plt.legend(
            handles=node_legend_items + edge_legend_items,
            loc='center left', bbox_to_anchor=(1, 0.5),
            title="Node and Edge Legend"
        )

def static_visualize_network(nxg=None, figsize=(8, 4), layout='spring', node_color='color', node_size='size', node_shape='shape',
    edge_color='color', edge_width='weight', edge_arrow='arrows', label=True, save_path=None, node_scale_factor=300,
    edge_scale_factor=5, default_node_color='gray', default_edge_color='black', frameon = True, **kwargs):

    if nxg is None:
        raise ValueError("No graph provided to visualize.")
    
    node_sizes = {
        n: nxg.nodes[n].get(node_size, 10)*20 if nxg.nodes[n].get(node_shape, 'circle') == 'circle' else nxg.nodes[n].get(node_size, 10)
        for n in nxg.nodes
    }
    node_shapes = {n: nxg.nodes[n].get(node_shape, 'circle') for n in nxg.nodes}
    node_colors = {n: nxg.nodes[n].get(node_color, default_node_color) for n in nxg.nodes}

    if isinstance(nxg, (nx.MultiGraph, nx.MultiDiGraph)):
        edge_colors = [nxg[u][v][k].get(edge_color, default_edge_color) for u, v, k in nxg.edges(keys=True)]
        edge_widths = [nxg[u][v][k].get(edge_width, 1) / edge_scale_factor for u, v, k in nxg.edges(keys=True)]
        edge_arrows = [nxg[u][v][k].get('arrows', 'none') == 'to' for u, v, k in nxg.edges(keys=True)]
        edge_list = list(nxg.edges(keys=True))
    else:
        edge_colors = [nxg.edges[e].get(edge_color, default_edge_color) for e in nxg.edges]
        edge_widths = [nxg.edges[e].get(edge_width, 1) / edge_scale_factor for e in nxg.edges]
        edge_arrows = [nxg.edges[e].get(edge_arrow, False) is True for e in nxg.edges]
        edge_arrows = [nxg.edges[e].get('arrows', 'none') == 'to' for e in nxg.edges]
        edge_list = list(nxg.edges)

    if layout == 'spring':
        pos = nx.spring_layout(nxg)
    elif layout == 'circular':
        pos = nx.circular_layout(nxg)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(nxg)
    elif layout == 'shell':
        pos = nx.shell_layout(nxg)
    elif layout == 'spectral':
        pos = nx.spectral_layout(nxg)
    elif layout == 'tree':
        try:
            root = [n for n, d in nxg.in_degree() if d == 0][0]
            pos, node_levels = compute_node_positions(nxg, root=root)
        except ValueError as e:
            print(f"Error: {e}. Attempting to convert the graph to a tree.")
            if isinstance(nxg, nx.DiGraph):
                root = next((n for n, d in nxg.in_degree() if d == 0), None)
                tree = nx.dfs_tree(nxg, source=root) if root else nxg
            else:
                tree = nx.minimum_spanning_tree(nxg)
            pos, node_levels = compute_node_positions(tree, root=root)

    plt.figure(figsize=figsize)
    ax = plt.gca()
    
    if layout != 'tree':
        for n, (x, y) in pos.items():
            shape = node_shapes[n]
            color = node_colors[n]
            size = node_sizes[n] / node_scale_factor

            if color == default_node_color:
                ax.text(
                    x, y, str(n),
                    fontsize=8,
                    ha='center',
                    va='center',
                    zorder=3,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.6)
                )
            
            if shape == 'circle':
                plt.scatter(x, y, s=node_sizes[n], c=color, zorder=2, edgecolor='k')
            elif shape == 'hexagon':
                hexagon = RegularPolygon((x, y), numVertices=6, radius=size, color=color, zorder=2)
                ax.add_patch(hexagon)
            elif shape == 'ellipse':
                ellipse = Ellipse((x, y), width=size * 1.5, height=size, color=color, zorder=2)
                ax.add_patch(ellipse)
            else:
                plt.scatter(x, y, s=size, c=color, zorder=2, edgecolor='k')
        
        for i, edge in enumerate(edge_list):
            u, v, *rest = edge
            if edge_arrows[i]: 
                nx.draw_networkx_edges(
                    nxg, pos,
                    edgelist=[(u, v)],
                    edge_color=edge_colors[i],
                    width=edge_widths[i],
                    arrowstyle='-|>',
                    arrowsize=10,
                )
            else:
                nx.draw_networkx_edges(
                    nxg, pos,
                    edgelist=[(u, v)],
                    edge_color=edge_colors[i],
                    width=edge_widths[i],
                    arrowstyle='-', 
                )
    else:
        draw_tree(ax, pos, node_levels, nxg, node_colors, node_shapes, node_sizes, default_node_color, edge_colors, node_scale_factor)
    
    if label:
        create_label(nxg, pos, default_node_color, node_shape, default_edge_color)
    
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, format=save_path.split(".")[-1])
        print(f"Graph saved to {save_path}")
    
    if not frameon:
        ax.set_axis_off()
    plt.show()