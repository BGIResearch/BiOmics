from pyvis.network import Network

def interact_visualize_network(nxg, node_label = 'name',
                               width = 500, height = 500, notebook=True, directed=True, save='nx.html'):
    """
    Visualizes a given NetworkX graph as an interactive network diagram.

    Args:
        nxg (networkx.Graph): The NetworkX graph object to be visualized.
        width (int, optional): The width of the visualization interface in pixels. Defaults to 500.
        height (int, optional): The height of the visualization interface in pixels. Defaults to 500.
        notebook (bool, optional): Whether to display the visualization in a Jupyter Notebook environment. Defaults to True.
        directed (bool, optional): Whether the graph is directed. Defaults to True.
        save (str, optional): The filename to save the visualization output. Defaults to 'nx.html'.

    Returns:
        object: The generated visualization object.
    """
    for n in nxg.nodes():
        nxg.nodes[n]['label'] = nxg.nodes[n][node_label]

    nt = Network(f'{height}px', f'{width}px', notebook=notebook, directed=directed)
    nt.from_nx(nxg)
    if notebook:
        nt.show(save)
    else:
        nt.write_html(save) 
    return nt
