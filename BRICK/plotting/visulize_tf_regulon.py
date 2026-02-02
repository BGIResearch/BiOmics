import networkx as nx
import matplotlib.pyplot as plt

def visualize_tf_regulon(TF, top_tf, verified_gene, save_path=None):
    """
    Visualizes the interaction network between transcription factors (TFs) and their regulated genes (Regulons).

    This function generates a network graph illustrating the relationship between transcription factors 
    and their target genes, highlighting verified regulatory interactions.

    Args:
        TF (list): A list of target transcription factors.
        top_tf (pd.DataFrame): A DataFrame containing the regulated genes for each transcription factor. 
            Must include a column named "Regulon".
        verified_gene (list): A list of experimentally verified regulated genes.
        save_path (str, optional): If provided, saves the figure to the specified path. Defaults to None.

    Returns:
        None: Displays the network graph or saves it to the specified path.
    """
    G = nx.Graph()

    regulon = top_tf.loc[TF, "Regulon"].explode().tolist()

    intersection = set(verified_gene) & set(regulon)

    for gene in regulon:
        G.add_edge(TF[0], gene, color="gray")

    edge_colors = [data["color"] for _, _, data in G.edges(data=True)]

    plt.figure(figsize=(8, 6))
    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, nodelist=TF, node_color="#ffb3c2", node_shape="s", node_size=800)
    nx.draw_networkx_nodes(G, pos, nodelist=regulon, node_color="#b8b8b8", node_shape="o", node_size=800)
    nx.draw_networkx_nodes(G, pos, nodelist=intersection, node_color="#b3ffb4", node_shape="o", node_size=800)

    nx.draw_networkx_edges(G, pos, edge_color=edge_colors)
    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.title("TF-Regulon Interaction Network")
    plt.gca().set_axis_off()

    if save_path is not None:
        plt.savefig(save_path, format=save_path.split(".")[-1])
        print(f"Graph saved to {save_path}")

    plt.show()


