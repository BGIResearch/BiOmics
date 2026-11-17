import ast
import networkx as nx
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend to support environments without a GUI
import matplotlib.pyplot as plt
import pandas as pd

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
    # Validate inputs
    if not isinstance(TF, list):
        raise TypeError(f"TF must be a list, got {type(TF)}")
    
    if len(TF) == 0:
        raise ValueError("TF list cannot be empty")
    
    if not isinstance(top_tf, pd.DataFrame):
        raise TypeError(f"top_tf must be a pandas DataFrame, got {type(top_tf)}")
    
    if "Regulon" not in top_tf.columns:
        raise ValueError(f"top_tf must contain 'Regulon' column. Available columns: {list(top_tf.columns)}")
    
    # Check that every TF is present in the DataFrame
    missing_tfs = [tf for tf in TF if tf not in top_tf.index]
    if missing_tfs:
        raise ValueError(f"TFs not found in top_tf index: {missing_tfs}. "
                        f"Available TFs: {list(top_tf.index)}")
    
    G = nx.Graph()

    # Collect regulon genes for all TFs
    all_regulon = set()
    for tf in TF:
        tf_regulon = top_tf.loc[tf, "Regulon"]
        # Handle lists, sets, and string representations of lists
        if isinstance(tf_regulon, (list, set)):
            genes = set(tf_regulon)
        elif isinstance(tf_regulon, str):
            # Try to parse string representations such as "['Gene1', 'Gene2']"
            try:
                parsed = ast.literal_eval(tf_regulon)
                if isinstance(parsed, (list, set, tuple)):
                    genes = set(parsed)
                else:
                    genes = {str(parsed)}
            except (ValueError, SyntaxError):
                # If parsing fails, treat it as a single string value
                genes = {tf_regulon}
        else:
            genes = {str(tf_regulon)}
        all_regulon.update(genes)
        # Add an edge from each TF to each of its regulon genes
        for gene in genes:
            G.add_edge(tf, gene, color="gray")

    regulon = list(all_regulon)
    intersection = set(verified_gene) & all_regulon

    edge_colors = [data["color"] for _, _, data in G.edges(data=True)]

    plt.figure(figsize=(8, 6))
    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, nodelist=TF, node_color="#ffb3c2", node_shape="s", node_size=800)
    nx.draw_networkx_nodes(G, pos, nodelist=regulon, node_color="#b8b8b8", node_shape="o", node_size=800)
    nx.draw_networkx_nodes(G, pos, nodelist=list(intersection), node_color="#b3ffb4", node_shape="o", node_size=800)

    nx.draw_networkx_edges(G, pos, edge_color=edge_colors)
    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.title("TF-Regulon Interaction Network")
    plt.gca().set_axis_off()

    if save_path is not None:
        plt.savefig(save_path, format=save_path.split(".")[-1])
        print(f"Graph saved to {save_path}")

