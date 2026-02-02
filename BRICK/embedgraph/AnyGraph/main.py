from collections import defaultdict

import torch
from torch import nn
import networkx as nx
import numpy as np
import scanpy as sc
import pandas as pd

from .data_handler import DataHandler
from .model import AnyGraph
from .params import args

class AnyGraph_model(nn.Module):
    """
    A model for learning node embeddings using AnyGraph.

    This class integrates a graph neural network (AnyGraph) with additional 
    transformations to process graph data and extract meaningful node embeddings.

    Args:
        nxg (networkx.Graph): The input graph in NetworkX format.
        model_path (str): Path to the pre-trained AnyGraph model.
        kg_dim (int, optional): Input feature dimension from the knowledge graph. Default is 1024.
        out_dim (int, optional): Output feature dimension. Default is 512.
        freeze (bool, optional): Whether to freeze the AnyGraph model's parameters. Default is True.
        n_layers (int, optional): Number of additional transformation layers. Default is 3.
        **kwargs: Additional arguments for configuring the model.

    Attributes:
        nxg (networkx.Graph): The input graph.
        model (AnyGraph): The AnyGraph model instance.
        data_handler (DataHandler): Handles graph data conversion and processing.
        fc (torch.nn.Linear): Fully connected layer for feature transformation.
        layers (torch.nn.ModuleList): Additional transformation layers.
    """
    def __init__(self, nxg, model_path, kg_dim=1024, out_dim=512,  freeze=True, n_layers=3, **kwargs):
        """
        Initializes the AnyGraph-based model.

        Args:
            nxg (networkx.Graph): The input graph.
            model_path (str): Path to the pre-trained AnyGraph model.
            kg_dim (int, optional): Knowledge graph feature dimension. Default is 1024.
            out_dim (int, optional): Output embedding dimension. Default is 512.
            freeze (bool, optional): Whether to freeze the AnyGraph model parameters. Default is True.
            n_layers (int, optional): Number of additional transformation layers. Default is 3.
            **kwargs: Additional arguments for model configuration.
        """
        super(AnyGraph_model, self).__init__()
        self.nxg = nxg
        for key, value in kwargs.items():
            if hasattr(args, key):
                setattr(args, key, value)
            else:
                print(f"Warning: {args} has no attribute named {key}")
                
        self.args = args
        self.model = self.init_anygraph(model_path)
        
        coo_mat = nx.to_scipy_sparse_array(nxg).tocoo()
        data_handler = DataHandler(coo_mat, self.args)
        self.data_handler = data_handler
        self.fc = nn.Linear(kg_dim, out_dim)
        self.layers = nn.ModuleList()
        for i in range(n_layers):
            self.layers.append(nn.Linear(out_dim, out_dim))
            self.layers.append(nn.ReLU())
            self.layers.append(nn.LayerNorm(out_dim))
        if freeze:
            for param in self.model.parameters():
                param.requires_grad = False

    #@staticmethod
    def init_anygraph(self, model_path):
        """
        Initializes the AnyGraph model and loads pre-trained weights if available.

        Args:
            model_path (str): Path to the pre-trained AnyGraph model.

        Returns:
            AnyGraph: The initialized AnyGraph model.
        """
        model = AnyGraph(self.args)
        if model_path:
            state_dict = torch.load(model_path, map_location='cpu')
            model.load_state_dict(state_dict)
        return model

    def forward(self):
        """
        Forward pass through the AnyGraph model to compute node embeddings.

        Returns:
            torch.Tensor: Node embeddings of shape (num_nodes, out_dim).
        """
        self.model.assign_experts([self.data_handler], reca=True, log_assignment=True)
        expert = self.model.summon(0)
        feats = self.data_handler.projectors
        try:
            final_embeds = expert.forward(feats)
        except Exception:
            final_embeds_list = []
            div = 256 * 3
            temlen = feats.shape[0] // div
            for i in range(temlen):
                st, ed = div * i, div * (i + 1)
                tem_projectors = feats[st: ed, :]
                final_embeds_list.append(expert.forward(tem_projectors))
            if temlen * div < feats.shape[0]:
                tem_projectors = feats[temlen * div:, :]
                final_embeds_list.append(expert.forward(tem_projectors))
            final_embeds = torch.concat(final_embeds_list, dim=0)
        return final_embeds

    def get_embbedding(self):
        """
        Computes and retrieves node embeddings.

        Returns:
            numpy.ndarray: The computed node embeddings.
        """
        # encoder_model.eval()
        z = self.forward()
        return z.to('cpu').detach().numpy()

    def get_result(self, umap=True):
        """
        Computes embeddings and returns them as an AnnData object.

        This function extracts node embeddings, constructs an AnnData object, 
        and optionally applies UMAP for visualization.

        Args:
            umap (bool, optional): Whether to compute UMAP visualization. Default is True.

        Returns:
            anndata.AnnData: An AnnData object containing node embeddings and metadata.
        """
        emb = self.get_embbedding()
        # embeddings = {x: y for x, y in zip(nxg.nodes, emb)}
        expr = []
        obs = defaultdict(list)

        for x, y in zip(self.nxg.nodes, emb):
            expr.append(y)
            obs['name'].append(x)
            obs['type'].append(self.nxg.nodes[x]['type'])

        expr_numpy = np.array(expr)
        adata_graph = sc.AnnData(expr_numpy, obs=pd.DataFrame(obs))
        adata_graph.obs_names = adata_graph.obs['name']
        # sc.pp.pca(adata_graph)
        if umap:
            sc.pp.neighbors(adata_graph)
            sc.tl.umap(adata_graph)
        return adata_graph