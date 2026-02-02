from tqdm import tqdm
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from torch.optim import Adam
from torch_geometric.nn import GCNConv
from torch_geometric.utils import from_networkx, to_networkx
import GCL.losses as L
import GCL.augmentors as A
from GCL.eval import get_split, LREvaluator
from GCL.models import DualBranchContrast


class GConv(torch.nn.Module):
    """
    Graph Convolutional Network (GCN) module.
    This class implements a multi-layer Graph Convolutional Network (GCN) using
    the GCNConv layer from PyTorch Geometric.

    Args:
        input_dim (int): The dimensionality of the input features.
        hidden_dim (int): The dimensionality of the hidden layers.
        activation (callable): The activation function to be applied after each 
            convolutional layer.
        num_layers (int): The number of GCN layers in the network.

    Attributes:
        activation (torch.nn.Module): The activation function instance.
        layers (torch.nn.ModuleList): A list of GCNConv layers.

    Methods:
        forward(x, edge_index, edge_weight=None):
            Performs the forward pass of the GCN.

    """
    def __init__(self, input_dim, hidden_dim, activation, num_layers):
        """
        Initializes the GConv module.

        Args:
            input_dim (int): The dimensionality of the input features.
            hidden_dim (int): The dimensionality of the hidden layers.
            activation (callable): The activation function to be applied.
            num_layers (int): The number of GCN layers.
        """
        super(GConv, self).__init__()
        self.activation = activation()
        self.layers = torch.nn.ModuleList()
        self.layers.append(GCNConv(input_dim, hidden_dim, cached=False))
        for _ in range(num_layers - 1):
            self.layers.append(GCNConv(hidden_dim, hidden_dim, cached=False))

    def forward(self, x, edge_index, edge_weight=None):
        """
        Performs a forward pass through the GCN.

        Args:
            x (torch.Tensor): The node feature matrix of shape (num_nodes, input_dim).
            edge_index (torch.Tensor): The edge indices of shape (2, num_edges),
                defining the graph structure.
            edge_weight (torch.Tensor, optional): Edge weights of shape (num_edges,).
                Defaults to None.

        Returns:
            torch.Tensor: The output feature matrix after applying GCN layers.
        """
        z = x
        for i, conv in enumerate(self.layers):
            z = conv(z, edge_index, edge_weight)
            z = self.activation(z)
        return z


class Encoder(torch.nn.Module):
    """
    Graph Encoder with Data Augmentation and Projection.

    This module applies graph data augmentation before encoding node features 
    using a given encoder model. It also includes a projection head for further 
    transformation of the encoded representations.

    Args:
        encoder (torch.nn.Module): A neural network module for encoding graph 
            node features.
        augmentor (tuple): A pair of augmentation functions applied to the input 
            graph before encoding.
        hidden_dim (int): The dimensionality of the hidden representation.
        proj_dim (int): The dimensionality of the projected representation.

    Attributes:
        encoder (torch.nn.Module): The encoder model for graph node features.
        augmentor (tuple): A pair of augmentation functions.
        fc1 (torch.nn.Linear): The first fully connected layer in the projection head.
        fc2 (torch.nn.Linear): The second fully connected layer in the projection head.

    Methods:
        forward(x, edge_index, edge_weight=None):
            Applies augmentations and encodes the input graph.
        project(z):
            Transforms the encoded representation using the projection head.
    """
    def __init__(self, encoder, augmentor, hidden_dim, proj_dim):
        """
        Initializes the Encoder module.

        Args:
            encoder (torch.nn.Module): A model for encoding graph node features.
            augmentor (tuple): A pair of data augmentation functions.
            hidden_dim (int): Hidden dimension of the encoded representations.
            proj_dim (int): Output dimension of the projection head.
        """
        super(Encoder, self).__init__()
        self.encoder = encoder
        self.augmentor = augmentor

        self.fc1 = torch.nn.Linear(hidden_dim, proj_dim)
        self.fc2 = torch.nn.Linear(proj_dim, hidden_dim)

    def forward(self, x, edge_index, edge_weight=None):
        """
        Performs forward propagation through the encoder.

        This function applies two augmentations to the input graph and encodes 
        both augmented versions along with the original graph.

        Args:
            x (torch.Tensor): The node feature matrix of shape (num_nodes, feature_dim).
            edge_index (torch.Tensor): The edge indices of shape (2, num_edges),
                defining the graph structure.
            edge_weight (torch.Tensor, optional): Edge weights of shape (num_edges,).
                Defaults to None.

        Returns:
            tuple: A tuple containing:
                - z (torch.Tensor): The encoded representation of the original graph.
                - z1 (torch.Tensor): The encoded representation of the first augmented graph.
                - z2 (torch.Tensor): The encoded representation of the second augmented graph.
        """
        aug1, aug2 = self.augmentor
        x1, edge_index1, edge_weight1 = aug1(x, edge_index, edge_weight)
        x2, edge_index2, edge_weight2 = aug2(x, edge_index, edge_weight)
        z = self.encoder(x, edge_index, edge_weight)
        z1 = self.encoder(x1, edge_index1, edge_weight1)
        z2 = self.encoder(x2, edge_index2, edge_weight2)
        return z, z1, z2

    def project(self, z: torch.Tensor) -> torch.Tensor:
        """
        Applies a projection head to the encoded representation.

        The projection head consists of an ELU activation followed by a linear transformation.

        Args:
            z (torch.Tensor): The input feature matrix.

        Returns:
            torch.Tensor: The projected feature representation.
        """
        z = F.elu(self.fc1(z))
        return self.fc2(z)


class GRACE_model():
    """
    Graph Contrastive Learning Model using GRACE.

    This class implements a self-supervised graph contrastive learning model based on 
    the GRACE framework. It applies graph augmentations, encodes node representations, 
    and optimizes contrastive loss.

    Args:
        nxg (networkx.Graph): The input graph in NetworkX format.
        hidden_dim (int, optional): The dimensionality of the hidden layer. Default is 32.
        num_layers (int, optional): The number of GCN layers. Default is 2.
        lr (float, optional): The learning rate for optimization. Default is 0.0001.
        device (str, optional): The computation device ('cpu' or 'cuda'). Default is 'cpu'.
        seed (int, optional): The random seed for reproducibility. Default is 20241204.

    Attributes:
        nxg (networkx.Graph): The original input graph.
        data (torch_geometric.data.Data): The graph data converted for PyTorch Geometric.
        gconv (GConv): The graph convolutional encoder.
        encoder_model (Encoder): The encoder that applies augmentations and encodes the graph.
        contrast_model (DualBranchContrast): The contrastive learning model.
        optimizer (torch.optim.Adam): The optimizer for training.

    Methods:
        train():
            Trains the model for one iteration using contrastive loss.
        test():
            Evaluates the model using a logistic regression classifier.
        get_embedding():
            Returns the learned node embeddings.
        fit(epoches=1000):
            Trains the model for a specified number of epochs.
        plot_loss():
            Plots the training loss curve.
        get_result(umap=True):
            Returns the graph embeddings as an AnnData object, with optional UMAP visualization.
    """
    def __init__(self, nxg, hidden_dim=32, num_layers=2, lr=0.0001, device='cpu', seed=20241204):
        """
        Initializes the GRACE_model.

        Args:
            nxg (networkx.Graph): The input graph.
            hidden_dim (int, optional): Hidden dimension size. Default is 32.
            num_layers (int, optional): Number of GCN layers. Default is 2.
            lr (float, optional): Learning rate. Default is 0.0001.
            device (str, optional): Computation device ('cpu' or 'cuda'). Default is 'cpu'.
            seed (int, optional): Random seed. Default is 20241204.
        """
        torch.manual_seed(seed)
        self.aug1 = A.Compose([A.EdgeRemoving(pe=0.3), A.FeatureMasking(pf=0.3)])
        self.aug2 = A.Compose([A.EdgeRemoving(pe=0.3), A.FeatureMasking(pf=0.3)])

        self.nxg = nxg
        data = from_networkx(nxg)
        data.x = torch.from_numpy(np.diag(np.ones(data.num_nodes))).type(torch.float32)
        self.data = data.to(device)
        self.gconv = GConv(input_dim=self.data.num_nodes, hidden_dim=hidden_dim, activation=torch.nn.ReLU,
                           num_layers=num_layers).to(device)
        self.encoder_model = Encoder(encoder=self.gconv, augmentor=(self.aug1, self.aug2), hidden_dim=hidden_dim,
                                     proj_dim=hidden_dim).to(device)
        self.contrast_model = DualBranchContrast(loss=L.InfoNCE(tau=0.2), mode='L2L', intraview_negs=True).to(device)
        self.optimizer = Adam(self.encoder_model.parameters(), lr=lr)

    def train(self):
        """
        Performs one training step using contrastive learning.

        This function computes node embeddings, applies a projection head, and 
        optimizes contrastive loss.

        Returns:
            float: The computed loss value for the training step.
        """
        self.encoder_model.train()
        self.optimizer.zero_grad()
        z, z1, z2 = self.encoder_model(self.data.x, self.data.edge_index, self.data.edge_attr)
        h1, h2 = [self.encoder_model.project(x) for x in [z1, z2]]
        loss = self.contrast_model(h1, h2)
        loss.backward()
        self.optimizer.step()
        return loss.item()

    def test(self):
        """
        Evaluates the learned embeddings using logistic regression.

        Returns:
            dict: Evaluation metrics, including accuracy.
        """
        self.encoder_model.eval()
        z, _, _ = self.encoder_model(self.data.x, self.data.edge_index, self.data.edge_attr)
        split = get_split(num_samples=z.size()[0], train_ratio=0.1, test_ratio=0.8)
        result = LREvaluator()(z, self.data.y, split)
        return result

    def get_embbedding(self):
        """
        Retrieves the learned node embeddings.

        Returns:
            numpy.ndarray: Node embeddings in NumPy format.
        """
        # encoder_model.eval()
        z, _, _ = self.encoder_model(self.data.x, self.data.edge_index, self.data.edge_attr)
        return z.to('cpu').detach().numpy()

    def fit(self, epoches=1000):
        """
        Trains the model for a specified number of epochs.

        Args:
            epoches (int, optional): The number of training epochs. Default is 1000.
        """
        with tqdm(total=epoches, desc='(T)') as pbar:
            self.losses = []
            for epoch in range(1, epoches + 1):
                loss = self.train()
                self.losses.append(loss)
                pbar.set_postfix({'loss': loss})
                pbar.update()

    def plot_loss(self):
        """
        Plots the training loss over epochs.
        """
        import matplotlib.pyplot as plt
        plt.plot(self.losses)

    def get_result(self, umap=True):
        """
        Generates the node embeddings and returns them as an AnnData object.

        This function extracts node embeddings and constructs an AnnData object 
        for downstream analysis. Optionally, it applies UMAP for visualization.

        Args:
            umap (bool, optional): Whether to compute UMAP embedding. Default is True.

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


