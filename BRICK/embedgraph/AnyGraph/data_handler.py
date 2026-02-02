import pickle
import numpy as np
from scipy.sparse import coo_matrix
import scipy.sparse as sp
from torch.utils import data
import torch as t

from .Utils.TimeLogger import log
from .model import Feat_Projector, Adj_Projector, TopoEncoder
#from . import params
#args = params.args


class MultiDataHandler:
    def __init__(self, datasets):
        all_datasets = datasets
        all_datasets.sort()
        self.handlers = []
        for data_name in all_datasets:
            handler = DataHandler(data_name)
            self.handlers.append(handler)

    def remake_initial_projections(self):
        for i in range(len(self.handlers)):
            trn_handler = self.handlers[i]
            trn_handler.make_projectors()


class DataHandler:
    def __init__(self, coo_mat, args, feat_path=None):
        #self.coo_path = coo_path
        self.args = args
        self.feat_path = feat_path
        self.topo_encoder = TopoEncoder()
        self.tst_input_adj = None
        self.asym_adj = None
        self.coo_mat = coo_mat
        self.coo_mat = self.load_data()
        self.trn_dataset = TrnData(self.coo_mat)

    def load_one_file(self):
        with open(self.coo_path, 'rb') as fs:
            ret = (pickle.load(fs) != 0).astype(np.float32)
        if type(ret) != coo_matrix:
            ret = sp.coo_matrix(ret)
        return ret
    
    def load_feats(self):
        if self.feat_path:
            try:
                with open(self.feat_path, 'rb') as fs:
                    feats = pickle.load(fs)
            except Exception as e:
                print(self.feat_path + str(e))
                exit()
        return feats

    def normalize_adj(self, mat, log=False):
        degree = np.array(mat.sum(axis=-1))
        d_inv_sqrt = np.reshape(np.power(degree, -0.5), [-1])
        d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.0
        d_inv_sqrt_mat = sp.diags(d_inv_sqrt)
        if mat.shape[0] == mat.shape[1]:
            return mat.dot(d_inv_sqrt_mat).transpose().dot(d_inv_sqrt_mat).tocoo()
        else:
            tem = d_inv_sqrt_mat.dot(mat)
            col_degree = np.array(mat.sum(axis=0))
            d_inv_sqrt = np.reshape(np.power(col_degree, -0.5), [-1])
            d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.0
            d_inv_sqrt_mat = sp.diags(d_inv_sqrt)
            return tem.dot(d_inv_sqrt_mat).tocoo()
    
    def unique_numpy(self, row, col):
        hash_vals = row * self.args.node_num + col
        hash_vals = np.unique(hash_vals).astype(np.int64)
        col = hash_vals % self.args.node_num
        row = (hash_vals - col).astype(np.int64) // self.args.node_num
        return row, col

    def make_torch_adj(self, mat, unidirectional_for_asym=False):
        if mat.shape[0] == mat.shape[1]:
            _row = mat.row
            _col = mat.col
            row = np.concatenate([_row, _col]).astype(np.int64)
            col = np.concatenate([_col, _row]).astype(np.int64)
            row, col = self.unique_numpy(row, col)
            data = np.ones_like(row)
            mat = coo_matrix((data, (row, col)), mat.shape)
            if self.args.selfloop == 1:
                mat = (mat + sp.eye(mat.shape[0])) * 1.0
            normed_asym_mat = self.normalize_adj(mat)
            row = t.from_numpy(normed_asym_mat.row).long()
            col = t.from_numpy(normed_asym_mat.col).long()
            idxs = t.stack([row, col], dim=0)
            vals = t.from_numpy(normed_asym_mat.data).float()
            shape = t.Size(normed_asym_mat.shape)
            asym_adj = t.sparse.FloatTensor(idxs, vals, shape)
            return asym_adj
        elif unidirectional_for_asym:
            mat = (mat != 0) * 1.0
            mat = self.normalize_adj(mat, log=True)
            idxs = t.from_numpy(np.vstack([mat.row, mat.col]).astype(np.int64))
            vals = t.from_numpy(mat.data.astype(np.float32))
            shape = t.Size(mat.shape)
            return t.sparse.FloatTensor(idxs, vals, shape)
        else:
            # make ui adj
            a = sp.csr_matrix((self.args.user_num, self.args.user_num))
            b = sp.csr_matrix((self.args.item_num, self.args.item_num))
            mat = sp.vstack([sp.hstack([a, mat]), sp.hstack([mat.transpose(), b])])
            mat = (mat != 0) * 1.0
            if self.args.selfloop == 1:
                mat = (mat + sp.eye(mat.shape[0])) * 1.0
            mat = self.normalize_adj(mat)

            # make cuda tensor
            idxs = t.from_numpy(np.vstack([mat.row, mat.col]).astype(np.int64))
            vals = t.from_numpy(mat.data.astype(np.float32))
            shape = t.Size(mat.shape)
            return t.sparse.FloatTensor(idxs, vals, shape)

    def load_data(self):
        #tst_mat = self.load_one_file()
        tst_mat = self.coo_mat

        if self.feat_path:
            self.feats = t.from_numpy(self.load_feats()).float()
            self.feats = self.feats
            self.args.featdim = self.feats.shape[1]
        else:
            self.feats = None
            self.args.featdim = self.args.latdim

        if tst_mat.shape[0] != tst_mat.shape[1]:
            self.args.user_num, self.args.item_num = tst_mat.shape
            self.args.node_num = self.args.user_num + self.args.item_num
            print('Dataset: User num: {user_num}, Item num: {item_num}, Node num: {node_num}'.format(user_num=self.args.user_num, item_num=self.args.item_num, node_num=self.args.node_num, edge_num=tst_mat.nnz))
        else:
            self.args.node_num = tst_mat.shape[0]
            print('Dataset: Node num: {node_num}, Edge num: {edge_num}'.format(node_num=self.args.node_num, edge_num=tst_mat.nnz))
        self.tst_input_adj = self.make_torch_adj(tst_mat)
        if tst_mat.shape[0] == tst_mat.shape[1]:
            self.asym_adj = self.tst_input_adj
        else:
            self.asym_adj = self.make_torch_adj(tst_mat, unidirectional_for_asym=True)
        self.make_projectors()
        return tst_mat

    def make_projectors(self):
        with t.no_grad():
            projectors = []
            if self.args.proj_method == 'adj_svd' or self.args.proj_method == 'both':
                tem = self.asym_adj.to(self.args.device)
                projectors = [Adj_Projector(tem)]
            if self.feats is not None and self.args.proj_method != 'adj_svd':
                tem = self.feats.to(self.args.device)
                projectors.append(Feat_Projector(tem))
            feats = projectors[0]()
            if len(projectors) == 2:
                feats2 = projectors[1]()
                feats = feats + feats2

            try:
                self.projectors = self.topo_encoder(self.tst_input_adj.to(self.args.device), feats.to(self.args.device)).detach().cpu()
            except Exception:
                mean, std = feats.mean(dim=-1, keepdim=True), feats.std(dim=-1, keepdim=True)
                tem_adj = self.tst_input_adj.to(self.args.device)
                mem_cache = 256
                projectors_list = []
                for i in range(feats.shape[1] // mem_cache):
                    st, ed = i * mem_cache, (i + 1) * mem_cache
                    tem_feats = (feats[:, st:ed] - mean) / (std + 1e-8)
                    tem_feats = self.topo_encoder(tem_adj, tem_feats.to(self.args.device), normed=True).detach().cpu()
                    projectors_list.append(tem_feats)
                self.projectors = t.concat(projectors_list, dim=-1)
            t.cuda.empty_cache()


class TstData(data.Dataset):
    def __init__(self, coomat, trn_mat):
        self.csrmat = (trn_mat.tocsr() != 0) * 1.0
        tstLocs = [None] * coomat.shape[0]
        tst_nodes = set()
        for i in range(len(coomat.data)):
            row = coomat.row[i]
            col = coomat.col[i]
            if tstLocs[row] is None:
                tstLocs[row] = list()
            tstLocs[row].append(col)
            tst_nodes.add(row)
        tst_nodes = np.array(list(tst_nodes))
        self.tst_nodes = tst_nodes
        self.tstLocs = tstLocs

    def __len__(self):
        return len(self.tst_nodes)

    def __getitem__(self, idx):
        return self.tst_nodes[idx]


class TrnData(data.Dataset):
    def __init__(self, coomat):
        self.ancs, self.poss = coomat.row, coomat.col
        self.negs = np.zeros(len(self.ancs)).astype(np.int32)
        self.cand_num = coomat.shape[1]
        self.neg_shift = 0 if coomat.shape[0] == coomat.shape[1] else coomat.shape[0]
        self.poss = coomat.col + self.neg_shift
        self.neg_sampling()

    def neg_sampling(self):
        self.negs = np.random.randint(self.cand_num + self.neg_shift, size=self.poss.shape[0])

    def __len__(self):
        return len(self.ancs)

    def __getitem__(self, idx):
        return self.ancs[idx], self.poss[idx], self.negs[idx]


class JointTrnData(data.Dataset):
    def __init__(self, dataset_list):
        self.batch_dataset_ids = []
        self.batch_st_ed_list = []
        self.dataset_list = dataset_list
        for dataset_id, dataset in enumerate(dataset_list):
            samp_num = len(dataset) // self.args.batch + (1 if len(dataset) % self.args.batch != 0 else 0)
            for j in range(samp_num):
                self.batch_dataset_ids.append(dataset_id)
                st = j * self.args.batch
                ed = min((j + 1) * self.args.batch, len(dataset))
                self.batch_st_ed_list.append((st, ed))

    def neg_sampling(self):
        for dataset in self.dataset_list:
            dataset.neg_sampling()

    def __len__(self):
        return len(self.batch_dataset_ids)

    def __getitem__(self, idx):
        st, ed = self.batch_st_ed_list[idx]
        dataset_id = self.batch_dataset_ids[idx]
        return *self.dataset_list[dataset_id][st: ed], dataset_id