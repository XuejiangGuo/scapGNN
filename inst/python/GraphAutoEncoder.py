import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
import scipy.sparse as sp
import numpy as np


def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):
        sparse_mx = sparse_mx.tocoo()
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape


def preprocess_graph(adj):
    adj = sp.coo_matrix(adj)
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return sparse_to_tuple(adj_normalized)


def mask_test_edges(adj, ratio_val):
    adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
    adj.eliminate_zeros()

    assert np.diag(adj.todense()).sum() == 0

    adj_triu = sp.triu(adj)
    adj_tuple = sparse_to_tuple(adj_triu)
    edges = adj_tuple[0]
    num_val = int(np.floor(edges.shape[0] * ratio_val))

    all_edge_idx = list(range(edges.shape[0]))
    np.random.shuffle(all_edge_idx)
    val_edge_idx = all_edge_idx[:num_val]

    val_edges = edges[val_edge_idx]
    train_edges = np.delete(edges, val_edge_idx, axis=0)

    data = np.ones(train_edges.shape[0])

    adj_train = sp.csr_matrix((data, (train_edges[:, 0], train_edges[:, 1])), shape=adj.shape)
    adj_train = adj_train + adj_train.T

    return adj_train, train_edges, val_edges


class VGAE(nn.Module):
    def __init__(self, adj, input_dim, hidden1_dim, hidden2_dim):
        super(VGAE, self).__init__()
        self.hidden2_dim = hidden2_dim
        self.base_gcn = GraphConvSparse(input_dim, hidden1_dim, adj)
        self.gcn_mean = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x: x)
        self.gcn_logstddev = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x: x)

    def encode(self, X):
        hidden = self.base_gcn(X)
        self.mean = self.gcn_mean(hidden)
        self.logstd = self.gcn_logstddev(hidden)
        gaussian_noise = torch.randn(X.size(0), self.hidden2_dim)
        sampled_z = gaussian_noise * torch.exp(self.logstd) + self.mean
        return sampled_z

    def forward(self, X):
        Z = self.encode(X)
        A_pred = dot_product_decode(Z)
        return A_pred, Z


class GraphConvSparse(nn.Module):
    def __init__(self, input_dim, output_dim, adj, activation=F.relu, **kwargs):
        super(GraphConvSparse, self).__init__(**kwargs)
        self.weight = glorot_init(input_dim, output_dim)
        self.adj = adj
        self.activation = activation

    def forward(self, inputs):
        x = inputs
        x = torch.mm(x, self.weight)
        x = torch.mm(self.adj, x)
        outputs = self.activation(x)
        return outputs


def dot_product_decode(Z):
    A_pred = torch.sigmoid(torch.matmul(Z, Z.t()))
    return A_pred


def glorot_init(input_dim, output_dim):
    init_range = np.sqrt(6.0 / (input_dim + output_dim))
    initial = torch.rand(input_dim, output_dim) * 2 * init_range - init_range
    return nn.Parameter(initial)


class GAE(nn.Module):
    def __init__(self, adj, input_dim, hidden1_dim, hidden2_dim):
        super(GAE, self).__init__()
        self.base_gcn = GraphConvSparse(input_dim, hidden1_dim, adj)
        self.gcn_mean = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x: x)

    def encode(self, X):
        hidden = self.base_gcn(X)
        z = self.mean = self.gcn_mean(hidden)
        return z

    def forward(self, X):
        Z = self.encode(X)
        A_pred = dot_product_decode(Z)
        return A_pred, Z


def GAE_function(net_m, feature_m, use_model=True, GAE_epochs=300, GAE_learning_rate=0.01, seed=1217, ratio_val=0.05, verbose=True):
    GAE_epochs = int(GAE_epochs)

    net_m = sp.coo_matrix(net_m)

    idx = np.argwhere(np.all(feature_m[..., :] == 0, axis=0))
    feature_m = np.delete(feature_m, idx, axis=1)

    for i in range(feature_m.shape[1]):
        feature_m[:, i] = feature_m[:, i] / feature_m[:, i].max()
    input_dim = feature_m.shape[1]
    hidden1_dim = int(input_dim / 2)
    hidden2_dim = int(hidden1_dim / 2)
    feature_m = sp.coo_matrix(feature_m)

    adj_orig = net_m
    adj_orig = adj_orig - sp.dia_matrix((adj_orig.diagonal()[np.newaxis, :], [0]), shape=adj_orig.shape)
    adj_orig.eliminate_zeros()

    adj_train, train_edges, val_edges = mask_test_edges(net_m, ratio_val)
    adj = net_m
    adj_norm = preprocess_graph(adj)

    features = sparse_to_tuple(feature_m.tocoo())
    # num_features = features[2][1]
    # features_nonzero = features[1].shape[0]

    pos_weight = float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()
    norm = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)

    adj_label = adj_train + sp.eye(adj_train.shape[0])
    adj_label = sparse_to_tuple(adj_label)

    adj_norm = torch.sparse.FloatTensor(torch.LongTensor(adj_norm[0].T),
                                        torch.FloatTensor(adj_norm[1]),
                                        torch.Size(adj_norm[2]))
    adj_label = torch.sparse.FloatTensor(torch.LongTensor(adj_label[0].T),
                                         torch.FloatTensor(adj_label[1]),
                                         torch.Size(adj_label[2]))
    features = torch.sparse.FloatTensor(torch.LongTensor(features[0].T),
                                        torch.FloatTensor(features[1]),
                                        torch.Size(features[2]))

    weight_mask = adj_label.to_dense().view(-1) == 1
    weight_tensor = torch.ones(weight_mask.size(0))
    weight_tensor[weight_mask] = pos_weight

    torch.manual_seed(int(seed))
    if use_model:
        model = VGAE(adj=adj_norm, input_dim=input_dim, hidden1_dim=hidden1_dim, hidden2_dim=hidden2_dim)
        optimizer = Adam(model.parameters(), lr=GAE_learning_rate)
    else:
        model = GAE(adj=adj_norm, input_dim=input_dim, hidden1_dim=hidden1_dim, hidden2_dim=hidden2_dim)
        optimizer = Adam(model.parameters(), lr=GAE_learning_rate)

    min_loss_val = 10000
    # broc = 0
    for epoch in range(GAE_epochs):

        A_pred, z = model(features)

        optimizer.zero_grad()

        try:
            loss = norm * F.binary_cross_entropy(A_pred.view(-1), adj_label.to_dense().view(-1), weight=weight_tensor)
        except:
            continue

        if use_model:
            kl_divergence = 0.5 / A_pred.size(0) * (
                        1 + 2 * model.logstd - model.mean ** 2 - torch.exp(model.logstd) ** 2).sum(1).mean()
            loss -= kl_divergence

        loss.backward()
        optimizer.step()

        if verbose:
             print("Progress (%): {}".format(round(epoch / GAE_epochs, 2) * 100), end="\r")

        if loss <= min_loss_val and epoch > 5:
            min_loss_val = loss.item()
            # broc = val_roc
            best_adj = A_pred
            best_epoch = epoch

    # test_roc, test_ap = get_scores(test_edges, test_edges_false, best_adj, adj_orig)
    if verbose:
        print("Minimum loss: {}".format(min_loss_val), end="\r")
        print("Minimum loss - epoch : {}".format(best_epoch), end="\r")

    best_adj = np.array(best_adj.detach())
    #z = np.array(z.detach())
    return best_adj, min_loss_val, best_epoch
