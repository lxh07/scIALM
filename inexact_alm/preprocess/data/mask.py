# %%
import os
import copy
import pandas as pd
import scanpy as sc
import numpy as np

from scipy import sparse


def maskPerCol(data_train, masked_prob):
    """
    将表达矩阵中每列非零的值随机置为0并返回，同时返回置为0的元素的坐标
    :param data_train: 表达矩阵
    :param masked_prob: 置0比例
    :return:
    """
    X_train = copy.deepcopy(data_train)
    rows = []
    cols = []
    for col in range(data_train.shape[1]):
        index_pair_train = np.where(data_train[:, col])
        if index_pair_train[0].shape[0] <= 3:
            continue
        masking_idx_train = np.random.choice(index_pair_train[0].shape[0],
                                             int(index_pair_train[0].shape[0] * masked_prob),
                                             replace=False)
        X_train[index_pair_train[0][masking_idx_train], [col] * masking_idx_train.shape[0]] = 0
        for i in index_pair_train[0][masking_idx_train]:
            rows.append(i)
            cols.append(col)

    return X_train, rows, cols
# %%

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--masked_prob', default=0.1, type=float)
parser.add_argument('--dataset', default='PBMC', type=str)

args = parser.parse_args()

adata = sc.read_h5ad('%s/processed/%s.h5ad' % (args.dataset, args.dataset))
sc.pp.normalize_total(adata)
adata.raw = adata

# %%

path = '%s/masked' % args.dataset

if not os.path.exists(path):
    os.makedirs(path)

pd.DataFrame(adata.raw.X.A).to_csv(path + '/un_%s.csv' % args.dataset)
masked, masking_row, masking_col = maskPerCol(adata.raw.X.A, args.masked_prob)

print(masked.shape)
pd.DataFrame(masked, index=adata.obs.index, columns=adata.var.index) \
    .to_csv(path + '/%s_%s.csv' % (args.dataset, str(args.masked_prob).replace('.', '')))

adata.X = sparse.csr_matrix(masked)
adata.write(path + '/%s_%s.h5ad' % (args.dataset, str(args.masked_prob).replace('.', '')))

# %%

maskIndex = sparse.coo_matrix(([1] * len(masking_col), (masking_row, masking_col)))
newindex = pd.concat([pd.Series(masking_row), pd.Series(masking_col)], axis=1)
print(newindex)
pd.DataFrame(newindex).to_csv(path + '/%s_maskIndex_%s.csv' % (args.dataset, str(args.masked_prob).replace('.', '')))
sparse.save_npz(path + '/%s_maskIndex_%s.csv' % (args.dataset, str(args.masked_prob).replace('.', '')), maskIndex)
