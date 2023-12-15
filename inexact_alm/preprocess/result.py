# %%
import numpy as np
import pandas as pd
import scanpy as sc

from scipy import sparse
from scipy.stats import pearsonr

from sklearn.cluster import KMeans
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score

import argparse

# %%

parser = argparse.ArgumentParser()
parser.add_argument('--masked_prob', default=0.1, type=float)
parser.add_argument('--dataset', default='Klein', type=str)

args = parser.parse_args()

adata = sc.read_h5ad(
    'data/%s/masked/%s_%s.h5ad' % (
        args.dataset, args.dataset, str(args.masked_prob).replace('.', '')))

maskIndex = sparse.load_npz(
    'data/%s/masked/%s_maskIndex_%s.csv.npz' % (
        args.dataset, args.dataset, str(args.masked_prob).replace('.', '')))

(masking_row, masking_col) = np.where(maskIndex.A > 0)

# %%

df = pd.read_csv(
    'data/%s/result/%s_%s.csv' % (
    args.dataset, args.dataset, str(args.masked_prob).replace('.', '')))
df[df<0] = 0
# %%
def pearsonr_error(y, h):
    res = []
    if len(y.shape) < 2:
        y = y.reshape((1, -1))
        h = h.reshape((1, -1))

    for i in range(y.shape[0]):
        res.append(pearsonr(y[i], h[i])[0])
    return np.mean(res)


def cosine_similarity_score(y, h):
    if len(y.shape) < 2:
        y = y.reshape((1, -1))
        h = h.reshape((1, -1))
    cos = cosine_similarity(y, h)
    res = []
    for i in range(len(cos)):
        res.append(cos[i][i])
    return np.mean(res)


# %%

y = adata.raw.X.A[masking_row, masking_col]
h = df.values[masking_row, masking_col]

mse = float('%.4f' % mean_squared_error(y, h))
mae = float('%.4f' % mean_absolute_error(y, h))
pear = float('%.4f' % pearsonr_error(y, h))
cos = float('%.4f' % cosine_similarity_score(y, h))

clus = {}
if 'cluster' in adata.obs:
    clusters = adata.obs.cluster.values

    adata_pred = sc.AnnData(df.values)

    sc.pp.normalize_total(adata_pred, target_sum=1e4)
    sc.pp.log1p(adata_pred)

    kmeans = KMeans(n_clusters=len(set(clusters))).fit(adata_pred.X)
    ari = float('%.4f' % adjusted_rand_score(clusters, kmeans.labels_))
    nmi = float('%.4f' % normalized_mutual_info_score(clusters, kmeans.labels_))

    clus = {
        'ari': ari,
        'nmi': nmi
    }
    print(
        '%s,%s,%s,%s,%s,%s' % (mse, mae, pear, cos, ari, nmi))

else:
    print('%s,%s,%s,%s' % (mse, mae, pear, cos))
