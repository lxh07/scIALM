# %%
import os
import scanpy as sc
from scipy import sparse

# %%

adataD0 = sc.read_csv('GSM1599494_ES_d0_main.csv.bz2')
adataD2 = sc.read_csv('GSM1599497_ES_d2_LIFminus.csv.bz2')
adataD4 = sc.read_csv('GSM1599498_ES_d4_LIFminus.csv.bz2')
adataD7 = sc.read_csv('GSM1599499_ES_d7_LIFminus.csv.bz2')

# %%

adata = sc.AnnData.concatenate(adataD0.T, adataD2.T, adataD4.T, adataD7.T, batch_key='cluster',
                               batch_categories=['d0', 'd2', 'd4', 'd7', ])
adata.X = sparse.csr_matrix(adata.X)

# %%

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.total_counts < 75000, :]

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.raw = adata

# %%
folder = os.path.exists('processed')

if not folder:
    os.makedirs('processed')

adata.write('processed/Klein.h5ad')
#text = pd.DataFrame(adata.raw.X.A, index=adata.obs.index, columns=adata.var.index)
#text.to_csv('processed/Klein.csv')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=40)
sc.tl.leiden(adata)

marker_genes = ['Rpl7', 'Rpl3', 'Tubb5', 'Pou5f1', 'Sox2', 'Zfp42', 'Dppa5a', 'Nanog', 'Ccnb1', 'Dnmt3b', 'Dppa3',
                'Zfp281', 'Ly6a', 'Krt18', 'Krt8', 'Col4a1', 'Col4a2', 'Plec', 'Peg10', 'Actb',
                'Actg1', 'Slc2a1', 'Lgmn', 'Ngfrap1', 'Fabp3']

sc.pl.heatmap(adata, marker_genes, groupby='leiden', dendrogram=True, swap_axes=True, use_raw=True, cmap='Blues')