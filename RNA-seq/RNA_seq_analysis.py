import numpy as np
import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import MinMaxScaler
import scanpy as sc
import anndata
from itertools import cycle, islice

color = np.array(list(islice(cycle(['#e6194b','#ffd8b1','#ffe119','#0082c8','#f58231','#911eb4','#46f0f0','#f032e6','#d2f53c','#fabebe','#808000','#e6beff','#aa6e28','#800000','#aaffc3','#008080','#3cb44b','#000080','#d1b26f','#fffac8','#650021','#808080','#000000', '#F4D03F', '#F39C12', '#FFAB91', '#FF00FF', '#CB4335','#27AE60','#76D7C4','#3498DB','#8E44AD','#6D4C41','#1A237E','#FFFF66','#004D40','#00FF33','#CCFFFF','#A1887F','#FFCDD2','#999966','#212121']), 100)))

#Fig. 1g
datasel = pd.read_csv('Fig1g.csv', index_col=0)
fig,ax = plt.subplots(figsize=(14,8))
ax.bar(np.arange(datasel.shape[0]), datasel, color=np.concatenate([['blue']*np.sum(datasel['Log2(FC)']<-0.5),['red']*np.sum(datasel['Log2(FC)']>0.5)]))
ax.set_xticks(np.arange(datasel.shape[0]))
ax.set_xticklabels(datasel.index, fontsize=6, rotation=90)
ax.set_xlim([-1,datasel.shape[0]])
ax.set_ylabel('log2(FoldChange)', fontsize=16)
ax.text(0, 0.2,'Pluripotent genes', color='blue', fontsize=12)
ax.text(105, -0.4,'Totipotent genes', color='red', fontsize=12)
ax.set_title('8_2.5d vs. mES \n Translational changes of pluripotent and totipotent genes', fontsize=16)
plt.savefig(filepath, transparent=True)


#Fig. 1h
#For up part
simdata = pd.read_csv('Fig1h_up.csv', index_col=0)
#For down part
simdata = pd.read_csv('Fig1h_down.csv', index_col=0)

corrx = np.concatenate(np.concatenate([[np.arange(simdata.shape[1])]*simdata.shape[0]]))
corry = np.concatenate([[x]*simdata.shape[1] for x in range(simdata.shape[0])])
scolor = np.concatenate(np.array(simdata[::-1]))

fig,ax = plt.subplots(figsize=(10,3))
plot = ax.scatter(corrx, corry, c = scolor, cmap='coolwarm',alpha=0.8,vmin=-1,vmax=1)
ax.set_xticks(np.arange(simdata.shape[1]))
ax.set_yticks(np.arange(simdata.shape[0]))
ax.set_xticklabels(datasc.index, rotation=60,rotation_mode='anchor', ha='right')
ax.set_yticklabels(datasc.columns[::-1])
ax.set_xlim([-0.5,simdata.shape[1]-0.5])
ax.set_ylim([-0.5,simdata.shape[0]-0.5])
cbar = plt.colorbar(plot,ax=ax,shrink=0.3)
vmin, vmax = -1,1
cbar.solids.set_clim([vmin, vmax])
cbar.set_ticks([vmin, vmax])
cbar.set_label('Normalized Gene expression')
cbar.draw_all()
plt.tight_layout()
plt.savefig(filepath, transparent=True)
plt.close()


#Fig. 1j
dataall = pd.read_csv('Fig1j.csv', index_col=0)
Z = linkage(dataall.T, method='complete', metric='minkowski', optimal_ordering=True)
fig,ax = plt.subplots(figsize=(8,4))
dn = dendrogram(Z, labels=dataall.columns, leaf_rotation=90, distance_sort='ascending')
plt.tight_layout()
plt.savefig(filepath, transparent=True)
plt.close()


#Fig. 2a
dataallsel = pd.read_csv('Fig2a.csv', index_col=0)
dataallsel = np.log2((dataallsel+1).astype(float))
samplesel = ['zygote_sc', '2C_sc','early2C_sc','mid2C_sc', 'late2C_sc', '4C_sc', '8C_sc', '16C_sc', 'ES_sc', 'midblast_sc','earlyblast_sc', 'lateblast_sc', 'D.EPSCs_sc', 'L.EPSCs_sc', 'TBLC_sc','ciTotiSC_P1', 'ciTotiSC_P2', 'ciTotiSC_P4', 'ciTotiSC_P8', '2CLC_bulk']

plotmarker3=['.','o','v','^','<','>','4','8','s','p','|','h','+','H','1','*','*','*','*', '2', '3', 'H', 'x', 'X', 'D', 'd', '_']

pca = PCA(n_components=5)
data_reduce = pca.fit_transform(dataallsel.T)

fig, ax = plt.subplots(figsize = (12,6))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .8, 1], elev=5, azim=200)
plt.cla()
ax.set_xlabel('PC-1', fontsize=15)
ax.set_ylabel('PC-2', fontsize=15)
ax.set_zlabel('PC-3', fontsize=15)
for i in range(len(samplesel)):
    ax.scatter(data_reduce[i][0],data_reduce[i][1],data_reduce[i][2], s=40, c=color5[i], label=samplesel[i], marker=plotmarker3[i], alpha=0.8)

ax.legend(markerscale=2, ncol = 2, prop={'size': 10}, bbox_to_anchor=(1,0.8), loc='upper left', fontsize=18)
plt.tight_layout()
plt.savefig(filepath, transparent=True)
plt.close()


#Fig. 2b
dataall = pd.read_csv('Fig2b.csv', index_col=0)
corr = np.array([[pearsonr(dataall.iloc[:,i], dataall.iloc[:,j])[0] for i in range(len(dataall.columns))] for j in range(len(dataall.columns))])
simdata = corr[[1,2,3,4,5,6,7,8,9]][:,[0,10,11,12,13,14]]

scaler = MinMaxScaler()
scaler.fit(zscore(simdata[::-1], axis=0))
scolor = np.concatenate(simdata[::-1])
ssize = np.concatenate(scaler.transform(zscore(simdata[::-1], axis=0)))
corrx = np.concatenate(np.concatenate([[np.arange(simdata.shape[1])]*simdata.shape[0]]))
corry = np.concatenate([[x]*simdata.shape[1] for x in range(simdata.shape[0])])

fig,ax = plt.subplots(figsize=(3.5,4))
plot = ax.scatter(corrx, corry, c = scolor, s = ssize*200, cmap='YlOrRd',vmin=-0.03,vmax=0.35,alpha=0.8)
ax.set_xticks(np.arange(simdata.shape[1]))
ax.set_yticks(np.arange(simdata.shape[0]))
ax.set_xticklabels(dataall.columns[[0,10,11,12,13,14]], rotation=60,rotation_mode='anchor', ha='right')
ax.set_yticklabels(dataall.columns[[1,2,3,4,5,6,7,8,9]][::-1])
ax.set_xlim([-0.5,simdata.shape[1]-0.5])
ax.set_ylim([-0.5,simdata.shape[0]-0.5])
cbar = plt.colorbar(plot,ax=ax,shrink=0.3)
vmin, vmax = -0.03, 0.35
cbar.solids.set_clim([vmin, vmax])
cbar.set_ticks([vmin, vmax])
cbar.set_label('Pearson Correlation')
cbar.draw_all()
plt.tight_layout()
plt.savefig(filepath, transparent=True)
plt.close()


#Fig. 2c
data = pd.read_csv('Fig2c.csv', index_col=0)
celltypeall1 = np.load('Fig2c_meta.npy')
cellmeta = pd.DataFrame(celltypeall1, index=data.columns)
expr = anndata.AnnData(X=data.T, var=pd.DataFrame(data.index, index=data.index), obs=cellmeta)
celltype = np.array(['8', 'zygote', 'early2C', 'mid2C', 'late2C', '2C', '4C', '8C', '16C','E3.5', 'E4.5', 'E5.5', 'E6.5', 'E6.75', 'ES', 'TBLC', 'L.EPSCs', 'D.EPSCs'])
exprfilter = [expr.obs[0][i] in celltype for i in range(len(expr))]
expr = expr[exprfilter]
sc.pp.filter_cells(expr, min_genes=300)
sc.pp.filter_genes(expr, min_cells=50)
expr.obs['n_counts'] = expr.X.sum(axis=1)
sc.pp.normalize_per_cell(expr, counts_per_cell_after=1e4)
sc.pp.log1p(expr)
expr.raw = expr
sc.pp.highly_variable_genes(expr, n_top_genes = 1000)
expr0 = expr[:, expr.var['highly_variable']]
sc.pp.regress_out(expr0, ['n_counts'])
sc.pp.scale(expr0, max_value=10)
expr1 = expr0[[expr0.obs.index[i] in dsample for i in range(len(expr0))]]
dsample = np.concatenate([np.random.choice(expr0[expr0.obs[0]=='8'].obs.index,size=100,replace=False),np.random.choice(expr0[expr0.obs[0]=='TBLC'].obs.index,size=100,replace=False),expr0[np.logical_and(expr0.obs[0]!='8', expr0.obs[0]!='TBLC')].obs.index])

# As for reproduce, here are cells we used by random choice
dsample = np.load('Fig2c_dsample.npy', allow_pickle=True)

expr1 = expr0[[expr0.obs.index[i] in dsample for i in range(len(expr0))]]
sc.pp.neighbors(expr1, n_neighbors=50, metric='correlation')
sc.tl.umap(expr1, min_dist=1)

fig, ax = plt.subplots(figsize=(6,5))
ax.set_xlabel('UMAP-1', fontsize=20)
ax.set_ylabel('UMAP-2', fontsize=20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', which='both', length=0)
ax.set_xticklabels([])
ax.set_yticklabels([])
for i in range(len(celltype)):
	ax.scatter(expr1.obsm['X_umap'][expr1.obs[0]==celltype[i]][:,0], expr1.obsm['X_umap'][expr1.obs[0]==celltype[i]][:,1], s=8, c=color[i], label=celltype[i], alpha=0.8)


ax.legend(markerscale=2, ncol = 1, prop={'size': 10}, bbox_to_anchor=(1,1), loc='upper left', fontsize=18)
plt.tight_layout()
plt.savefig(filepath, transparent=True, bbox_inches="tight")
plt.close()


#Fig. 2d
#For left part
data = pd.read_csv('Fig2d_left.csv', index_col=0)
#For right part
data = pd.read_csv('Fig2d_right.csv', index_col=0)

col_label = data.columns
color_dict2 = {celltype[i]:color[i] for i in range(len(celltype))}
col_color = np.array([color_dict2[col_label[i]] for i in range(len(col_label))])

cg = clustermap(exprsel, vmin= np.percentile(exprsel, 10), vmax = np.percentile(exprsel, 90), cmap='coolwarm', xticklabels=[], yticklabels = data.index, metric='cosine', figsize=(20,20), row_cluster=False, col_cluster=False, col_colors=col_color)
handles2 = [Patch(facecolor=color_dict2[name]) for name in color_dict2]
plt.legend(handles2, color_dict2, title='Cell type', bbox_to_anchor=(1, 0.5), bbox_transform=plt.gcf().transFigure, loc='right')
cg.savefig(filepath, transparent=True)






























