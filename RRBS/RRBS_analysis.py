import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import scanpy as sc
from itertools import cycle, islice
import matplotlib.pyplot as plt


plotmarker=['.','o','v','^','<','>','4','8','s','p','|','h','+','1','*','*','*','*', '2', '3', 'H', 'x', 'X', 'D', 'd', '_']
color4 = np.array(list(islice(cycle(['#800000', '#cc0000', '#ff0000', '#ff4d4d', '#ffc300', '#ffad60', '#293462','#ff5733', '#d9534f']), 100)))



#Fig. 3e PCA
rateb = pd.read_csv('Fig3e_data.csv', index_col=0)
pca = PCA(n_components=5)
data_reduce = pca.fit_transform(rateb.T)
color4 = np.array(list(islice(cycle(['#800000', '#cc0000', '#ff0000', '#ff4d4d', '#ffc300', '#ffad60', '#293462','#ff5733', '#d9534f']), 100)))

fig, ax = plt.subplots(figsize = (4,3))
ax.set_xlabel('PC-1', fontsize=15)
ax.set_ylabel('PC-2', fontsize=15)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', which='both', length=0)
for i,x in enumerate(['Zygote', '2C', '4C', '8C', 'E6.5', 'E7.5', 'mES2iL', '8-rep1','8-rep2']):
    ax.scatter(data_reduce[i, 0], data_reduce[i, 1], c=color4[i], s=80, edgecolors='none', alpha=0.8, label=x, marker=plotmarker[i])

ax.legend(markerscale=1, ncol = 1, prop={'size': 10}, bbox_to_anchor=(1,1), loc='upper left', fontsize=18)
plt.tight_layout()
plt.savefig(path_fig3e, transparent=True)
plt.close()


#Fig. 3f hist
samplename = np.array(['Zygote', '2C', '4C', '8C', 'E6.5', 'E7.5', 'mES-2iL', '8-rep1', '8-rep2'])
rrbsdata = rateb
fig,axes = plt.subplots(1,9,figsize=(36,4))
for i in range(len(samplename)):
	ax = axes[i]
	ax.hist(rrbsdata.iloc[:,i][rrbsdata.iloc[:,i]!='na'].dropna().astype(float), 20, facecolor=color4[i], alpha=0.75)
	ax.set_xlabel('CpG Methylation', fontsize=12)
	ax.set_ylabel('Count', fontsize=12)
	ax.set_title(rrbsdata.columns[i] + ' (n = ' + str(len(rrbsdata.iloc[:,i][rrbsdata.iloc[:,i]!='na'].dropna())) + ')')
	ax.set_xticklabels([0,0.5,1])
	ax.set_xticks([0,0.5,1])

plt.tight_layout()
plt.savefig(path_fig3f, transparent=True)
plt.close()

























