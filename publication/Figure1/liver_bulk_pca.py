# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:27:52 2022

PCA on the obtained liver tissue gene expression profile

@author: I.Azuma
"""
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

Base_dir = 'C:/github/LiverDeconv' # cloning repository

import sys
sys.path.append(Base_dir)
import liver_deconv as ld
from _utils import processing as pc

#%% load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
ans_df = pd.read_csv('C:/github/LiverDeconv/Data/processed/facs_true_population.csv',index_col=0).T
common = sorted(set(df_mix.columns) & set(ans_df.columns))
remove_list=["CIV_7","CIV_8","CIP_7","CIP_8"]
target_mix = df_mix.drop(columns=remove_list)
#%% normalize with Ctrl expression
print(target_mix.columns.tolist())
ctrl_df = target_mix[['Ctrl_1', 'Ctrl_10', 'Ctrl_12', 'Ctrl_15', 'Ctrl_16', 'Ctrl_17', 'Ctrl_18', 'Ctrl_2', 'Ctrl_3', 'Ctrl_4', 'Ctrl_7', 'Ctrl_8', 'Ctrl_9']]
ctrl_m = ctrl_df.T.mean()
ctrl_v = ctrl_df.T.var()

norm_df = (((target_mix.T-ctrl_m))/np.sqrt(ctrl_v)).T # normalization 
print(norm_df[['Ctrl_1', 'Ctrl_10', 'Ctrl_12', 'Ctrl_15', 'Ctrl_16', 'Ctrl_17', 'Ctrl_18', 'Ctrl_2', 'Ctrl_3', 'Ctrl_4', 'Ctrl_7', 'Ctrl_8', 'Ctrl_9']].T.sum())
print(norm_df[['Ctrl_1', 'Ctrl_10', 'Ctrl_12', 'Ctrl_15', 'Ctrl_16', 'Ctrl_17', 'Ctrl_18', 'Ctrl_2', 'Ctrl_3', 'Ctrl_4', 'Ctrl_7', 'Ctrl_8', 'Ctrl_9']].T.var())

#%% PCA (3D)
from mpl_toolkits.mplot3d import Axes3D

norm_df = norm_df.replace([np.inf, -np.inf], np.nan)
final_df = norm_df.dropna().T

target_names = [t.split('_')[0] for t in target_mix.columns.tolist()]
target_dic = dict(zip(list(set(target_names)), [i for i in range(len(set(target_names)))]))
targets = [target_dic.get(k) for k in target_names]

# perform PCA
pca3 = PCA(n_components=3)
pca3.fit(final_df)
transformed = pca3.fit_transform(final_df)

fig = plt.figure(figsize = (10, 8), dpi=300)
ax = Axes3D(fig)
for label in np.unique(targets):
    p = ax.scatter(transformed[targets == label, 0],
                   transformed[targets == label, 1],
                   transformed[targets == label, 2],
                   marker = 'o', s = 200, label=np.unique(target_names)[label],alpha=0.8,linewidths=0.8,edgecolor='k')
ax.set_xlabel("First Principal Component", fontsize=14)
ax.set_ylabel("Second Principal Component", fontsize=14)
ax.set_zlabel("Third Principal Component", fontsize=14)
plt.legend(shadow=True,loc='right')
ax.view_init(elev=15, azim=60)
plt.title('3D PCA',fontsize = 21)
plt.show()

#%% PCA (2D)
norm_df = norm_df.replace([np.inf, -np.inf], np.nan)
final_df = norm_df.dropna().T

pca3 = PCA()
pca3.fit(final_df)
transformed = pca3.fit_transform(final_df)

pd.DataFrame(transformed, columns=["PC{}".format(x + 1) for x in range(len(transformed))]).head()

fig,ax = plt.subplots(figsize = (10, 8), dpi=300)

for label in np.unique(targets):
    p = ax.scatter(transformed[targets == label, 0],
                   transformed[targets == label, 1],
                   marker = 'o', s = 200, label=np.unique(target_names)[label],alpha=0.8,linewidths=0.8,edgecolor='k')

ax.set_xlabel("First Principal Component", fontsize=14)
ax.set_ylabel("Second Principal Component", fontsize=14)
plt.legend(shadow=True)
plt.title('2D PCA',fontsize = 21)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.show()

#%% heatmap
import seaborn as sns
sns.clustermap(target_mix,col_cluster=False)
plt.show()

z_target_mix = pc.standardz_sample(target_mix.T) # sample-wide z score
z_target_mix.sum()

sns.clustermap(z_target_mix.T,col_cluster=False)
plt.show()

#%% DEG definition
#df_mix = pd.read_csv('C:/github/LiverDeconv/data/input/mix_processed.csv',index_col=0)
#df_all = pd.read_csv('C:/github/LiverDeconv/data/input/ref_13types.csv',index_col=0)

dat = ld.LiverDeconv()
dat.set_data(df_mix=target_mix, df_all=target_mix)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=False,do_quantile=False,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.2,log2=False,verbose=True)
final_ref = dat.final_ref
deg_dic = dat.deg_dic

# plot
sns.clustermap(final_ref,col_cluster=False)
plt.show()

z_final_ref = pc.standardz_sample(final_ref.T) # sample-wide
z_final_ref.sum()

sns.clustermap(z_final_ref.T,col_cluster=False)
plt.show()

#%% Treatment correlation
treat_df = target_mix.T
treat_df['Treat'] = [t.split('_')[0] for t in treat_df.index]
treat_avg = treat_df.groupby('Treat').mean()
treat_cor = treat_avg.T.corr()

import seaborn as sns
fig,ax = plt.subplots(dpi=300)
sns.heatmap(treat_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":8}, linewidths=.5)
plt.yticks(rotation=0)
plt.xlabel("")
plt.ylabel("")
plt.show()