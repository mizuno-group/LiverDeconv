# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:53:49 2023

Relationship between multicollinearity and estimation performance

@author: I.Azuma
"""
Base_dir = 'C:/github/LiverDeconv' # cloning repository

import sys
sys.path.append(Base_dir)
import liver_deconv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% correlation of explanatory variables (detected DEGs)
# load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_all)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
print(final_ref.shape)
final_ref.head()
final_dic = dat.deg_dic
print(final_dic['Neutrophil'])
"""
['MMP9', 'IL1F9', 'S100A9', 'GCNT2', 'CSF3R', 'ACTA2', 'LCN2', 'IL1B', 'C5AR1', 'SLPI', 'CD33', 'MIRT2', 'ELP1', 'GSDME', 'DHRS7', 'GM19696', 'NLRP12', 'SCNN1A', 'FBXL5', 'RAB11FIP1', 'MXD1', 'IGFBP6', 'MPZL3', 'FPR2', 'TCTEX1D2', 'GM4518', 'ARRDC3', 'CDK2AP2', 'TNFRSF23', 'CLEC7A', 'PFKFB4', 'ZFP119B', 'LRG1', 'ADIPOR1', 'TLR6', 'CD44', '4833418N02RIK', 'TLR5', 'FGL2', 'GPCPD1', 'HDAC4', 'GM15832', 'SFXN5', 'PLA2G7', 'RALGPS1', 'RESF1', 'SLC35A5', 'TIRAP', 'ZFP119A', 'ANTXR2']
"""

exp_cor = final_ref.corr()

fig,ax = plt.subplots(dpi=300)
sns.heatmap(exp_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

#%% correlation of explanatory variables (detected DEGs)
# load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_all)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
print(final_ref.shape)
final_ref.head()

# gene-wide
from sklearn.preprocessing import MinMaxScaler
mm_scaler = MinMaxScaler()
mm_df = (pd.DataFrame(mm_scaler.fit_transform(final_ref.T),index=final_ref.T.index, columns=final_ref.T.columns)).T

exp_cor = mm_df.corr()

fig,ax = plt.subplots(dpi=300)
sns.heatmap(exp_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

#%% correlation of expalanatory variables (without Kupffer cells)
# load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

target_cells = ["Neutrophil","Monocyte","B","CD8","CD4","NK","Eosinophil","Basophil","Hepatocyte",'Cholangiocyte',"Stellate","LSEC"]
use_samples = []
for t in df_all.columns.tolist():
    if t.split("_")[0] in target_cells:
        use_samples.append(t)
df_target = df_all[use_samples]


dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_target)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
print(final_ref.shape)
final_ref.head()
final_dic = dat.deg_dic
print(final_dic['Neutrophil'])
"""
['MMP9', 'IL1B', 'IL1F9', 'CD300LD', 'S100A9', 'GCNT2', 'CSF3R', 'FCGR4', 'CLEC4N', 'PILRA', 'ACTA2', 'LCN2', 'PILRB2', 'SIGLECE', 'TREML4', 'SLPI', 'C5AR1', 'FPR2', 'ENTPD1', 'CD33', 'MIRT2', 'ELP1', 'GSDME', 'DHRS7', 'PILRB1', 'GM19696', 'NLRP12', 'PLA2G7', 'CD300C2', 'SCNN1A', 'MXD1', 'FBXL5', 'RAB11FIP1', 'ARRDC3', 'IGFBP6', 'MPZL3', 'GM4518', 'TCTEX1D2', 'CDK2AP2', 'CLEC7A', 'TNFRSF23', 'PFKFB4', 'WFDC17', 'LRG1', 'ADIPOR1', 'ZFP119B', 'TLR6', 'TNFAIP2', 'CD44', 'SLC11A1']
"""

# correlation
exp_cor = final_ref.corr()
fig,ax = plt.subplots(dpi=300)
sns.heatmap(exp_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

# vif
import pandas as pd
import copy
import numpy as np
from patsy import dmatrices
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

df = copy.deepcopy(final_ref)
#gather features
vif_data = pd.DataFrame()
vif_data["feature"] = df.columns

# calculating VIF for each feature
vif_data["VIF"] = [variance_inflation_factor(df.values, i) for i in range(len(df.columns))]
  
print(vif_data)


fig,ax = plt.subplots(dpi=300)
plt.bar([i for i  in range(len(vif_data))],vif_data['VIF'])
plt.yticks(rotation=0)
plt.xticks([i for i  in range(len(vif_data))],vif_data['feature'],rotation=90)
plt.xlabel("")
plt.ylabel("Variance Inflation Factor (VIF)")

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.show()

# common degs
set1 = set(['MMP9', 'IL1F9', 'S100A9', 'GCNT2', 'CSF3R', 'ACTA2', 'LCN2', 'IL1B', 'C5AR1', 'SLPI', 'CD33', 'MIRT2', 'ELP1', 'GSDME', 'DHRS7', 'GM19696', 'NLRP12', 'SCNN1A', 'FBXL5', 'RAB11FIP1', 'MXD1', 'IGFBP6', 'MPZL3', 'FPR2', 'TCTEX1D2', 'GM4518', 'ARRDC3', 'CDK2AP2', 'TNFRSF23', 'CLEC7A', 'PFKFB4', 'ZFP119B', 'LRG1', 'ADIPOR1', 'TLR6', 'CD44', '4833418N02RIK', 'TLR5', 'FGL2', 'GPCPD1', 'HDAC4', 'GM15832', 'SFXN5', 'PLA2G7', 'RALGPS1', 'RESF1', 'SLC35A5', 'TIRAP', 'ZFP119A', 'ANTXR2'])
set2 = set(['MMP9', 'IL1B', 'IL1F9', 'CD300LD', 'S100A9', 'GCNT2', 'CSF3R', 'FCGR4', 'CLEC4N', 'PILRA', 'ACTA2', 'LCN2', 'PILRB2', 'SIGLECE', 'TREML4', 'SLPI', 'C5AR1', 'FPR2', 'ENTPD1', 'CD33', 'MIRT2', 'ELP1', 'GSDME', 'DHRS7', 'PILRB1', 'GM19696', 'NLRP12', 'PLA2G7', 'CD300C2', 'SCNN1A', 'MXD1', 'FBXL5', 'RAB11FIP1', 'ARRDC3', 'IGFBP6', 'MPZL3', 'GM4518', 'TCTEX1D2', 'CDK2AP2', 'CLEC7A', 'TNFRSF23', 'PFKFB4', 'WFDC17', 'LRG1', 'ADIPOR1', 'ZFP119B', 'TLR6', 'TNFAIP2', 'CD44', 'SLC11A1'])
common = set1 & set2

without_kup = set2 - common
kup_set = set(['CLEC4F', 'C1QB', 'CD207', 'FOLR2', 'VSIG4', 'C1QC', 'CD5L', 'CXCL13', 'TIMD4', 'B3GALNT1', 'C1QA', 'CCL24', 'EPB41L3', 'AIF1', 'FABP7', 'MARCO', 'TANC2', 'HSPA1A', 'ZDHHC14', 'BATF3', 'CEBPB', 'MAFB', 'SPIC', 'IL18BP', 'DMPK', 'FTL1-PS1', 'SPIRE1', 'CD14', 'HSPA1B', 'NR1H3', 'GNGT2', 'CTSC', 'SIGLEC1', 'CFP', 'LPL', 'TMEM26', 'NEURL1A', 'CSF1R', 'SERPINB8', 'ENGASE', 'ZFP90', 'KLF13', 'IRF7', 'GPX3', 'ABHD12', 'CD68', 'CREG1', 'WFDC17', 'HPGD', 'NTPCR'])

kup_common = without_kup & kup_set

#%% correlation of explanatory variables (whole profile)
df = df_all.T
df['cell'] = [t.split('_')[0] for t in df.index]
med_df = df.groupby('cell').median()

med_cor = med_df.T.corr()

fig,ax = plt.subplots(dpi=300)
sns.heatmap(med_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

#%% variance inflation factor (detected DEGs)
# load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_all)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
print(final_ref.shape)
final_ref.head()

import pandas as pd
import copy
import numpy as np
from patsy import dmatrices
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

df = copy.deepcopy(final_ref)
#gather features
vif_data = pd.DataFrame()
vif_data["feature"] = df.columns

# calculating VIF for each feature
vif_data["VIF"] = [variance_inflation_factor(df.values, i) for i in range(len(df.columns))]
  
print(vif_data)


fig,ax = plt.subplots(dpi=300)
plt.bar([i for i  in range(len(vif_data))],vif_data['VIF'])
plt.yticks(rotation=0)
plt.xticks([i for i  in range(len(vif_data))],vif_data['feature'],rotation=90)
plt.xlabel("")
plt.ylabel("Variance Inflation Factor (VIF)")

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.show()

#%% https://zenn.dev/wsuzume/articles/fab5000ed8bafd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

def collinearity_check(Xc, model=None, alpha=1.0, emph=False):
    """
    Parameters
    ----------
    Xc : np.ndarray(m, n)
        Input data (each data stored in a row).
    model
        Regression model (default Ridge).
    alpha : float
        Hyper parameter of Ridge (default 1.0),
        ignored if model is not None.
    emph : bool
        Emphasize the result by R2 score or not.

    Returns
    -------
    rc : np.ndarray(n, n)
        Regression coefficient, emphasized by R2 score if emph is True.
    scores : np.ndarray(n)
        R2 scores.
    """
    
    if model is None:
        model = Ridge(alpha=alpha)

    m, n = Xc.shape
    if n < 2:
        raise ValueError()

    # 戻り値
    rc = np.empty((n, n)) # 回帰係数
    scores = np.empty(n) # R2スコア

    X = np.copy(Xc)
    for i in range(n):
        y = np.copy(X[:,i]) # 自分を他の変数で回帰させる
        X[:,i] = 1 # 自分は多重共線性の定数項に対応させる

        model.fit(X, y)
        y_calc = model.predict(X)

        score = r2_score(y, y_calc)
        if score < 0:
            # R2 スコアが 0 以下なら線形性なしとみなす
            scores[i] = 0
            rc[i] = 0
        else:
            scores[i] = score
            if emph:
                # 係数が大きくても R2 スコアが 0 に近ければ 0 になるように加工
                rc[i] = model.coef_ * score
            else:
                rc[i] = model.coef_

        X[:,i] = y
    
    return rc, scores

from sklearn.preprocessing import StandardScaler

Xsc = StandardScaler().fit_transform(final_ref)
cells = final_ref.columns.tolist()

# visualize
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

cc, cc_scores = collinearity_check(Xsc, emph=True)
plt.figure(figsize=(10,8))
sns.heatmap(pd.DataFrame(cc), vmin=-1, vmax=1, cmap='bwr', cbar=True)

plt.yticks([i+0.5 for i  in range(len(cc))],cells,rotation=0)
plt.xticks([i+0.5 for i  in range(len(cc))],cells,rotation=90)
plt.ylabel('target variable', fontsize=16)
plt.xlabel('used variables to explain the target', fontsize=16)
plt.show()


from matplotlib import pyplot as plt

plt.figure(figsize=(8,5))
plt.plot(cc_scores)
#plt.xticks(np.arange(0, len(cc_scores), 1))
plt.xticks([i for i  in range(len(cc))],cells,rotation=90)
plt.ylabel('R2 score', fontsize=16)
plt.xlabel('selected variables', fontsize=16)
plt.show()

#%%
# load data
df_mix = pd.read_csv(Base_dir+'/data/processed/mix_processed.csv',index_col=0)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

# neutrophil
dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_all)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
deg_dic = dat.deg_dic

neu_deg = deg_dic['Neutrophil']
wo_neu_df = final_ref.loc[sorted(list(set(final_ref.index) - set(neu_deg)))]

# correlation
exp_cor = wo_neu_df.corr()
fig,ax = plt.subplots(dpi=300)
sns.heatmap(exp_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

cells = wo_neu_df.columns.tolist()
cor_list = []
for cell in cells:
    target_deg = deg_dic.get(cell)
    wo_neu_target_df = wo_neu_df.loc[sorted(list(set(wo_neu_df.index) - set(target_deg)))]
    cor = wo_neu_target_df.corr()['Neutrophil'][cell]
    cor_list.append(cor)
cor_dic = dict(zip(cells,cor_list))


# eosinophil
dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_mix, df_all=df_all)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

final_ref = dat.final_ref
deg_dic = dat.deg_dic

neu_deg = deg_dic['Eosinophil']
wo_neu_df = final_ref.loc[sorted(list(set(final_ref.index) - set(neu_deg)))]

# correlation
exp_cor = wo_neu_df.corr()
fig,ax = plt.subplots(dpi=300)
sns.heatmap(exp_cor,annot=True, fmt="1.2f",annot_kws={"fontsize":5}, linewidths=.5)
plt.yticks(rotation=0,size=5)
plt.xticks(rotation=90,size=5)
plt.xlabel("")
plt.ylabel("")
plt.show()

cells = wo_neu_df.columns.tolist()
cor_list = []
for cell in cells:
    target_deg = deg_dic.get(cell)
    wo_neu_target_df = wo_neu_df.loc[sorted(list(set(wo_neu_df.index) - set(target_deg)))]
    cor = wo_neu_target_df.corr()['Eosinophil'][cell]
    cor_list.append(cor)
cor_dic = dict(zip(cells,cor_list))

#%%
neu_df = pd.read_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/data/neu_fluctuate_quartile_comb2ref.csv',index_col=0)
