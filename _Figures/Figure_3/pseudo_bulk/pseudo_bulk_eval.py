#!/usr/bin/env python3
"""
Created on 2023-10-16 (Mon) 17:24:56

evaluation with pseudo-bulk from scRNA-Seq

@author: I.Azuma
"""
# %%
Base_dir = 'C:/github/LiverDeconv' # cloning repository

import sys
sys.path.append(Base_dir)
import liver_deconv
from Evaluation import evaluator

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %% deconvolution with bulk derived reference
# load data
df_mix = pd.read_csv(Base_dir + '/_Figures/Figure_3/data/pseudo_bulk/27361x1000_ann.csv',index_col=0)
fxn = lambda x : np.log2(x+1)
df_log = df_mix.applymap(fxn)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

comb_df = pd.read_pickle(Base_dir + '/_Figures/Figure_3/data/bulk/nk_fluctuate_quartile_comb2ref.pkl')
comb_df = comb_df.sort_values("NK",ascending=False)
comb = comb_df['pair'].tolist()

n = 5
score_list = []
for i in range(n):
    target_cells = ['NK']+(list(comb[i])) # FIXME: -i-1 when consider bad combinations
    use_samples = []
    for t in df_all.columns.tolist():
        if t.split("_")[0] in target_cells:
            use_samples.append(t)
    df_target = df_all[use_samples]

    # prep
    dat = liver_deconv.LiverDeconv()
    dat.set_data(df_mix=df_log, df_all=df_target)
    dat.pre_processing(do_ann=False,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
    dat.narrow_intersec()
    dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

    #  deconvolution
    dat.do_fit(method="ElasticNet",alpha=1,l1_ratio=0.05)
    res = dat.get_res()

    # evaluation
    val_df = pd.read_csv(Base_dir + '/_Figures/Figure_3/data/pseudo_bulk/sc10cell_proportion_df.csv',index_col=0).T
    # validation cell type selection
    #val_df = val_df.loc[['Neutrophil','Monocyte','B','CD8','NK','Kupffer']]
    #val_norm = val_df/(val_df.sum())

    ev = evaluator.Evaluator()
    ev.set_deconv_res(res_df=res,z_score=True)
    ev.set_validation_ref(val_df=val_df)
    ev.process_validation_ref(z_score=True)

    ev.evaluate(dec_names=[["NK"]],
                val_names=[["NK"]],title="NK",do_plot=True,simple=True,eval_all=False,dpi=100)
    score_list.append(ev.cor_res.get('NK'))
print(score_list)
# %%
"""
- Neutrophils
neu_good = [0.939, 0.9043, 0.8771, 0.9309, 0.9008]
neu_bad = [0.8575, 0.8172, 0.8284, 0.8232, 0.828]
neu_baseline = 0.926

- Monocytes
mon_good = [0.7317, 0.8167, 0.7857, 0.8096, 0.7755]
mon_bad = [0.75, 0.6622, 0.7029, 0.727, 0.6771]
mon_baseline = 0.941

- NK cells
nk_good = [0.9, 0.8949, 0.8897, 0.9048, 0.896]
nk_bad = [0.878, 0.7771, 0.8856, 0.87, 0.6468]
nk_baseline = 0.911
"""

# %% barplot
neu_list = [[0.939, 0.9043, 0.8771, 0.9309, 0.9008],[0.8575, 0.8172, 0.8284, 0.8232, 0.828]]
mon_list = [[0.7317, 0.8167, 0.7857, 0.8096, 0.7755],[0.75, 0.6622, 0.7029, 0.727, 0.6771]]
nk_list = [[0.9, 0.8949, 0.8897, 0.9048, 0.896],[0.878, 0.7771, 0.8856, 0.87, 0.6468]]

# plot
data = [neu_list,mon_list,nk_list]
cells = ['Neutrophils', 'Monocytes', 'NK']
names = ['Top 5', 'Bottom 5']
color = ['tab:blue','tab:orange']
baseline = [0.926,0.941,0.911]

fig = plt.figure(figsize=(8*4,6*1),dpi=300)
for idx,data_list in enumerate(data):
    ax = fig.add_subplot(1,4,idx+1)
    df = pd.DataFrame({'Top 5':data_list[0],'Bottom 5':data_list[1]})
    error_bar_set = dict(lw=1,capthick=1,capsize=20)
    ax.bar([0,1],df.mean(),yerr=df.std(),tick_label=df.columns,error_kw=error_bar_set,color=color,width=0.5)
    # jitter plot
    df_melt = pd.melt(df)
    sns.stripplot(x='variable', y='value', data=df_melt, jitter=True, color='black', ax = ax)
    # baseline
    plt.hlines(0.,-0.4,1.4,ls="--",color="black")
    plt.hlines(baseline[idx],-0.4,1.4,ls="-.",color="red")
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel('Correlation',fontsize=18)
    plt.xlabel('')
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    ax.set_axisbelow(True)
    ax.grid(color="#ababab",linewidth=1.0, axis='y')
    plt.title(cells[idx],fontsize=18)
    
plt.show()
# %%
