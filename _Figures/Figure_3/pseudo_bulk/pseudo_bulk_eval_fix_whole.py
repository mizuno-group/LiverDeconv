#!/usr/bin/env python3
"""
Created on 2023-10-17 (Tue) 17:11:06

evaluate with pseudo-bulk from scRNA-Seq

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
cell_name = "NK"
# load data
df_mix = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/fixwhole_27361x1000_ann.csv',index_col=0)
fxn = lambda x : np.log2(x+1)
df_log = df_mix.applymap(fxn)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

comb_df = pd.read_pickle(Base_dir + '/_Figures/Figure_3/data/bulk/nk_fluctuate_quartile_comb2ref.pkl')
comb_df = comb_df.sort_values(cell_name,ascending=False)
comb = comb_df['pair'].tolist()

n = 10
score_list = []
for i in range(n):
    target_cells = [cell_name]+(list(comb[])) # FIXME: -i-1 when consider bad combinations
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
    val_df = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/sc10cell_proportion_fixwhole.csv',index_col=0).T
    """
    # validation cell type selection
    val_df = val_df.loc[['Neutrophil','Monocyte','B','CD8','NK','Kupffer']]
    val_norm = val_df/(val_df.sum())
    """

    ev = evaluator.Evaluator()
    ev.set_deconv_res(res_df=res,z_score=True)
    ev.set_validation_ref(val_df=val_df)
    ev.process_validation_ref(z_score=True)

    ev.evaluate(dec_names=[[cell_name]],
                val_names=[[cell_name]],title=cell_name,do_plot=True,simple=True,eval_all=False,dpi=100)
    score_list.append(ev.cor_res.get(cell_name))
print(score_list)

"""
neu_good: [0.9252, 0.8983, 0.821, 0.927, 0.8563, 0.9426, 0.862, 0.9165, 0.9369, 0.9462]
neu_bad: [0.8814, 0.777, 0.7847, 0.8078, 0.8182, 0.8476, 0.8473, 0.8292, 0.8365, 0.7905]
neu base: 0.882

mon_good: [0.7514, 0.8675, 0.8352, 0.8116, 0.8286, 0.732, 0.7682, 0.8426, 0.8263, 0.7587]
mon_bad: [0.8053, 0.6648, 0.8033, 0.8173, 0.7425, 0.9161, 0.9023, 0.8778, 0.8837, 0.896]
mon base: 0.902

nk_good: [0.9047, 0.8986, 0.8936, 0.9049, 0.9074, 0.9069, 0.9004, 0.9054, 0.9021, 0.9009]
nk_bad: [0.9026, 0.8617, 0.9061, 0.8996, 0.684, 0.8987, 0.8873, 0.7669, 0.7778, 0.7544]
nk base: 0.923

"""

# %% baseline
# load data
df_mix = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/fixwhole_27361x1000_ann.csv',index_col=0)
fxn = lambda x : np.log2(x+1)
df_log = df_mix.applymap(fxn)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

comb_df = pd.read_pickle('C:/github/LiverDeconv/_Figures/Figure_3/data/bulk/neu_fluctuate_quartile_comb2ref.pkl')
comb_df = comb_df.sort_values("Neutrophil",ascending=False)
comb = comb_df['pair'].tolist()

# prep
dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_log, df_all=df_all)
dat.pre_processing(do_ann=False,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

#  deconvolution
dat.do_fit(method="ElasticNet",alpha=1,l1_ratio=0.05)
res = dat.get_res()

# evaluation
val_df = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/sc10cell_proportion_fixwhole.csv',index_col=0).T
# validation cell type selection
"""
val_df = val_df.loc[['Neutrophil','Monocyte','B','CD8','NK','Kupffer']]
val_norm = val_df/(val_df.sum()) # normalize with the selected cell types (FACS measured scinario)
"""

ev = evaluator.Evaluator()
ev.set_deconv_res(res_df=res,z_score=True)
ev.set_validation_ref(val_df=val_df)
ev.process_validation_ref(z_score=True)

ev.evaluate(dec_names=[["Neutrophil"]],
            val_names=[["Neutrophil"]],title="Neutrophils",do_plot=True,simple=True,eval_all=False,dpi=100)
ev.evaluate(dec_names=[["Monocyte"]],
            val_names=[["Monocyte"]],title="Monocytes",do_plot=True,simple=True,eval_all=False,dpi=100)
ev.evaluate(dec_names=[["NK"]],
            val_names=[["NK"]],title="NK cells",do_plot=True,simple=True,eval_all=False,dpi=100)


# %% barplot
neu_list = [[0.9252, 0.8983, 0.821, 0.927, 0.8563, 0.9426, 0.862, 0.9165, 0.9369, 0.9462],[0.8814, 0.777, 0.7847, 0.8078, 0.8182, 0.8476, 0.8473, 0.8292, 0.8365, 0.7905]]
mon_list = [[0.7514, 0.8675, 0.8352, 0.8116, 0.8286, 0.732, 0.7682, 0.8426, 0.8263, 0.7587],[0.8053, 0.6648, 0.8033, 0.8173, 0.7425, 0.9161, 0.9023, 0.8778, 0.8837, 0.896]]
nk_list = [[0.9047, 0.8986, 0.8936, 0.9049, 0.9074, 0.9069, 0.9004, 0.9054, 0.9021, 0.9009],[0.9026, 0.8617, 0.9061, 0.8996, 0.684, 0.8987, 0.8873, 0.7669, 0.7778, 0.7544]]

# plot
data = [neu_list,mon_list,nk_list]
cells = ['Neutrophils', 'Monocytes', 'NK']
names = ['Top 10', 'Bottom 10']
color = ['tab:blue','tab:orange']
baseline = [0.882,0.902,0.923]

fig = plt.figure(figsize=(8*4,6*1),dpi=300)
for idx,data_list in enumerate(data):
    ax = fig.add_subplot(1,4,idx+1)
    df = pd.DataFrame({'Top 10':data_list[0],'Bottom 10':data_list[1]})
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
