# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 21:33:54 2023

@author: I.Azuma
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
Base_dir = 'C:/github/LiverDeconv' # cloning repository
import sys
sys.path.append('C:/github/LiverDeconv')
from Optimization import ctrl_norm_comb2ref_optim as cnco

#%%
Base_dir = 'C:/github/LiverDeconv' # cloning repository
# load data
mix_df = pd.read_csv(Base_dir + '/data/processed/mix_processed.csv',index_col=0)
raw_ref = pd.read_csv(Base_dir + '/data/processed/ref_13types.csv',index_col=0)
val_df = pd.read_csv(Base_dir + '/data/processed/facs_true_population.csv',index_col=0)
sig_df1 = pd.read_csv('C:/github/LiverDeconv/Data/processed/503_13_signature_ref.csv',index_col=0)

optim = cnco.Ctrlnorm_comb2ref()
optim.set_datasets(mix_df=mix_df,val_df=val_df,ignore_samples=['CIV_1','CIV_4','CIV_7','CIV_8','CIP_1','CIP_2','CIP_7','CIP_8'])
optim.set_target(target="Eosinophil")
optim.select_samples(method="outerquartile",ctrl="Ctrl",threshold=1.5,q=[0, 0.25, 0.5, 0.75, 1],ignore_ctrl=True)
optim.optim_base(base_sig_df=sig_df1,dec_names=[["Eosinophil"]],val_names=[["Eosinophil"]],do_plot=True)
optim.optim_loop(raw_ref=raw_ref,candi=['B','Basophil','CD4','CD8','Cholangiocyte','Eosinophil','Hepatocyte','Kupffer',
                                            'LSEC','Monocyte','NK','Neutrophil','Stellate'],
                 target=["Eosinophil"],dec_names=[["Eosinophil"]],val_names=[["Eosinophil"]],do_plot=False)
cor_df = optim.cor_df

cor_df.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.csv')
pd.to_pickle(cor_df,'C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')

#%% 230209 add all combination
import itertools
cor_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
add_df = pd.DataFrame([0.131,('B','Basophil','CD4','CD8','Cholangiocyte','Eosinophil','Hepatocyte','Kupffer',
                                            'LSEC','Monocyte','NK','Neutrophil','Stellate')],index=cor_df.columns.tolist()).T
update_df = pd.concat([cor_df,add_df])
pd.to_pickle(update_df,'C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
update_df.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.csv')

#%% Good
# top 10 good combinations
cor_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
# bottom 10 good combinations
#eos_res = pd.read_csv(Base_dir+'/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.csv',index_col=0)
eos_res = pd.read_pickle(Base_dir+'/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
eos_res = eos_res.sort_values('Eosinophil',ascending=True)
print(eos_res.head())
bad_combination = eos_res['pair'].tolist()[0:10]

# load data
mix_df = pd.read_csv(Base_dir + '/data/processed/mix_processed.csv',index_col=0)
raw_ref = pd.read_csv(Base_dir + '/data/processed/ref_13types.csv',index_col=0)
val_df = pd.read_csv(Base_dir + '/data/processed/facs_true_population.csv',index_col=0)

# each result of good combinations
ensemble = []
for b in bad_combination:
    optim = cnco.Ctrlnorm_comb2ref()
    optim.set_datasets(mix_df=mix_df,val_df=val_df,ignore_samples=['CIV_1','CIV_4','CIV_7','CIV_8','CIP_1','CIP_2','CIP_7','CIP_8'])
    optim.set_target(target="Eosinophil")
    optim.select_samples(method="outerquartile",ctrl="Ctrl",threshold=1.5,q=[0, 0.25, 0.5, 0.75, 1],ignore_ctrl=True)
    optim.optimized_single(raw_ref=raw_ref,combination=b,target=["Eosinophil"],dec_names=[["Eosinophil"]],val_names=[["Eosinophil"]],do_plot=True,dpi=100)
    res = optim.final_dec
    val = optim.final_val
    ensemble.append(res["Eosinophil"].tolist())

# ensemble the estimated values for each iterations
ensemble_df = pd.DataFrame(ensemble).T
ensemble_df.index = res.index.tolist()
x = optim.final_val["Eosinophil"].tolist()
pd.Series(ensemble[0]).corr(pd.Series(x))
x_df = pd.DataFrame(x)
x_df.index = res.index.tolist()

#%% plot
# Plotting functions
dec_min = 0
dec_max = 0
drugs = sorted(list(set([t.split("_")[0] for t in ensemble_df.index.tolist()])))
colorlist = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray"]
colordic = dict(zip(drugs,colorlist))

handles = []
fig,ax = plt.subplots(figsize=(5,5),dpi=300)
for i,d in enumerate(drugs):
    res1 = ensemble_df.filter(regex="^"+d+"_",axis=0).T
    res2 = x_df.filter(regex="^"+d+"_",axis=0).T
    std_df = res1.std().T
    
    #tmp_cor = round(np.corrcoef(res1,res2)[0][1],3)
    for n,j in enumerate(res1.columns.tolist()):
        plt.errorbar(x=res1[j].mean(),y=res2[j].tolist(),xerr=std_df[j],capsize=10,marker='_',fmt='o',ms=8,mew=1,color=colordic.get(d),zorder=-1)
        sc = plt.scatter(res1[j],[res2[j]]*len(res1),label=d,alpha=0.95,color=colordic.get(d),s=100,linewidths=0.4,edgecolor='k')
        if n == 0:
            handles.append(sc)
    
    if res1.min().min()<dec_min:
        dec_min = res1.min().min()
    if res1.max().max()>dec_max:
        dec_max = res1.max().max()
plt.plot([dec_min,dec_max],[dec_min,dec_max],linewidth=2,color='black',linestyle='dashed',zorder=-1)
plt.legend(handles,drugs,shadow=True,loc='upper left', borderaxespad=0, bbox_to_anchor=(1, 1))
plt.xlabel("Deconvolution estimated value",fontsize=18)
plt.ylabel("FACS value",fontsize=18)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.title("Non-optimized for Eosinophils",fontsize=18)
plt.show()
#%% GSE111818 optimization
import copy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('C:/github/LiverDeconv')
import liver_deconv as ld
from _utils import processing as pc
from _utils import plot4deconv

#%%
Base_dir = 'C:/github/LiverDeconv' # cloning repository
# load data
mix_df = pd.read_csv('C:/github/LiverDeconv/_Figures_new/Figure_5/data/gse111828_exp_matrix.csv',index_col=0)
ref = pd.read_csv(Base_dir + '/data/processed/ref_13types.csv',index_col=0)

comb_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
comb_df = comb_df.sort_values("Eosinophil",ascending=False)
comb = comb_df['pair'].tolist()

#%% good mean vs bad mean
n = 10
# good
good_summary = pd.DataFrame()
good_size = []
for i in range(n):
    use_cell = ['Eosinophil']
    target_cells = copy.deepcopy(use_cell)
    #target_cells.extend([t.split("'")[0] for t in comb[i][2:-2].split("', '")])
    target_cells.extend(list(comb[i]))
    good_size.append(len(target_cells))
    # select use sample
    use_sample = []
    for t in ref.columns.tolist():
        if t.split("_")[0] in target_cells:
            use_sample.append(t)
    use_ref = ref[use_sample]
    
    # define reference
    dat = ld.LiverDeconv()
    dat.set_data(df_mix=mix_df,df_all=use_ref)
    dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
    dat.narrow_intersec()
    dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=False)
    
    df_mix = dat.df_mix
    df_all = dat.df_all
    
    final_ref = dat.final_ref
    
    dat.do_fit(file_dat=mix_df,file_ref=final_ref,max_iter=1e6,number_of_repeats=3,alpha=1,l1_ratio=0.05)
    res = dat.get_res()
    z_good = pc.standardz_sample(res)
    z_good.index = [i.split(" ")[0] for i in z_good.index]
    ctrl_mean = z_good.loc["0hr"]["Eosinophil"].mean()
    fxn = lambda x : x-ctrl_mean
    z_good['Eosinophil'] = z_good['Eosinophil'].apply(fxn)
    good_summary[i] = z_good['Eosinophil']

good_summary['name'] = good_summary.index.tolist()
good_melt = pd.melt(good_summary.groupby('name').mean().T)
good_melt['status']=['good']*len(good_melt)
good_melt.columns = ['variable','value','status']

good_min = min(good_size)
pd.to_pickle(good_melt,'C:/github/LiverDeconv/_Figures_new/Figure_5/data/eos_good_10_ave.pkl')

#%%
# bad
bad_summary = pd.DataFrame()
c = 0
bad_size = []
for i in range(100):
    if c == n:
        break
    use_cell = ['Eosinophil']
    target_cells = copy.deepcopy(use_cell)
    #target_cells.extend([t.split("'")[0] for t in comb[-i][2:-2].split("', '")]) # -i th
    target_cells.extend(list(comb[-i]))
    if len(target_cells)<good_min:
        # skip if the target cell combination size is maller than good one
        continue
    bad_size.append(len(target_cells))
    # select use sample
    use_sample = []
    for t in ref.columns.tolist():
        if t.split("_")[0] in target_cells:
            use_sample.append(t)
    use_ref = ref[use_sample]
    
    # define reference
    dat = ld.LiverDeconv()
    dat.set_data(df_mix=mix_df,df_all=use_ref)
    dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
    dat.narrow_intersec()
    dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=False)
    
    df_mix = dat.df_mix
    df_all = dat.df_all
    
    final_ref = dat.final_ref
    
    dat.do_fit(file_dat=mix_df,file_ref=final_ref,max_iter=1e6,number_of_repeats=3,alpha=1,l1_ratio=0.05)
    res = dat.get_res()
    z_bad = pc.standardz_sample(res)
    z_bad.index = [i.split(" ")[0] for i in z_bad.index]
    ctrl_mean = z_bad.loc["0hr"]["Eosinophil"].mean()
    fxn = lambda x : x-ctrl_mean
    z_bad['Eosinophil'] = z_bad['Eosinophil'].apply(fxn)
    bad_summary[i] = z_bad['Eosinophil']
    c += 1
bad_summary['name'] = bad_summary.index.tolist()
bad_melt = pd.melt(bad_summary.groupby('name').mean().T)
bad_melt['status']=['bad']*len(bad_melt)
bad_melt.columns = ['variable','value','status']

pd.to_pickle(bad_melt,'C:/github/LiverDeconv/_Figures_new/Figure_5/data/eos_bad_10_ave.pkl')
#%%
# plot
merge_melt = pd.concat([good_melt,bad_melt])
sort_index = ['12hr','24hr','36hr','48hr','72hr']
merge_melt = merge_melt[merge_melt['variable'].isin(sort_index)]
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1,1,1) # identify the plotting 
# violin
sns.violinplot(x='variable', y='value', data=merge_melt, hue="status", inner=None, cut=0, scale="count",linewidth=1.5)
plt.setp(ax.collections, alpha=.6)
# strip
sns.stripplot(x='variable', y='value', data=merge_melt, hue="status", jitter=True, linewidth=1.5, size=9, dodge=True)
handles = ax.legend_.legendHandles
plt.legend(handles[2:4],["Top %d average"% n,"Bottom %d average"% n])
# ctrl base line
plt.scatter(x=[-1],y=[0],marker="s",color='firebrick',s=100)
plt.hlines(y=0,xmin=-1.5,xmax=5-0.5,color="red",ls="dashed",linewidth=2,zorder=-2)
# mean value plotting
good_means = good_melt.groupby('variable').mean()['value'].tolist()
bad_means = bad_melt.groupby('variable').mean()['value'].tolist()
plt.plot([i-1 for i in range(6)],good_means,linewidth=2,alpha=0.9,zorder=-3)
plt.plot([i-1 for i in range(6)],bad_means,ls='-.',linewidth=2,alpha=0.9,zorder=-3)
plt.xticks([i-1 for i in range(6)],['0hr','12hr','24hr','36hr','48hr','72hr'])
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
plt.title("Eosinophils trafficking estimation")
plt.xlabel("Hours post acetaminophen")
plt.ylabel("deconvolution values")
plt.show()

#%% FACS
good_melt = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_5/data/eos_good_10_ave.pkl')
bad_melt = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_5/data/eos_bad_10_ave.pkl')
# 1st
# add ctrl information
facs = pd.read_csv('C:/github/LiverDeconv/_Figures_new/Figure_5/data/FACS_1st.csv',index_col=0)

cell = 'Eosinophil'
df = facs[[cell]].T
for i in range(4):
    df['APAP#0hr_'+str(i)] = [1]
df = df.T
df.index = [t.split('_')[0] for t in df.index]

z_df = pc.standardz_sample(df)
ctrl_mean = z_df.loc["APAP#0hr"]["Eosinophil"].mean()
fxn = lambda x : x-ctrl_mean
z_df['Eosinophil'] = z_df['Eosinophil'].apply(fxn)

z_df['variable'] = [t.split("#")[1] for t in z_df.index.tolist()]
z_df.columns = ['value','variable']
z_df['status'] = ['facs']*len(z_df)
print(z_df)

# merge and plot
sort_index = ['12hr','24hr','48hr']
good_target = good_melt[good_melt['variable'].isin(sort_index)]
bad_target = bad_melt[bad_melt['variable'].isin(sort_index)]
facs_target = z_df[z_df['variable'].isin(sort_index)]

merge_melt = pd.concat([good_target,bad_target,facs_target])

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1,1,1) # identify the plotting 
# violin
sns.violinplot(x='variable', y='value', data=merge_melt, hue="status", inner=None, cut=0, scale="count",linewidth=1.5)
plt.setp(ax.collections, alpha=.6)
# strip
sns.stripplot(x='variable', y='value', data=merge_melt, hue="status", jitter=True, linewidth=1.5, size=10, dodge=True)
handles = ax.legend_.legendHandles
plt.legend(handles[3:6],["Optimized","Not-optimized","FACS (True)"],fontsize=12,shadow=True)

# ctrl base line
plt.scatter(x=[-1],y=[0],marker="s",color='firebrick',s=100)
plt.hlines(y=0,xmin=-1.5,xmax=3-0.5,color="red",ls="dashed",linewidth=2,zorder=-2)
# mean value plotting
good_means = good_target.groupby('variable').mean()['value'].tolist()
bad_means = bad_target.groupby('variable').mean()['value'].tolist()
facs_means = facs_target.groupby('variable').mean()['value'].tolist()
plt.plot([i-1.3 for i in range(4)],[0]+good_means,ls='--',linewidth=2,alpha=0.9,zorder=-3)
plt.plot([i-1 for i in range(4)],[0]+bad_means,ls='-.',linewidth=2,alpha=0.9,zorder=-3)
plt.plot([i-0.7 for i in range(4)],[0]+facs_means,ls='-',linewidth=2,alpha=0.9,zorder=-3)
plt.xticks([i-1 for i in range(4)],['0hr','12hr','24hr','48hr'])
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=0.5)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
plt.title("Eosinophils trafficking estimation",fontsize=18)
plt.xlabel("Hours post acetaminophen",fontsize=18)
plt.ylabel("Deconvolution values",fontsize=18)
plt.show()