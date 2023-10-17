#!/usr/bin/env python3
"""
Created on 2023-10-17 (Tue) 17:59:11

reference cell types impact on pseudo-bulk dataset

- reference: bulk
- mix: pseudo-bulk (from scRNA-Seq)
- fix the number of whole cell types

@author: I.Azuma
"""
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
Base_dir = 'C:/github/LiverDeconv' # cloning repository
import sys
sys.path.append('C:/github/LiverDeconv')
import liver_deconv
from Evaluation import evaluator

#%% LM6
df_mix = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/fixwhole_27361x1000_ann.csv',index_col=0)
fxn = lambda x : np.log2(x+1)
df_log = df_mix.applymap(fxn)
df_all = pd.read_csv(Base_dir+'/data/processed/ref_13types.csv',index_col=0)

#target_cells = ["Neutrophil","Monocyte","B","CD8","CD4","NK"]
#target_cells = ["Neutrophil","Monocyte","B","CD8","CD4","NK","Eosinophil","Basophil","Kupffer"]
target_cells = ["Neutrophil","Monocyte","B","CD8","CD4","NK","Eosinophil","Basophil","Kupffer","Hepatocyte","Cholangiocyte","Stellate","LSEC"]

use_samples = []
for t in df_all.columns.tolist():
    if t.split("_")[0] in target_cells:
        use_samples.append(t)
df_target = df_all[use_samples]

# reference prep
dat = liver_deconv.LiverDeconv()
dat.set_data(df_mix=df_log, df_all=df_target)
dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True)
dat.narrow_intersec()
dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=True,do_plot=True)

# deconvolution
dat.do_fit(method="ElasticNet",alpha=1,l1_ratio=0.05)
res = dat.get_res()

# evaluation
val_df = pd.read_csv('D:/GdriveSymbol/notebook/dry/Deconvolution/230906_revise/results/231002_sc_pseudo/sc10cell_proportion_fixwhole.csv',index_col=0).T

ev = evaluator.Evaluator()
ev.set_deconv_res(res_df=res,z_score=True)
ev.set_validation_ref(val_df=val_df)
ev.process_validation_ref(z_score=True)

ev.evaluate(dec_names=[["Neutrophil"]],
            val_names=[["Neutrophil"]],title="Neutrophils",do_plot=True,simple=True,eval_all=False,dpi=300)
ev.evaluate(dec_names=[["Monocyte"]],
            val_names=[["Monocyte"]],title="Monocytes",do_plot=True,simple=True,eval_all=False,dpi=300)
ev.evaluate(dec_names=[["NK"]],
            val_names=[["NK"]],title="NK cells",do_plot=True,simple=True,eval_all=False,dpi=300)

# %%
lm6 = [0.921, 0.760, 0.900]
lm9 = [0.958, 0.896, 0.914]
lm13 = [0.882, 0.902, 0.923]
data = [lm6,lm9,lm13]
methods = ['6 cells', '9 cells', '13 cells']

cells = ['Neutrophils', 'Monocytes', 'NK']

hatch_list = ['/', '|', '-', '+', 'x', 'o', 'O', '.', '*']
fig,ax = plt.subplots(figsize=(10,6),dpi=300)
x = [k*len(data) for k in range(len(cells))]
for i in range(len(data)):
    v = data[i]
    x2 = [t+i*0.8 for t in x]
    plt.bar(x2,v,label=methods[i],width=0.7,hatch=hatch_list[i]*3)
plt.xticks([t+0.8 for t in x],cells,fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel('Correlation',fontsize=18)
plt.legend(loc='upper right',shadow=True,fontsize=15,bbox_to_anchor=(1.2,1))
#plt.legend(loc='upper center',shadow=True,fontsize=13,ncol=6,bbox_to_anchor=(.5, 1.12))
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=1.0, axis='y')
plt.show()
# %%
