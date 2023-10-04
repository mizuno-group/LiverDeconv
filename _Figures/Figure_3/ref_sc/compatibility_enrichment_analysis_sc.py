#!/usr/bin/env python3
"""
Created on 2023-09-23 (Sat) 17:30:25

@author: I.Azuma
"""
# %% loop (Neutrophil)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

neu_df = pd.read_pickle('C:/github/LiverDeconv/_Figures/Figure_3/data/sc/neu_fluctuate_quartile_comb2ref_sc.pkl')
neu_df.index = neu_df['pair']
pairs = neu_df.index.tolist()

cell_types = ['B','CD8', 'Cholangiocyte', 'Hepatocyte', 'Kupffer', 'LSEC', 'Monocyte', 'NK','Stellate']
cell_name = ['B cells', 'CD8 T cells','Cholangiocytes','Hepatocytes','Kupffer cells', 'LSEC', 'Monocytes', 'NK cells', 'Stellate cells']
concat_res = pd.DataFrame()
for idx, cell in enumerate(cell_types):
    contain = []
    for t in pairs:
        if cell in t:
            contain.append(t)
        else:
            pass
    dat = GSEA() # generate an instance
    ref_dic = {cell:set(contain)}
    obj = neu_df[['Neutrophil']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="Neutrophil",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating Neutrophils"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures/Figure_3/enrichment_result/sc/neu_enrichmet_summary.csv')
# %% loop (Monocytes)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

neu_df = pd.read_pickle('C:/github/LiverDeconv/_Figures/Figure_3/data/sc/mon_fluctuate_quartile_comb2ref_sc.pkl')
neu_df.index = neu_df['pair']
pairs = neu_df.index.tolist()

cell_types = ['B','CD8', 'Cholangiocyte', 'Hepatocyte', 'Kupffer', 'LSEC', 'Neutrophil', 'NK','Stellate']
cell_name = ['B cells', 'CD8 T cells','Cholangiocytes','Hepatocytes','Kupffer cells', 'LSEC', 'Neutrophils', 'NK cells', 'Stellate cells']
concat_res = pd.DataFrame()
for idx, cell in enumerate(cell_types):
    contain = []
    for t in pairs:
        if cell in t:
            contain.append(t)
        else:
            pass
    dat = GSEA() # generate an instance
    ref_dic = {cell:set(contain)}
    obj = neu_df[['Monocyte']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="Monocyte",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating Monocytes"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures/Figure_3/enrichment_result/sc/mon_enrichmet_summary.csv')

# %% loop (NK)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

neu_df = pd.read_pickle('C:/github/LiverDeconv/_Figures/Figure_3/data/sc/nk_fluctuate_quartile_comb2ref_sc.pkl')
neu_df.index = neu_df['pair']
pairs = neu_df.index.tolist()

cell_types = ['B','CD8', 'Cholangiocyte', 'Hepatocyte', 'Kupffer', 'LSEC', 'Monocyte', 'Neutrophil','Stellate']
cell_name = ['B cells', 'CD8 T cells','Cholangiocytes','Hepatocytes','Kupffer cells', 'LSEC', 'Monocytes', 'Neutrophils', 'Stellate cells']
concat_res = pd.DataFrame()
for idx, cell in enumerate(cell_types):
    contain = []
    for t in pairs:
        if cell in t:
            contain.append(t)
        else:
            pass
    dat = GSEA() # generate an instance
    ref_dic = {cell:set(contain)}
    obj = neu_df[['NK']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="NK",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating NK cells"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures/Figure_3/enrichment_result/sc/nk_enrichmet_summary.csv')
# %% summary
import glob
import pandas as pd

l = glob.glob('C:/github/LiverDeconv/_Figures/Figure_3/enrichment_result/sc/*.csv')

summary_df = pd.DataFrame()
for p in l:
    cell = p.split('\\')[1].split('_')[0]
    tmp_df = pd.read_csv(p,index_col=0)
    summary_df = pd.concat([summary_df,tmp_df],axis=1)
    
import seaborn as sns
import matplotlib.pyplot as plt

fig,ax = plt.subplots(dpi=300)
sns.heatmap(summary_df,annot=True, fmt="1.3f",annot_kws={"fontsize":8}, linewidths=.5, cmap='PiYG_r')
plt.yticks(rotation=0,size=10)
plt.xticks(rotation=0,size=10)
plt.xlabel("")
plt.ylabel("")
plt.show()
# %%
