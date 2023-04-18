# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 17:23:21 2023

@author: I.Azuma
"""

#%%
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import numpy as np
from matplotlib import gridspec
        
dat = GSEA() # generate an instance
ref,obj = dat.generate_test_data() # generate test data
dat.fit(ref) # load reference

ref['EDPWo'] = {7167, 1138, 1388, 1600, 8223, 8900}

res = dat.calc(obj,method='gsva') # analyze data of interest
res.head()

dat.plot_running(sample_name="xxxx",fterm="EDPWo",heatmap=True)

#%% Kupffer and stellate
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import numpy as np
from matplotlib import gridspec
import pandas as pd

neu_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/neu_fluctuate_quartile_comb2ref.pkl')
neu_df.index = neu_df['pair']
pairs = neu_df.index.tolist()

contain = []
for t in pairs:
    if 'Kupffer' in t:
        if 'Stellate' in t:
            contain.append(t)
    else:
        pass

ref_dic = {'Kupffer and Stellate':set(contain)}
obj = neu_df[['Neutrophil']]

obj = obj.astype(float)

dat = GSEA() # generate an instance
dat.fit(ref_dic) # load reference
res = dat.calc(obj,method='standard') # analyze data of interest
res.head()
dat.plot_running(sample_name="Neutrophil",fterm="Kupffer and Stellate",heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Kupffer and Stellate in estimating Neutrophils")

#%% loop (Neutrophil)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

neu_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/neu_fluctuate_quartile_comb2ref.pkl')
neu_df.index = neu_df['pair']
pairs = neu_df.index.tolist()

cell_types = ['B', 'Basophil', 'CD4', 'CD8', 'Cholangiocyte', 'Eosinophil', 'Hepatocyte', 'Kupffer', 'LSEC', 'Monocyte', 'NK','Stellate']
cell_name = ['B cells', 'Basohpils', 'CD4 T cells', 'CD8 T cells','Cholangiocytes', 'Eosinophils', 'Hepatocytes','Kupffer cells', 'LSEC', 'Monocytes', 'NK cells', 'Stellate cells']
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

concat_res.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/enrichment_result/neu_enrichmet_summary.csv')

#%% loop (Monocytes)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

mon_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.pkl')
mon_df.index = mon_df['pair']
pairs = mon_df.index.tolist()

cell_types = ['B', 'Basophil', 'CD4', 'CD8', 'Cholangiocyte', 'Eosinophil', 'Hepatocyte', 'Kupffer', 'LSEC', 'Neutrophil', 'NK','Stellate']
cell_name = ['B cells', 'Basohpils', 'CD4 T cells', 'CD8 T cells','Cholangiocytes', 'Eosinophils', 'Hepatocytes','Kupffer cells', 'LSEC', 'Neutrophils', 'NK cells', 'Stellate cells']
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
    obj = mon_df[['Monocyte']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="Monocyte",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating Monocytes"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/enrichment_result/mon_enrichmet_summary.csv')

#%% loop (NK)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

nk_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/nk_fluctuate_quartile_comb2ref.pkl')
nk_df.index = nk_df['pair']
pairs = nk_df.index.tolist()

cell_types = ['B', 'Basophil', 'CD4', 'CD8', 'Cholangiocyte', 'Eosinophil', 'Hepatocyte', 'Kupffer', 'LSEC', 'Neutrophil', 'Monocyte','Stellate']
cell_name = ['B cells', 'Basohpils', 'CD4 T cells', 'CD8 T cells','Cholangiocytes', 'Eosinophils', 'Hepatocytes','Kupffer cells', 'LSEC', 'Neutrophils', 'Monocytes', 'Stellate cells']
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
    obj = nk_df[['NK']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="NK",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating NK cells"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/enrichment_result/nk_enrichmet_summary.csv')

#%% loop (Eosinophil)
import sys
sys.path.append('C:/github/enan')

from enan import GSEA
import pandas as pd

eos_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/eos_fluctuate_quartile_comb2ref.pkl')
eos_df.index = eos_df['pair']
pairs = eos_df.index.tolist()

cell_types = ['B', 'Basophil', 'CD4', 'CD8', 'Cholangiocyte', 'NK', 'Hepatocyte', 'Kupffer', 'LSEC', 'Neutrophil', 'Monocyte','Stellate']
cell_name = ['B cells', 'Basohpils', 'CD4 T cells', 'CD8 T cells','Cholangiocytes', 'NK cells', 'Hepatocytes','Kupffer cells', 'LSEC', 'Neutrophils', 'Monocytes', 'Stellate cells']
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
    obj = eos_df[['Eosinophil']]

    obj = obj.astype(float)

    dat.fit(ref_dic) # load reference
    res = dat.calc(obj,method='standard') # analyze data of interest
    concat_res = pd.concat([concat_res,res])
    
    dat.plot_running(sample_name="Eosinophil",fterm=cell,heatmap=True,barcode_params={'color':'black','size':0.15},dpi=300,title="Reference including %s in estimating Eosinophils"%cell_name[idx])

concat_res.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/enrichment_result/eos_enrichmet_summary.csv')

#%% summary
import glob
import pandas as pd

l = glob.glob('C:/github/LiverDeconv/_Figures_new/Figure_4/enrichment_result/*.csv')

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