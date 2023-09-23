# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 00:15:02 2023

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
optim.set_target(target="Monocyte")
optim.select_samples(method="outerquartile",ctrl="Ctrl",threshold=1.5,q=[0, 0.25, 0.5, 0.75, 1],ignore_ctrl=True)
optim.optim_base(base_sig_df=sig_df1,dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=True)
optim.optim_loop(raw_ref=raw_ref,candi=['B','Basophil','CD4','CD8','Cholangiocyte','Eosinophil','Hepatocyte','Kupffer',
                                            'LSEC','Monocyte','NK','Neutrophil','Stellate'],
                 target=["Monocyte"],dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=False)
cor_df = optim.cor_df
cor_df.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.csv')
pd.to_pickle(cor_df,'C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.pkl')


#%% 230209 add all combination
import itertools
cor_df = pd.read_pickle('C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.pkl')
add_df = pd.DataFrame([0.322,('B','Basophil','CD4','CD8','Cholangiocyte','Eosinophil','Hepatocyte','Kupffer',
                                            'LSEC','Monocyte','NK','Neutrophil','Stellate')],index=cor_df.columns.tolist()).T
update_df = pd.concat([cor_df,add_df])
pd.to_pickle(update_df,'C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.pkl')
update_df.to_csv('C:/github/LiverDeconv/_Figures_new/Figure_4/data/mon_fluctuate_quartile_comb2ref.csv')

