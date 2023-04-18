# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 17:54:40 2023

@author: I.Azuma
"""
import pandas as pd

import sys
sys.path.append('C:/github/LiverDeconv')
from _utils import plot4deconv as p4d

#%% all samples
facs_res = pd.read_csv('C:/github/LiverDeconv/Data/processed/facs_true_population.csv',index_col=0)
remove_list = ['CIV_1','CIV_4','CIV_7','CIV_8','CIP_1','CIP_2','CIP_7','CIP_8']
common_list = list(set(facs_res.columns.tolist()) & set(remove_list))
removed_df = facs_res.drop(columns=common_list) # remove

p4d.plot_violin(target_df=removed_df.T,sep="_",sort_index=["APAP","MDA","ANIT","TAA","ConA","GAL","CCl4"],
                target_cell=[],ctrl="Ctrl",
                x_doc="",y_doc="score",row_n=3,col_n=3)

#%% all samples
facs_res = pd.read_csv('C:/github/LiverDeconv/Data/processed/facs_true_population.csv',index_col=0)
remove_list = ['APAP_8','APAP_11','CIV_1','CIV_4','CIV_7','CIV_8','CIP_1','CIP_2','CIP_7','CIP_8']
common_list = list(set(facs_res.columns.tolist()) & set(remove_list))
removed_df = facs_res.drop(columns=common_list) # remove

p4d.plot_violin(target_df=removed_df.T,sep="_",sort_index=["APAP","MDA","ANIT","TAA","ConA","GAL","CCl4"],
                target_cell=[],ctrl="Ctrl",
                x_doc="",y_doc="score",row_n=3,col_n=3)