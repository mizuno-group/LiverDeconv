# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 11:58:16 2022

@author: I.Azuma
"""
import copy
import pandas as pd
from collections import defaultdict

from _utils import processing
from Evaluation import plot4eval as p4v

class Evaluator():
    def __init__(self):
        self.deconv_res = None
        self.__processing = processing
    
    def set_deconv_res(self,res_df,z_score=False):
        """
        set the deconvolution result
        
                     B       CD4       CD8  Monocytes        NK  Neutrophils
        AN_1  0.059349 -0.032014  0.074860   0.093262 -0.073984     0.056436
        AN_2  0.064035 -0.034487  0.089953   0.075578 -0.064471     0.028734
        AN_3  0.052030 -0.003279  0.083852   0.076546 -0.087147     0.052131
        AN_4  0.060684 -0.035129  0.089692   0.072977 -0.073421     0.043747
        AP_1  0.056954 -0.012693  0.116242   0.054742 -0.098874     0.073787
        
        """
        if z_score:
            self.deconv_res = self.__processing.standardz_sample(res_df)
            if len(self.deconv_res)<1:
                raise ValueError('!! z-score normalization is not appropriate !!')
        else:
            self.deconv_res = res_df
    
    def remove_samples(self,remove_list=["APAP_8","APAP_11","CIV_7","CIV_8","CIP_7","CIP_8"]):
        """
        remove samples not included in the analysis target
        e.g. Samples with extremely high ALT values or  Ctrl samples with inflammation.
        """
        df = copy.deepcopy(self.deconv_res)
        common_list = list(set(df.index.tolist()) & set(remove_list))
        removed_df = df.drop(index=common_list) # remove
        print("original :",len(df))
        print("after removing :",len(removed_df))
        self.deconv_res = removed_df
        
    def set_validation_ref(self,val_df=None):
        """Ground Truth data (e.g. FACS measured value)

        Args:
            val_df (DataFramr): Cell types are in rows and sample names are in columns. Defaults to None.
                        Ctrl_1    Ctrl_2    Ctrl_3    Ctrl_4    Ctrl_7    Ctrl_8  
            abT         0.144218  0.354034  0.128803  0.324972  0.220206  0.237484   
            gdT         0.015701  0.044586  0.017606  0.031974  0.077327  0.081104   
            NKT         0.207992  0.356555  0.173388  0.291621  0.057071  0.056809   
            NK          0.073829  0.156051  0.084659  0.148015  0.069933  0.064143   
            Monocyte    0.042073  0.037323  0.033776  0.038917  0.070734  0.088048  
        """
        if val_df is None:
            self.val_df = pd.read_csv('C:/github/LiverDeconv/Data/processed/facs_true_population.csv',index_col=0)
        else:
            self.val_df = val_df
    
    def process_validation_ref(self,z_score=True):
        """
        process validation reference data
        """
        val_df = copy.deepcopy(self.val_df)
        val_df.columns = [str(t) for t in val_df.columns.tolist()] # str
        
        common_sample = list(set(val_df.columns.tolist()) & set(self.deconv_res.index.tolist()))
        target_res = self.deconv_res.loc[common_sample]
        self.deconv_res = target_res
        target_val = val_df[common_sample].T # select samples in deconv_res
        if z_score:
            self.target_val_ref = self.__processing.standardz_sample(target_val)
        else:
            self.target_val_ref = target_val
    
    def evaluate(self,dec_names=[["Neutrophil"],["Monocyte"],["NK"],["CD4","CD8"]],val_names=[["Neutrophil"],["Monocyte"],["NK"],["abT"]],sort_index=["Ctrl","APAP","MDA","ANIT","TAA","ConA","GAL","CCl4"],title=None,do_plot=True,simple=False,eval_all=False,dpi=300):
        print("----------")
        self.cor_res = defaultdict(float)
        self.incorrect_res = defaultdict(list)

        # select common sample
        common = sorted(list(set(self.deconv_res.index.tolist()) & set(self.target_val_ref.index.tolist())))
        deconv_res = self.deconv_res.loc[common]
        target_val_ref = self.target_val_ref.loc[common]

        if simple:
            X1 = []
            X2 = []
            labels = []
            for i in range(len(dec_names)):
                cor,x1,x2,label = p4v.plot_simple_corr(deconv_res,target_val_ref,dec_name=dec_names[i],val_name=val_names[i],do_plot=do_plot,dpi=dpi)
                X1.append(x1)
                X2.append(x2)
                labels.append(label)
                self.cor_res[val_names[i][0]] += cor
            self.total_cor = p4v.plot_multi_corr(X1=X1,X2=X2,labels=labels,dpi=dpi)
        else:
            for i in range(len(dec_names)):
                x,y,cor = p4v.plot_value_corr(deconv_res,target_val_ref,dec_name=dec_names[i],val_name=val_names[i],
                                              sort_index=sort_index,do_plot=do_plot,dpi=dpi,title=title)
                #print(x,y)
                self.cor_res[val_names[i][0]] += cor
                if eval_all:
                    incorrect = p4v.plot_mean_change(deconv_res,target_val_ref,dec_name=dec_names[i],val_name=val_names[i],do_plot=do_plot,dpi=dpi)
                    self.incorrect_res[val_names[i][0]].append(incorrect)
                    if do_plot:
                        p4v.plot_box_change(deconv_res,target_val_ref,dec_name=dec_names[i],val_name=val_names[i],sort_index = sort_index,dpi=dpi)
                else:
                    pass
    
    def get_res(self):
        dec_idx = self.deconv_res.index.tolist()
        self.target_val_ref = self.target_val_ref.loc[dec_idx]
        return self.deconv_res, self.target_val_ref
            