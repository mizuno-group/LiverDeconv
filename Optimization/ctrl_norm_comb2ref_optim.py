# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:13:55 2022

1. determine reference cell combination
2. define marker gene expression
3. evaluate correlation after normalizing with Ctrl sample

@author: I.Azuma
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

from _utils import processing
from Optimization import sample_selector as ss
import liver_deconv
from Evaluation import evaluator, plot4eval

target_candi = ["Neutrophil","Monocyte","NK","Eosinophil"]
method_candi = ["z-score","outerquartile"]
#%%
class Ctrlnorm_comb2ref():
    def __init__(self):
        self.mix_df = None # deconvolution target samples
        self.val_df = None # immune population validation data
        self.raw_ref = None # cell type specific reference data
        self.__processing = processing

#%% sample preparation functions
    def set_datasets(self,mix_df=None,val_df=None,ignore_samples=['CIV_1','CIV_4','CIV_7','CIV_8','CIP_1','CIP_2','CIP_7','CIP_8']):
        """
        mix_df: (genes, samples)
        val_df: (cells, samples)
        """
        self.common_samples = list(sorted((set(mix_df.columns.tolist()) & set(val_df.columns.tolist()))-set(ignore_samples)))
        
        self.mix_df = mix_df[self.common_samples]
        self.val_df = val_df[self.common_samples]
        
    def set_target(self,target="Neutrophil"):
        print("target: ",target)
        if target in target_candi:
            self.target_cell = target
        else:
            raise ValueError("!! Inappropriate target !!")
    
    def select_samples(self,method="z-score",ctrl="Ctrl",**kwargs):
        if method == "z-score":
            self.target_sample, self.norm_val = ss.zscore_sampling(val_df=self.val_df,
                                                                   ctrl=ctrl,
                                                                   target_cell=self.target_cell,
                                                                   upper_z=kwargs.get("upper_z"),
                                                                   lower_z=kwargs.get("lower_z"),
                                                                   outlier=kwargs.get("outlier"),
                                                                   centering=kwargs.get("centering"),
                                                                   ignore_ctrl=kwargs.get("ignore_ctrl"))
        elif method == "outerquartile":
            self.target_sample, self.norm_val = ss.quartile_sampling(val_df=self.val_df,
                                                                     ctrl=ctrl,
                                                                     target_cell=self.target_cell,
                                                                     ignore_ctrl=kwargs.get("ignore_ctrl"),
                                                                     threshold=kwargs.get("threshold"),
                                                                     q=kwargs.get("q"))
        else:
            raise ValueError("!! Inappropriate method !!")
        
    
    def select_samples_legacy(self,ctrl="C",upper_z=1.0,lower_z=-1.0,outlier=3.0,centering=False,ignore_ctrl=True):
        """
        -outlier < x < lower_z
        upper_z < x < outlier
        """
        # normalize the validation data (FACS) with ctrl values
        self.norm_val = self.__processing.ctrl_norm(self.val_df.T,ctrl=ctrl)
        if ignore_ctrl:
            use = []
            for t in self.norm_val.index.tolist():
                if t.split("_")[0] == ctrl:
                    pass
                else:
                    use.append(t)
            wo_norm_val = self.norm_val.loc[use]
            self.z_norm_val = self.__processing.standardz_sample(wo_norm_val) # remove ctrl samples
        else:
            self.z_norm_val = self.__processing.standardz_sample(self.norm_val)
            # centering with Ctrl value
            if centering:
                self.z_norm_val = self.__processing.ctrl_centering(self.z_norm_val)
        # plot the current situation
        plot4select(self.z_norm_val,upper_z=upper_z,lower_z=lower_z)
        # re add ctrl samples (use in deconvolution output normalization)
        self.z_norm_val = self.__processing.standardz_sample(self.norm_val)
        target_df = self.z_norm_val[[self.target_cell]]
        #self.target_sample = target_df[((target_df[self.target_cell]>upper_z)|(target_df[self.target_cell]<lower_z))].index.tolist()
        #self.target_sample = target_df.query('3.0 > Monocyte > 0.5 | -3.0 < Monocyte < -0.5')
        self.target_sample = target_df.query(str(outlier)+' > '+self.target_cell+' > '+str(upper_z)+' | '+str(-outlier)+' < '+self.target_cell+' < '+str(lower_z)).index.tolist()
        print("--- evaluation ---")
        print("target: ",self.target_cell)
        print(len(self.target_sample),"/",len(self.z_norm_val),"samples are analytical target")
        
        if len(self.target_sample) < 20:
            print("!! Attention! Sample size may be too small. !!")
        
        # reflect target sample
        self.target_mix = self.mix_df[self.target_sample]
        self.target_val = self.z_norm_val.T[self.target_sample]

#%% optimization functions
    def optim_loop(self,raw_ref=None,candi=['B','Basophil','CD4','CD8','Cholangiocyte','Eosinophil','Hepatocyte','Kupffer',
                                            'LSEC','Monocyte','NK','Neutrophil','Stellate'],
                   target=["Monocyte"],dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=False):
        self.raw_ref = raw_ref
        other_cells = sorted(list(set(candi)-set(target)))
        self.all_cor = []
        names = []
        for i in range(1,len(other_cells)+1):
            l = list(itertools.combinations(other_cells,i))
            for factor in l:
                print(factor)
                names.append(factor)
                # extract the target combination
                self.selected_ref = ref_comb_selector(ref=raw_ref,target_cell=target[0],target_comb=factor)
                
                dat = liver_deconv.LiverDeconv()
                dat.set_data(df_mix=self.mix_df,df_all=self.selected_ref)
                dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=False)
                dat.narrow_intersec()
                dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=False,do_plot=False)
                self.final_ref = dat.final_ref
                
                cor,self.final_dec,self.final_val = self.single_run(mix_df=self.mix_df, sig_df=self.final_ref, val_df=self.norm_val,dec_names=dec_names,val_names=val_names,do_plot=do_plot)
                """
                mix_df and val_df are both before sample narrowing
                """
                self.all_cor.append(cor)
        cor_df = pd.DataFrame(self.all_cor,columns=target,index=names)
        cor_df["pair"]=names
        self.cor_df = cor_df
    
    def optimized_single(self,raw_ref=None,combination=('B', 'CD4', 'NK'),target=["Monocyte"],dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=False,dpi=100):
        """
        run with the specific reference cell combination
        """
        self.raw_ref = raw_ref
        # extract the target combination
        self.selected_ref = ref_comb_selector(ref=raw_ref,target_cell=target[0],target_comb=combination)
        
        dat = liver_deconv.LiverDeconv()
        dat.set_data(df_mix=self.mix_df,df_all=self.selected_ref)
        dat.pre_processing(do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=False)
        dat.narrow_intersec()
        dat.create_ref(sep="_",number=50,limit_CV=10,limit_FC=1.5,log2=False,verbose=False,do_plot=False)
        self.final_ref = dat.final_ref
        
        cor,self.final_dec,self.final_val = self.single_run(mix_df=self.mix_df, sig_df=self.final_ref, val_df=self.norm_val,dec_names=dec_names,val_names=val_names,do_plot=do_plot,dpi=dpi)      

    def optim_base(self,base_sig_df=None,dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=True):
        # set data before sorting
        mix_df = self.mix_df
        norm_val = self.norm_val
        self.single_run(mix_df, base_sig_df, norm_val, dec_names=dec_names,val_names=val_names,do_plot=do_plot)
    
    def single_run(self,mix_df,sig_df,val_df,dec_names=[["Monocyte"]],val_names=[["Monocyte"]],do_plot=False,dpi=100):
        dec = liver_deconv.LiverDeconv()
        dec.do_fit(file_dat=mix_df,file_ref=sig_df)
        self.res = dec.get_res()
        self.norm_res = self.__processing.ctrl_norm(self.res) # normalize with ctrl sample results
        # reflect targe samples
        self.target_res = self.res.loc[self.target_sample]
        self.target_val = val_df.loc[self.target_sample]
        # remove ctrl samples
        remove_list=["Ctrl_1","Ctrl_10","Ctrl_12","Ctrl_15","Ctrl_16","Ctrl_17","Ctrl_18","Ctrl_2","Ctrl_3","Ctrl_4","Ctrl_7","Ctrl_8","Ctrl_9"]
        self.remove_res = remove_ctrl(self.target_res,remove_list,axis=1) # remove ctrl samples
        self.remove_val = remove_ctrl(self.target_val,remove_list,axis=1) # remove ctrl samples
        #z_res = self.__processing.standardz_sample(norm_res)
        #plot4deconv.plot_immune_box(z_res,sort_index = ["C","AP","MDA","AN","TA","CA","GAL","C4"],ctrl=["C"],row_n=3,col_n=5)
        if len(self.norm_res)==0:
            print("deconvolution running error")
            return np.nan
        else:
            ev = evaluator.Evaluator()
            ev.set_deconv_res(self.remove_res,z_score=True) # remove ctrl samples and calc z-score (samples, cell types)
            #ev.remove_samples(remove_list=["C_1","C_10","C_12","C_15","C_16","C_17","C_18","C_2","C_3","C_4","C_7","C_8","C_9"]) # remove ctrl samples
            ev.set_validation_ref(val_df=self.remove_val.T) # (cell types, samples)
            ev.process_validation_ref(z_score=True)
            ev.evaluate(dec_names=dec_names,val_names=val_names,sort_index=[],do_plot=do_plot,simple=False,eval_all=False,dpi=dpi)
            cor_res = ev.cor_res
            final_dec, final_val = ev.get_res()
            return list(cor_res.values()), final_dec, final_val
        
  

#%% global functions
def plot4select(df,upper_z=1.0,lower_z=-1.0,outlier=3):
    n = len(df.T)
    for i in range(n):
        tmp = df.iloc[:,i].tolist()
        plt.scatter([i]*len(tmp),tmp,s=1)
    plt.xticks([i for i in range(n)],df.columns.tolist(),rotation=45)
    plt.hlines(upper_z,-0.1,n-0.9,ls="dashed",color="black")
    plt.hlines(lower_z,-0.1,n-0.9,ls="dashed",color="black")
    if df.max().max() > outlier:
        plt.hlines(outlier,-0.1,n-0.9,ls="-.",color="red")
    if df.min().min() < -outlier:
        plt.hlines(-outlier,-0.1,n-0.9,ls="-.",color="red")
    plt.show()
        
def ref_comb_selector(ref,target_cell="Neutrophil",target_comb=('NK', 'CD8', 'Stellate', 'Kupffer')):
    use_samples = []
    for t in ref.columns.tolist():
        if t.split("_")[0] in target_comb:
            use_samples.append(t)
        elif t.split("_")[0] == target_cell:
            use_samples.append(t)
        else:
            pass
    selected_ref = ref[use_samples]
    return selected_ref

def remove_ctrl(df,remove_list=["C_1","C_10","C_12","C_15","C_16","C_17","C_18","C_2","C_3","C_4","C_7","C_8","C_9"],axis=0):
    """
    process validation reference data
    """
    if axis == 0:
        common_sample = list(set(remove_list) & set(df.columns.tolist()))
        target_df = df.drop(columns=common_sample)
    elif axis == 1:
        common_sample = list(set(remove_list) & set(df.index.tolist()))
        target_df = df.drop(index=common_sample)
    return target_df