# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 11:16:33 2022

@author: I.Azuma
"""
import pandas as pd
import copy
from scipy import stats as st
import statsmodels.stats.multitest as sm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import random

from _utils import processing

class Deg_Multi_FC():
    def __init__(self):
        self.df_mix=pd.DataFrame()
        self.df_all=pd.DataFrame()
        self.__processing=processing
        
    ### main ###
    def set_data(self,df_mix,df_all,verbose=True):
        """
        set directly
        """
        self.df_mix = df_mix
        self.df_all = df_all
        
        if verbose:
            print(self.df_mix.shape)
            print(self.df_all.shape)
    
    def pre_processing(self,do_ann=False,ann_df=None,do_log2=True,do_quantile=True,do_trimming=False,do_drop=True):
        if do_ann:
            self.df_all = self.__processing.annotation(self.df_all,ann_df)
        else:
            pass
        if do_log2:
            df_c = copy.deepcopy(self.df_all)
            self.df_all = self.__processing.log2(df_c)
            print("log2 conversion")
        else:
            pass
        if do_quantile:
            df_c = copy.deepcopy(self.df_all)
            self.df_all = self.__processing.quantile(df_c)
            print("quantile normalization")
        else:
            pass
        if do_trimming:
            df_c = copy.deepcopy(self.df_all)
            raw_batch = [t.split("_")[0] for t in df_c.columns.tolist()]
            batch_list = pd.Series(raw_batch).astype("category").cat.codes.tolist()
            self.df_all = self.__processing.array_imputer(df_c,threshold=0.9,strategy="median",trim=1.0,batch=True,lst_batch=batch_list, trim_red=False)
            print("trimming")
        else:
            pass
        if do_drop:
            df_c = copy.deepcopy(self.df_all)
            replace = df_c.replace(0,np.nan)
            drop = replace.dropna(how="all")
            drop = drop.fillna(0)
            self.df_all = drop
            print("drop nan")
    
    def narrow_intersec(self):
        """
        df_mix : deconvolution target data (for which we want to know the population of each immune cell)
        df_all : reference data for estimate the population of the immune cells (e.g. LM6)
        
        Note that mix_data is already processed (trimmed) in general (log2 --> trim+impute --> batch norm --> QN).
        This is because the robustness of the analysis is reduced if the number of genes to be analyzed is not narrowed down to a certain extent.
        """   
        # trimming
        mix_data = copy.deepcopy(self.df_mix)
        reference_data = copy.deepcopy(self.df_all)
        
        self.df_mix, self.df_all = self.__intersection_index(mix_data,reference_data) # update
        print("narrowd gene number :",len(self.df_all))
    
    def deg_extraction(self,sep_ind="_",number=150,limit_CV=0.1,limit_FC=1.5,verbose=True):
        """
        Define DEGs between the target and other one for the cells that make up the REFERENCE.
        e.g. B cell vs CD4, B cell vs CD8, ...

        Parameters
        ----------
        sep_ind : str
            Assume the situation that the columns name is like "CellName_GSEXXX_n". The default is "_".
        number : int
            Number of top genes considered as DEGs. The default is 150.
        limit_CV : float
            Coefficient of Variation threshold. The default is 0.3.
        limit_FC : TYPE, float
            Minimum threshold for logFC detection. The default is 1.5.

        """
        df_c = copy.deepcopy(self.df_all)
        #cluster, self.samples = self.sepmaker(df=df_c,delimiter=sep_ind)
        #print(cluster)
        if sep_ind == "":
            immunes = df_c.columns.tolist()
        else:
            immunes = [t.split(sep_ind)[0] for t in df_c.columns.tolist()]
        df_c.columns = immunes
        self.min_FC = pd.DataFrame()
        self.pickup_genes_list = []
        self.__pickup_genes = []
        for c in sorted(list(set(immunes))):
            self.df_target = df_c[c]
            self.tmp_summary = pd.DataFrame()
            for o in sorted(list(set(immunes))):
                if o == c:
                    pass
                else:
                    self.df_else = df_c[o]
                    self.__logFC()
                    df_logFC = self.df_logFC
                    df_logFC.columns = [o]
                    self.tmp_summary = pd.concat([self.tmp_summary,df_logFC],axis=1)
            tmp_min = self.tmp_summary.T.min()
            self.df_minFC = pd.DataFrame(tmp_min)
            self.__calc_CV()
            
            pickup_genes = self.__selection(number=number,limit_CV=limit_CV,limit_FC=limit_FC,verbose=verbose)
            self.pickup_genes_list.append(pickup_genes)
            self.min_FC = pd.concat([self.min_FC,tmp_min],axis=1)
        self.min_FC.columns = sorted(list(set(immunes)))
        self.pickup_genes_df=pd.DataFrame(self.pickup_genes_list).T.dropna(how="all")
        self.pickup_genes_df.columns = sorted(list(set(immunes)))
        curate = [[i for i in t if str(i)!='nan'] for t in self.pickup_genes_list]
        self.deg_dic = dict(zip(sorted(list(set(immunes))),curate))
    
    def create_ref(self,**kwargs):
        """
        create reference dataframe which contains signatures for each cell

        """
        ref_inter_df = copy.deepcopy(self.df_all)
        df2 = copy.deepcopy(self.df_all)
        if kwargs["log2"]:
            df2 = np.log2(df2+1)
        else:
            pass
        self.df_all = df2
        
        # DEG extraction
        sep = kwargs["sep"]
        self.deg_extraction(sep_ind=sep,number=kwargs["number"],limit_CV=kwargs["limit_CV"],limit_FC=kwargs["limit_FC"],verbose=kwargs["verbose"])
        
        signature = self.get_res() # union of each reference cell's signatures
        print("signature genes :",len(signature))
        self.sig_ref = ref_inter_df.loc[signature]
        final_ref = self.__df_median(self.sig_ref,sep=sep)
        #sns.clustermap(final_ref,col_cluster=False,z_score=0)
        #plt.show()
        self.final_ref = final_ref
    
    def create_ref_legacy(self,sep="_",number=200,limit_CV=1,limit_FC=1.5,log2=False):
        """
        create reference dataframe which contains signatures for each cell

        """
        ref_inter_df = copy.deepcopy(self.df_all)
        df2 = copy.deepcopy(self.df_all)
        if log2:
            df2 = np.log2(df2+1)
        else:
            pass
        self.df_all = df2
        
        # DEG extraction
        self.deg_extraction(sep_ind=sep,number=number,limit_CV=limit_CV,limit_FC=limit_FC)
        
        signature = self.get_res() # union of each reference cell's signatures
        print("signature genes :",len(signature))
        sig_ref = ref_inter_df.loc[signature]
        final_ref = self.__df_median(sig_ref,sep=sep)
        sns.clustermap(final_ref,col_cluster=False,z_score=0)
        plt.show()
        self.final_ref = final_ref
    
    def create_random_ref(self,sep="_",seed=123,high_cut=100.0,low_cut=0.0,do_plot=False):
        """
        create randomized reference which is the same size to the correct reference used in deconvolution
        """
        if self.final_ref is None:
            raise ValueError("!! Conduct create_ref() at fist !!")
        n = len(self.final_ref)
        
        tmp_df = copy.deepcopy(self.df_all)
        s_max = pd.DataFrame(tmp_df.T.max()).sort_values(0)
        selected = s_max[(s_max[0]>low_cut)&(s_max[0]<high_cut)]
        if len(selected)<n:
            print(len(selected),"<",n)
            raise ValueError("!! Threshold setting is too strict !!")
        
        random.seed(seed)
        random_sig = random.sample(selected.index.tolist(),n)
        ref_inter_df = copy.deepcopy(self.df_all)
        random_tmp = ref_inter_df.loc[random_sig]
        random_ref = self.__df_median(random_tmp,sep=sep)
        if do_plot:
            try:
                sns.clustermap(random_ref,col_cluster=False,z_score=0)
            except:
                sns.clustermap(random_ref,col_cluster=False)
                print("without Z-score due to the inf values error")
            plt.show()
        else:
            pass
        self.random_ref = random_ref
        
    def get_res(self):
        self.__pickup_genes=[i for i in self.__pickup_genes if str(i)!='nan']
        self.__pickup_genes=list(set(self.__pickup_genes))
        return self.__pickup_genes
    
    def deg_exterior(self):
        """
        Overlook the DEGs definition condition with heatmap plotting
        """
        df = copy.deepcopy(self.df_all)
        df.index = [t.upper() for t in df.index.tolist()]
        df.fillna(0)
        for i,sample in enumerate(self.samples):
            tmp_df = df.loc[[t for t in self.pickup_genes_df[sample] if str(t)!='nan']]
            sns.heatmap(tmp_df)
            plt.title(sample)
            plt.show()
    
    ### processing ###
    def __intersection_index(self,df,df2):
        ind1 = df.index.tolist()
        ind2 = df2.index.tolist()
        df.index = [i.upper() for i in ind1]
        df2.index = [i.upper() for i in ind2]
        ind = list(set(df.index) & set(df2.index))
        df = df.loc[ind,:]
        df2 = df2.loc[ind,:]
        return df,df2
    
    def sepmaker(self,df=None,delimiter='.'):
        if df is None:
            df = self.df_all
        else:
            pass
        samples = list(df.columns)
        sample_unique = []
        seps=[]
        ap1 = sample_unique.append
        ap2 = seps.append
        for i in samples:
            if i.split(delimiter)[0] in sample_unique:
                number = sample_unique.index(i.split(delimiter)[0])
                ap2(number)
            else:
                ap1(i.split(delimiter)[0])
                ap2(len(sample_unique)-1)
        return seps, sample_unique
    
    def __df_median(self,df,sep="_"):
        df_c = copy.deepcopy(df)
        df_c.columns=[i.split(sep)[0] for i in list(df_c.columns)]
        df_c = df_c.groupby(level=0,axis=1).median()
        return df_c
    
    ### calculation ###
    def __logFC(self):
        # calculate df_target / df_else logFC
        df_logFC = self.df_target.T.median() - self.df_else.T.median()
        df_logFC = pd.DataFrame(df_logFC)
        self.df_logFC = df_logFC
    
    def __calc_CV(self):
        """
        CV : coefficient of variation
        """
        df_CV = np.std(self.df_target.T) / np.mean(self.df_target.T)
        df_CV.index = self.df_target.index
        df_CV = df_CV.replace(np.inf,np.nan)
        df_CV = df_CV.replace(-np.inf,np.nan)
        df_CV = df_CV.dropna()
        self.df_CV=pd.DataFrame(df_CV)
    
    def __selection(self,number=50,limit_CV=0.1,limit_FC=1.5,verbose=False):
        self.__intersection()
        df_minFC=self.df_minFC
        df_CV=self.df_CV
        df_minFC=df_minFC.sort_values(0,ascending=False)
        genes=df_minFC.index.tolist()
    
        pickup_genes=[]
        ap = pickup_genes.append
        i=0
        while len(pickup_genes)<number:
            if len(genes)<i+1:
                pickup_genes = pickup_genes+[np.nan]*number
                if  verbose:
                    print('not enough genes picked up')
            elif df_CV.iloc[i,0] < limit_CV and df_minFC.iloc[i,0] > limit_FC:
                ap(genes[i])
            i+=1
        else:
            self.__pickup_genes = self.__pickup_genes + pickup_genes
            return pickup_genes
    
    def __intersection(self):
        lis1 = list(self.df_minFC.index)
        lis2 = list(self.df_CV.index)
        self.df_minFC = self.df_minFC.loc[list(set(lis1)&set(lis2))]
        self.df_CV = self.df_CV.loc[list(set(lis1)&set(lis2))]
    
