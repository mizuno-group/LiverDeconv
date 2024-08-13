# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 19:48:32 2022

@author: I.Azuma
"""
import copy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet

from _utils import processing,plot4deconv
import fit


class LiverDeconv():
    def __init__(self):
        self.df_mix=pd.DataFrame()
        self.df_all=pd.DataFrame()
        self.__processing=processing
        self.method_dict={'ElasticNet':fit.fit_ElasticNet,'NuSVR':fit.fit_NuSVR,'NNLS':fit.fit_NNLS,'OLS':fit.fit_OLS,
                          'Ridge':fit.fit_Ridge,'Lasso':fit.fit_Lasso,'LR':fit.fit_LR,'RLR':fit.fit_RLR}
    
    def set_data(self,df_mix,df_all,verbose=True):
        """
        set input data
        ----------
        df_mix : DataFrame
            deconvolution target data (for which we want to know the population of each immune cell). 
            Set the genes and samples in rows and columns, respectively.
            e.g. /data/input/mix_processed.csv
                                AN_1       AN_2       AN_3  ...     GAL_4     GAL_7     GAL_8
            0610005C13Rik  11.875867  12.148424  11.252824  ...  8.431683  7.700966  8.962926
            0610009B22Rik   8.734242   8.457720   8.184438  ...  8.663929  8.546301  8.989601
            0610010F05Rik   6.490900   6.421802   5.806085  ...  6.207227  6.382519  6.357380
            0610010K14Rik   8.564556   8.673014   8.348932  ...  8.512630  8.407845  8.454941
            0610012G03Rik   8.097187   7.933714   7.698964  ...  8.108420  7.913200  8.247160
            
        df_all : DataFrame
            Reference data for estimate the population of the immune cells
            e.g. /data/input/ref_13types.csv
                           NK_GSE114827_1  NK_GSE114827_2  ...  NK_GSE103901_1  NK_GSE103901_2                                        ...                                
            0610005C13Rik       11.216826       23.474148  ...       15.364877       31.139537
            0610006L08Rik        0.000000        0.000000  ...        0.000000        0.000000
            0610009B22Rik     3801.734228     3838.271735  ...     1673.625997     1265.782913
            0610009E02Rik       38.492335       26.822363  ...       31.225440        5.380941
            0610009L18Rik      260.189214      144.910984  ...      639.627003      664.538873
            
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
        Note that mix_data is already processed (trimmed) in general (log2 --> trim+impute --> batch norm --> QN).
        This is because the robustness of the analysis is reduced if the number of genes to be analyzed is not narrowed down to a certain extent.
        """   
        # trimming
        mix_data = copy.deepcopy(self.df_mix)
        reference_data = copy.deepcopy(self.df_all)
        
        self.df_mix, self.df_all = self.__intersection_index(mix_data,reference_data) # update
        print("narrowd gene number :",len(self.df_all))
    
    def __intersection_index(self,df,df2):
            ind1 = df.index.tolist()
            ind2 = df2.index.tolist()
            df.index = [i.upper() for i in ind1]
            df2.index = [i.upper() for i in ind2]
            ind = list(set(df.index) & set(df2.index))
            df = df.loc[ind,:]
            df2 = df2.loc[ind,:]
            return df,df2
#%% 1. detect differentially expressed genes and create cell specific marker gene matrix (reference)
    # signature definition
    def create_ref(self,do_plot=False,**kwargs):
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
        
        signature = self.get_genes() # union of each reference cell's signatures
        print("signature genes :",len(signature))
        self.sig_ref = ref_inter_df.loc[signature]
        if sep == "":
            final_ref = self.sig_ref
        else:
            final_ref = self.__df_median(self.sig_ref,sep=sep)
        if do_plot:
            sns.clustermap(final_ref,col_cluster=False,z_score=0,figsize=(6, 6))
            plt.show()
        # FIXME:sort gene name
        s_gene = sorted(final_ref.index.tolist())
        final_ref = final_ref.loc[s_gene]
        self.final_ref = final_ref
        
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
        if sep_ind=="":
            immunes = df_c.columns.tolist()
        else:
            immunes = [t.split(sep_ind)[0] for t in df_c.columns.tolist()]
        df_c.columns = immunes
        self.min_FC = pd.DataFrame()
        self.pickup_genes_list = []
        self.__pickup_genes = []
        for c in sorted(list(set(immunes))):
            if sep_ind=="":
                self.df_target = pd.DataFrame(df_c[[c]])  # FIXME: 240813 
            else:
                self.df_target = pd.DataFrame(df_c[c]) # FIXME: 230117 
            self.tmp_summary = pd.DataFrame()
            for o in sorted(list(set(immunes))):
                if o == c:
                    pass
                else:
                    if sep_ind =="":
                        self.df_else = df_c[[o]]
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
                if verbose:
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
        
    def get_genes(self):
        self.__pickup_genes=[i for i in self.__pickup_genes if str(i)!='nan']
        self.__pickup_genes=list(set(self.__pickup_genes))
        return self.__pickup_genes
    
#%% 2. conduct deconvolution
    def do_fit(self,method="ElasticNet",file_dat=None,file_ref=None,prints=True,number_of_repeats=3,
               alpha=1,l1_ratio=0.05,max_iter=100000,nu=0.5,C=1,epsilon=1.35):
        
        # data input
        if file_dat is None:
            self.__set_mix_data(file_dat=self.df_mix)
        else:
            self.__set_mix_data(file_dat=file_dat)
        if file_ref is None:
            self.__set_reference_data(file_ref=self.final_ref)
        else:
            self.__set_reference_data(file_ref=file_ref)
        
        # data processing
        self.__mix_data=self.__calc_median_same_gene(self.__mix_data)
        self.__reference_data=self.__calc_median_same_gene(self.__reference_data)
        
        self.__mix_data = self.__processing.drop_all_missing(self.__mix_data)
        self.__reference_data = self.__processing.drop_all_missing(self.__reference_data)
        
        self.__gene_intersection(prints=prints)
        
        # fitting
        self.__fit(method=method,number_of_repeats=number_of_repeats,alpha=alpha,l1_ratio=l1_ratio,max_iter=max_iter,prints=prints,nu=nu,C=C,epsilon=epsilon)
        
        return
    
    def __fit(self,method='ElasticNet',number_of_repeats=1,alpha=1,l1_ratio=0.05,nu=0.5,C=1,max_iter=100000,epsilon=1.35,prints=True):
        if method in self.method_dict.keys():
            for i in range(number_of_repeats):
                res_mat = self.method_dict[method](self.__reference_data,self.__mix_data,alpha=alpha,l1_ratio=l1_ratio,nu=nu,max_iter=max_iter,epsilon=epsilon)
                # sum up
                if i == 0:
                    res = res_mat
                else:
                    res = res + res_mat
            if prints:
                print('fitting method : {}'.format(method))
            res = res / number_of_repeats
        else:
            print(self.method_dict.keys())
            raise ValueError("set appropriate method")
        res = res.astype(float)
        self.__res=res
        return
    
    def __fit_legacy(self,number_of_repeats=1,alpha=1,l1_ratio=0.05,max_iter=100000,prints=True):
        try:
            for i in range(number_of_repeats):
                res_mat = self.fit_ElasticNet(self.__reference_data,self.__mix_data,alpha=alpha,l1_ratio=l1_ratio,max_iter=max_iter)
                # sum up
                if i == 0:
                    res = res_mat
                else:
                    res = res + res_mat
            if prints:
                print('fitting method : ElasticNet')
            res = res / number_of_repeats
        except:
            res=np.nan
            print('fitting error')
            print('confirm fitting method name')
        self.__res=res
        return
    
    def __calc_median_same_gene(self,df):
            df = df.dropna()
            df2 = pd.DataFrame()
            dup = df.index[df.index.duplicated(keep="first")]
            gene_list = pd.Series(dup).unique().tolist()
            if len(gene_list) != 0:
                for gene in gene_list:
                    new = df.loc[gene].median()
                    df2[gene] = new
                df = df.drop(gene_list)
                df = pd.concat([df,df2.T])
            return df
    
    def __gene_intersection(self,prints=True):
            # delete nan / inf
            ref = self.__reference_data.replace(np.inf,np.nan)
            ref = ref.replace(-np.inf,np.nan)
            ref = ref.dropna()
            
            dat = self.__mix_data.replace(np.inf,np.nan)
            dat = dat.replace(-np.inf,np.nan)
            dat = dat.dropna()
            
            # upper gene name
            ref.index = [i.upper() for i in list(ref.index)]
            dat.index = [i.upper() for i in list(dat.index)]
            
            # intersection
            marker = set(ref.index)
            genes = set(dat.index)
            marker = list(marker&genes)
            self.__reference_data = ref.loc[marker,:]
            self.__mix_data = dat.loc[marker,:]
            if prints:
                print("number of used genes = {}".format(len(marker)))
            return
    
    ### in/out put control ###
    def __set_mix_data(self,file_dat):
        try:
            self.__mix_data = pd.read_csv(file_dat,index_col=0)
        except:
            self.__mix_data = file_dat

            
    def __set_reference_data(self,file_ref):
        try:
            self.__reference_data = pd.read_csv(file_ref,index_col=0)
        except:
            self.__reference_data = file_ref
    
    def get_res(self):
        return self.__res
    
    def get_data(self):
        return self.__mix_data, self.__reference_data

#%% 3. summarize the output
    def summarize(self,remove_list=["APPAP_8","APAP_11","CIV_7","CIV_8","CIP_7","CIP_8"],
                  sort_index=["APAP","MDA","ANIT","TAA","ConA","GAL","CCl4"],
                  x_doc="",y_doc="score",ctrl="Ctrl",row_n=0,col_n=0,z_score=True):
        """plot the estimated cell populations"""
        res = self.get_res()
        df = copy.deepcopy(res)
        common_list = list(set(df.index.tolist()) & set(remove_list))
        removed_df = df.drop(index=common_list) # remove
        if z_score:
            z_res = self.__processing.standardz_sample(removed_df)
            plot4deconv.plot_violin(z_res,sep="_",sort_index=sort_index,target_cell=[],ctrl=ctrl,x_doc=x_doc,y_doc=y_doc,row_n=row_n,col_n=col_n)
        else:
            plot4deconv.plot_violin(removed_df,sep="_",sort_index=sort_index,target_cell=[],ctrl=ctrl,x_doc=x_doc,y_doc=y_doc,row_n=row_n,col_n=col_n)
