# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 10:58:22 2022

@author: I.Azuma
"""
import copy
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import chain
import matplotlib.pyplot as plt
from matplotlib import patches

dpi = 300

def plot_value_corr(deconv_df,val_df,dec_name=["CD4","CD8"],val_name=["abT"],sort_index=[],do_plot=True,sep=True,title=None,do_print=False,dpi=300):
    """
    Correlation Scatter Plotting
    Format of both input dataframe is as follows
    
                 B       CD4       CD8  Monocytes        NK  Neutrophils
    AN_1 -0.327957 -0.808524 -0.768420   0.311360  0.028878     0.133660
    AN_2  0.038451 -0.880116 -0.278970  -1.039572  0.865344    -0.437588
    AN_3 -0.650633  0.574758 -0.498567  -0.796406 -0.100941     0.035709
    AN_4 -0.479019 -0.005198 -0.675028  -0.787741  0.343481    -0.062349
    AP_1 -1.107050  0.574758  0.858366  -1.503722 -1.053643     1.010999
    
    """
    if title is None:
        title = str(dec_name)+" vs "+str(val_name)
    
    if len(sort_index)>0:
        drugs = sort_index
    elif sep:
        drugs = sorted(list(set([t.split("_")[0] for t in deconv_df.index.tolist()])))
    else:
        drugs = sorted(deconv_df.index.tolist())
    
    # align the index
    val_df = val_df.loc[deconv_df.index.tolist()]
    
    dec_min = 100
    dec_max = 0
    total_x = deconv_df[dec_name].sum(axis=1).tolist()
    total_y = val_df[val_name].sum(axis=1).tolist()
    total_cor = round(np.corrcoef(total_x,total_y)[0][1],4)
    
    if do_print:
        print(str(dec_name)+" vs "+str(val_name))
        print(total_cor)
    
    if do_plot:
        fig,ax = plt.subplots(figsize=(6,6),dpi=dpi)
        for i,d in enumerate(drugs):
            tmp1 = deconv_df.filter(regex="^"+d+"_",axis=0)
            tmp2 = val_df.filter(regex="^"+d+"_",axis=0)
            
            res1 = tmp1[dec_name].sum(axis=1).tolist()
            res2 = tmp2[val_name].sum(axis=1).tolist()
            tmp_cor = round(np.corrcoef(res1,res2)[0][1],3)
        
            plt.scatter(res1,res2,label=d+" : "+str(tmp_cor),alpha=1.0)
            
            if min(res1)<dec_min:
                dec_min = min(res1)
            if max(res1)>dec_max:
                dec_max = max(res1)
        
        plt.plot([dec_min,dec_max],[dec_min,dec_max],linewidth=2,color='black',linestyle='dashed',zorder=-1)
        
        plt.text(0.3,0.05,'Cor = {}'.format(str(round(total_cor,3))), transform=ax.transAxes, fontsize=15)
        
        plt.legend(shadow=True)
        plt.xlabel("Deconvolution estimated value",fontsize=15)
        plt.ylabel("Flow cytometry value",fontsize=15)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        ax.set_axisbelow(True)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.title(title,fontsize=15)
        plt.show()
    else:
        pass
    return total_x,total_y,total_cor

def plot_simple_corr(deconv_df,val_df,dec_name=["CD4"],val_name=["CD4T"],do_plot=True,sep=True,do_print=False,dpi=300):
    """
    Correlation Scatter Plotting
    Format of both input dataframe is as follows
    Note that the targe data contains single treatment group (e.g. APAP treatment only)
    
                 B       CD4       CD8      Monocytes        NK  Neutrophils
    Donor_1 -0.327957 -0.808524 -0.768420   0.311360  0.028878     0.133660
    Donor_2  0.038451 -0.880116 -0.278970  -1.039572  0.865344    -0.437588
    Donor_3 -0.650633  0.574758 -0.498567  -0.796406 -0.100941     0.035709
    Donor_4 -0.479019 -0.005198 -0.675028  -0.787741  0.343481    -0.062349
    
    """
    dec_min = 100
    dec_max = 0
    total_cor = round(np.corrcoef(deconv_df[dec_name].sum(axis=1).tolist(),val_df[val_name].sum(axis=1).tolist())[0][1],4)
    
    if do_print:
        print(str(dec_name)+" vs "+str(val_name))
        print(total_cor)
    
    res1 = deconv_df[dec_name].sum(axis=1).tolist()
    res2 = val_df[val_name].sum(axis=1).tolist()
    tmp_cor = round(np.corrcoef(res1,res2)[0][1],3)
    label = dec_name[0]+" : "+str(tmp_cor)
    
    if min(res1)<dec_min:
        dec_min = min(res1)
    if max(res1)>dec_max:
        dec_max = max(res1)
    
    if do_plot:
        fig,ax = plt.subplots(figsize=(6,6),dpi=dpi)
        plt.scatter(res1,res2,label=label,alpha=1.0)
        plt.plot([dec_min,dec_max],[dec_min,dec_max],linewidth=2,color='black',linestyle='dashed',zorder=-1)
        plt.text(0.3,0.05,'Cor = {}'.format(str(round(total_cor,3))), transform=ax.transAxes, fontsize=15)
        
        plt.legend(shadow=True)
        plt.xlabel("Deconvolution estimated value")
        plt.ylabel("FACS value")
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        ax.set_axisbelow(True)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.title(str(dec_name)+" vs "+str(val_name))
        plt.show()
    else:
        pass
    return total_cor,res1,res2,label

def plot_multi_corr(X1=[[1,2,3],[2,3,4]],X2=[[2,4,5],[3,5,7]],labels=["B","NK"],dpi=300):
    """
    Plot comprehensive results obtained from plot_simple_corr()
    """
    dec_min = 100
    dec_max = 0
    fig,ax = plt.subplots(figsize=(6,6),dpi=dpi)
    for i in range(len(X1)):
        x1 = X1[i]
        x2 = X2[i]
        label = labels[i]
        plt.scatter(x1,x2,label=label,alpha=1.0)
        
        if min(x1)<dec_min:
            dec_min = min(x1)
        if max(x1)>dec_max:
            dec_max = max(x1)
    
        plt.legend(shadow=True)
        plt.xlabel("Deconvolution estimated value")
        plt.ylabel("FACS value")
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        ax.set_axisbelow(True)
        ax.grid(color="#ababab",linewidth=0.5)
    
    x1chain = list(chain.from_iterable(X1))
    x2chain = list(chain.from_iterable(X2))
    total_cor = round(np.corrcoef(x1chain,x2chain)[0][1],5)
    plt.text(0.3,0.05,'Total Corr = {}'.format(str(round(total_cor,3))), transform=ax.transAxes, fontsize=15) # show the total Pearson Correlation
    plt.plot([dec_min,dec_max],[dec_min,dec_max],linewidth=2,color='black',linestyle='dashed',zorder=-1) # plot diagonal
    return total_cor
        

def plot_mean_change(deconv_df,val_df,dec_name=["CD4","CD8"],val_name=["abT"],do_plot=True,dpi=300):
    """
    Correlation Scatter Plotting of each change mean value to Ctrl samples
    Format of both input dataframe is as follows
    
                 B       CD4       CD8  Monocytes        NK  Neutrophils
    AN_1 -0.327957 -0.808524 -0.768420   0.311360  0.028878     0.133660
    AN_2  0.038451 -0.880116 -0.278970  -1.039572  0.865344    -0.437588
    AN_3 -0.650633  0.574758 -0.498567  -0.796406 -0.100941     0.035709
    AN_4 -0.479019 -0.005198 -0.675028  -0.787741  0.343481    -0.062349
    AP_1 -1.107050  0.574758  0.858366  -1.503722 -1.053643     1.010999
    
    """
    df = copy.deepcopy(deconv_df)
    df["drug"] = [t.split("_")[0] for t in df.index.tolist()]
    mean_df = df.groupby("drug").mean()
    c = mean_df.loc["C"]
    dec_mean_change = mean_df - c
    
    df2 = copy.deepcopy(val_df)
    df2["drug"] = [t.split("_")[0] for t in df2.index.tolist()]
    dec_mean_df = df2.groupby("drug").mean()
    c2 = dec_mean_df.loc["C"]
    mean_change = dec_mean_df - c2
    
    res1 = dec_mean_change[dec_name].sum(axis=1).tolist()
    res2 = mean_change[val_name].sum(axis=1).tolist()
    
    dec_min = 100
    dec_max = 0
    val_min=100
    val_max=0
    incorrect_drugs = []
    for i in range(len(res1)):
        if res1[i] < dec_min:
            dec_min = res1[i]
        elif res1[i] > dec_max:
            dec_max = res1[i]
            
        if res2[i] < val_min:
            val_min = res2[i]
        elif res2[i] > val_max:
            val_max = res2[i]
        # detect if the estimated value is reasonable
        if res1[i]*res2[i]>=0: # First or Third quadrant
            pass
        else:
            incorrect_drugs.append(sorted(mean_change.index.tolist())[i])
    
    if do_plot:
        fig,ax = plt.subplots(figsize=(6,6),dpi=dpi)
        for i in range(len(res1)):
            plt.scatter(res1[i],res2[i],label=sorted(mean_change.index.tolist())[i],alpha=1.0)
        
        f_min = min(dec_min,val_min)
        f_max = max(dec_max,val_max)
            
        plt.plot([f_min,f_max],[f_min,f_max],linewidth=2,color='black',linestyle='dashed',zorder=-1)
        mergin=0.1
        r1 = patches.Rectangle(xy=(0,0),width=f_max+mergin,height=f_max+mergin,zorder=-2,color="paleturquoise",alpha=1.0) # rectangle object
        r2 = patches.Rectangle(xy=(0,0),width=f_min-mergin,height=f_min-mergin,zorder=-2,color="paleturquoise",alpha=1.0)
        ax.add_patch(r1)
        ax.add_patch(r2)
        
        plt.legend(shadow=True)
        
        plt.xlabel("Deconvolution estimated mean value")
        plt.ylabel("FACS mean value")
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        ax.set_axisbelow(True)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.title(str(dec_name)+" vs "+str(val_name))
        plt.show()
        print("False estimation :",incorrect_drugs)
    else:
        pass
    return incorrect_drugs

def plot_box_change(deconv_df,val_df,dec_name=["CD4","CD8"],val_name=["abT"],sort_index=["Ctrl","APAP","MDA","ANIT","TAA","ConA","GAL","CCl4"]):
    """
    Display deconvolution and FACS boxes side-by-side
    """
    # plot deconvolution result
    df = copy.deepcopy(deconv_df)
    df.index = [i.split("_")[0] for i in df.index]
    df = df.loc[sort_index]
    
    fig = plt.figure(figsize=(10,5))
    df_melt = pd.melt(df[dec_name].T)
    final_melt = pd.DataFrame()
    for t in sort_index:
        final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
    my_pal = {val: "yellow" if val in ["C","CIV","CIP"] else "lightgreen" for val in final_melt.variable.unique()}
    ax = fig.add_subplot(1,2,1)
    sns.boxplot(x='variable', y='value', data=final_melt, width=0.6, showfliers=False, notch=False, ax=ax, boxprops=dict(alpha=0.8), palette=my_pal)
    sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, color='black', ax=ax)
    ax.set_title("Deconvolution :"+str(dec_name))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.grid(False)
    plt.tight_layout()
    
    # plot FACS reference result
    df2 = copy.deepcopy(val_df)
    df2.index = [i.split("_")[0] for i in df2.index]
    df2 = df2.loc[sort_index]
    df_melt2 = pd.melt(df2[val_name].T)
    final_melt2 = pd.DataFrame()
    for t in sort_index:
        final_melt2 = pd.concat([final_melt2,df_melt2[df_melt2["variable"]==t]])
    my_pal = {val: "yellow" if val in ["C","CIV","CIP"] else "lightgreen" for val in final_melt2.variable.unique()}
    ax = fig.add_subplot(1,2,2)
    sns.boxplot(x='variable', y='value', data=final_melt2, width=0.6, showfliers=False, notch=False, ax=ax, boxprops=dict(alpha=0.8), palette=my_pal)
    sns.stripplot(x='variable', y='value', data=final_melt2, jitter=True, color='black', ax=ax)
    ax.set_title("FACS :"+str(val_name))
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.grid(False)
    plt.tight_layout()
    plt.show()

def plot_box(target_df,sort_index = ["C","AP","MDA","AN","TA","CA","GAL","C4"],row_n=2,col_n=3):
    df = copy.deepcopy(target_df)
    immunes = df.columns.tolist()
    df.index = [i.split("_")[0] for i in df.index]
    df = df.loc[sort_index]
    fig = plt.figure(figsize=(5*col_n,5*row_n))
    for i,immune in enumerate(immunes):
        df_melt = pd.melt(df[[immune]].T)
        final_melt = pd.DataFrame()
        for t in sort_index:
            final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
        my_pal = {val: "yellow" if val in ["C","CIV","CIP"] else "lightgreen" for val in final_melt.variable.unique()}
        ax = fig.add_subplot(row_n,col_n,i+1)
        sns.boxplot(x='variable', y='value', data=final_melt, width=0.6, showfliers=False, notch=False, ax=ax, boxprops=dict(alpha=0.8), palette=my_pal)
        sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, color='black', ax=ax)
        ax.set_title(immune)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
    plt.grid(False)
    plt.tight_layout()
    plt.show()

