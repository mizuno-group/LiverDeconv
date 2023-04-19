# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:25:03 2022

Plotting modules for deconvolution

@author: I.Azuma
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import copy

def plot_violin(target_df,sep="_",sort_index=["AP","MDA","AN","TA","CA","GAL","C4"],
                target_cell=[],ctrl="C",
                x_doc="",y_doc="score",row_n=0,col_n=0):
    """
    Violon plotting module for population analysis
    ----------
    target_df : DataFrame
        Samples are in rows and cell names are in columns.
                     B  Basophil       CD4  ...        NK  Neutrophil  Stellate
        AN_1 -0.029284 -0.000000 -0.043389  ... -0.114273    0.030369  0.083063
        AN_2 -0.022052  0.015056 -0.046337  ... -0.091490    0.022222  0.146941
        AN_3  0.000000  0.002150 -0.046066  ... -0.099084    0.009675  0.075869
        AN_4 -0.012023 -0.000000 -0.052731  ... -0.093115    0.021021  0.155818
        AP_1 -0.046533 -0.022409 -0.052776  ... -0.116688    0.076448  0.161657
        
    sep : str
        Key to distinguish sample names. The default is "_".
    sort_index : list
        Samples in interest. The default is ["AP","MDA","AN","TA","CA","GAL","C4"].
    target_cell : list
        Cells to be analyzed. The default is ["Monocyte","Neutrophil","NK","Eosinophil"].
    ctrl : str
        Sample name to determine the baseline. The default is "C".
    x_doc : str, optional
        Documents displayed in x axis. The default is "".
    y_doc : str, optional
        Documents displayed in y axis.. The default is "score".
    row_n : int, optional
        The number of subplots in row. The default is 0.
    col_n : int, optional
        The number of subplots in colums. The default is 0.

    """
    if len(target_cell) == 0:
        target_cell = target_df.columns.tolist()
    
    # optimize the subplots number
    if row_n == 0 or col_n == 0:
        row_n = int(np.sqrt(len(target_cell)))
        col_n =  len(target_cell)//row_n + 1
    
    df = copy.deepcopy(target_df)
    df.index = [i.split(sep)[0] for i in df.index] # rename the samples with sep condition
    fig = plt.figure(figsize=(5*col_n,5*row_n))
    for i,cell in enumerate(target_cell):
        df_melt = pd.melt(df[[cell]].T)
        # collect samples in interest (which are registered in sort_index)
        final_melt = pd.DataFrame()
        for t in sort_index:
            if t == ctrl:
                pass
            else:
                final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
        
        ax = fig.add_subplot(row_n,col_n,i+1) # identify the plotting location
        sns.violinplot(x='variable', y='value', data=final_melt, inner=None, cut=0, scale="count",linewidth=1.5)
        plt.setp(ax.collections, alpha=.55)
        sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, linewidth=1.5, size=9)
        #sns.swarmplot(x='variable', y='value', data=final_melt, ax=ax, linewidth=1, size=8)
        
        # plot the baseline of control samples
        ctrl_mean = df_melt[df_melt["variable"]==ctrl]["value"].mean() # mean value of ctrl samples (baseline)
        plt.hlines(y=ctrl_mean,xmin=-0.5,xmax=len(sort_index)-0.5,color="red",ls="dashed",linewidth=2)
        
        ax.set_title(cell,fontsize=15)
        ax.set_xlabel(x_doc,fontsize=13)
        ax.set_ylabel(y_doc,fontsize=13)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()

def plot_box(target_df,sep="_",sort_index=["AP","MDA","AN","TA","CA","GAL","C4"],
                target_cell=[],ctrl="C",
                x_doc="",y_doc="score",row_n=0,col_n=0):
    """
    Box plotting module for population analysis. Under construction...
    """
    if len(target_cell) == 0:
        target_cell = target_df.columns.tolist()
    
    # optimize the subplots number
    if row_n or col_n == 0:
        row_n = int(np.sqrt(len(target_cell)))
        col_n =  len(target_cell)//row_n + 1
        
    df = copy.deepcopy(target_df)
    df.index = [i.split(sep)[0] for i in df.index]
    fig = plt.figure(figsize=(5*col_n,5*row_n))
    for i,cell in enumerate(target_cell):
        df_melt = pd.melt(df[[cell]].T)
        final_melt = pd.DataFrame()
        for t in sort_index:
            final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
        my_pal = {val: "yellow" if val in ctrl else "lightgreen" for val in final_melt.variable.unique()}
        ax = fig.add_subplot(row_n,col_n,i+1)
        sns.boxplot(x='variable', y='value', data=final_melt, width=0.6, showfliers=False, notch=False, ax=ax, boxprops=dict(alpha=0.8), palette=my_pal)
        sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, color='black', ax=ax)
        ax.set_title(cell)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
    plt.grid(False)
    plt.tight_layout()
    plt.show()
    

def plot_errorbar(target_df,sep="_",sort_index=["AP","MDA","AN","TA","CA","GAL","C4"],
                target_cell=[],ctrl="C",
                x_doc="",y_doc="score",row_n=0,col_n=0):
    """
    Stripp plotting with errorbar for population analysis
    """
    if len(target_cell) == 0:
        target_cell = target_df.columns.tolist()
    
    # optimize the subplots number
    if row_n or col_n == 0:
        row_n = int(np.sqrt(len(target_cell)))
        col_n =  len(target_cell)//row_n + 1
        
    df = copy.deepcopy(target_df)
    df.index = [i.split(sep)[0] for i in df.index]
    fig = plt.figure(figsize=(5*col_n,5*row_n))
    for i,cell in enumerate(target_cell):
        df_melt = pd.melt(df[[cell]].T)
        
        ctrl_mean = df_melt[df_melt["variable"]==ctrl]["value"].mean() # mean value of ctrl samples (baseline)
        
        final_melt = pd.DataFrame()
        for t in sort_index:
            if t == ctrl:
                pass
            else:
                final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
        std = final_melt.groupby("variable").std().loc[sort_index] # standard deviation
        m = final_melt.groupby("variable").mean().loc[sort_index]
        ax = fig.add_subplot(row_n,col_n,i+1)
        sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, ax=ax, linewidth=1.5, size=9)
        ax.errorbar(x=[i for i in range(len(std))],y=m["value"].tolist(),yerr=std["value"].tolist(),capsize=10,marker='_',fmt='o',ms=20,mew=2,color="black")
        plt.hlines(y=ctrl_mean,xmin=-0.2,xmax=len(sort_index)-0.8,color="red",ls="dashed")
        
        ax.set_title(cell,fontsize=15)
        ax.set_xlabel(x_doc)
        ax.set_ylabel(y_doc)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
    plt.tight_layout()
    plt.show()

#%% Versatile Functions
def versatile_violin(target_df,sep="_",sort_index=["AP","MDA","AN","TA","CA","GAL","C4"],
                target_cell=[],x_doc="",y_doc="score",row_n=0,col_n=0):

    if len(target_cell) == 0:
        target_cell = target_df.columns.tolist()
    
    # optimize the subplots number
    if row_n or col_n == 0:
        row_n = int(np.sqrt(len(target_cell)))
        col_n =  len(target_cell)//row_n + 1
    
    df = copy.deepcopy(target_df)
    df.index = [i.split(sep)[0] for i in df.index] # rename the samples with sep condition
    fig = plt.figure(figsize=(5*col_n,5*row_n))
    for i,cell in enumerate(target_cell):
        df_melt = pd.melt(df[[cell]].T)
        # collect samples in interest (which are registered in sort_index)
        final_melt = pd.DataFrame()
        for t in sort_index:
            final_melt = pd.concat([final_melt,df_melt[df_melt["variable"]==t]])
        
        ax = fig.add_subplot(row_n,col_n,i+1) # identify the plotting location
        sns.violinplot(x='variable', y='value', data=final_melt, inner=None, scale="count",linewidth=1.5)
        plt.setp(ax.collections, alpha=.55)
        sns.stripplot(x='variable', y='value', data=final_melt, jitter=True, linewidth=1.5, size=9)
        #sns.swarmplot(x='variable', y='value', data=final_melt, ax=ax, linewidth=1, size=8)
        
        ax.set_title(cell,fontsize=15)
        ax.set_xlabel(x_doc,fontsize=13)
        plt.xticks(rotation=90)
        ax.set_ylabel(y_doc,fontsize=13)
        ax.grid(color="#ababab",linewidth=0.5)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()