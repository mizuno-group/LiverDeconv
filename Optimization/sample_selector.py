# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:20:29 2022

@author: I.Azuma
"""
import matplotlib.pyplot as plt
import seaborn as sns


from _utils import processing

def zscore_sampling(val_df,ctrl="C",target_cell=None,upper_z=1.0,lower_z=-1.0,outlier=3.0,centering=False,ignore_ctrl=True):
    """
    -outlier < x < lower_z
    upper_z < x < outlier

    """
    if target_cell is None:
        raise ValueError("!! Set target cells !!")
    # normalize the validation data (FACS) with ctrl values
    norm_val = processing.ctrl_norm(val_df.T,ctrl=ctrl)
    if ignore_ctrl:
        use = []
        for t in norm_val.index.tolist():
            if t.split("_")[0] == ctrl:
                pass
            else:
                use.append(t)
        wo_norm_val = norm_val.loc[use]
        z_norm_val = processing.standardz_sample(wo_norm_val) # remove ctrl samples
    else:
        z_norm_val = processing.standardz_sample(norm_val)
        # centering with Ctrl value
        if centering:
            z_norm_val = processing.ctrl_centering(z_norm_val)
    # plot the current situation
    #plot4select(z_norm_val,upper_z=upper_z,lower_z=lower_z)
    # re add ctrl samples (use in deconvolution output normalization)
    z_norm_val = processing.standardz_sample(norm_val)
    target_df = z_norm_val[[target_cell]]
    #self.target_sample = target_df[((target_df[self.target_cell]>upper_z)|(target_df[self.target_cell]<lower_z))].index.tolist()
    #self.target_sample = target_df.query('3.0 > Monocyte > 0.5 | -3.0 < Monocyte < -0.5')
    target_sample = target_df.query(str(outlier)+' > '+target_cell+' > '+str(upper_z)+' | '+str(-outlier)+' < '+target_cell+' < '+str(lower_z)).index.tolist()
    print("--- evaluation ---")
    print("target: ",target_cell)
    print(len(target_sample),"/",len(z_norm_val),"samples are analytical target")
    
    if len(target_sample) < 20:
        print("!! Attention! Sample size may be too small. !!")
    
    # reflect target sample
    #target_mix = mix_df[target_sample]
    #target_val = z_norm_val.T[target_sample]
    
    return target_sample, z_norm_val

def quartile_sampling(val_df,ctrl="C",target_cell=None,ignore_ctrl=True,threshold=1.5,q=[0, 0.25, 0.5, 0.75, 1]):
    if target_cell is None:
        raise ValueError("!! Set target cells !!")
    # normalize the validation data (FACS) with ctrl values
    norm_val = processing.ctrl_norm(val_df.T,ctrl=ctrl)
    if ignore_ctrl:
        use = []
        for t in norm_val.index.tolist():
            if t.split("_")[0] == ctrl:
                pass
            else:
                use.append(t)
        norm_val = norm_val.loc[use]
        #z_norm_val = processing.standardz_sample(wo_norm_val) # remove ctrl samples
    else:
        pass
    
    # plot box
    #sns.boxplot(data=norm_val)
    #plt.xticks(rotation=45)
    #plt.show()
    
    describe = norm_val.quantile(q=q)
    target_desc = describe[[target_cell]]
    second = float(target_desc.loc[0.25])
    third = float(target_desc.loc[0.75])
    interquartile = third-second
    upper_outer = threshold*interquartile + third       # upper outlier
    lower_outer = -threshold*interquartile + second     # lower outlier
    
    target_sample = norm_val.query(str(upper_outer)+' > '+target_cell+' > '+str(third)+' | '+str(lower_outer)+' < '+target_cell+' < '+str(second)).index.tolist()
    
    print("--- evaluation ---")
    print("target: ",target_cell)
    print(len(target_sample),"/",len(norm_val),"samples are analytical target")
    
    if len(target_sample) < 20:
        print("!! Attention! Sample size may be too small. !!")
    
    return target_sample, norm_val

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