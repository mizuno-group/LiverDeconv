#!/usr/bin/env python3
"""
Created on 2023-10-16 (Mon) 15:30:29

@author: I.Azuma
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
liverdeconv = [np.mean([0.7549, 0.7554, 0.7614, 0.7625, 0.7625, 0.7656, 0.7691, 0.7785, 0.779, 0.7833]), np.mean([0.805, 0.8104, 0.8194, 0.8209, 0.8234, 0.8289, 0.8302, 0.8324, 0.8328, 0.8563]), np.mean([0.8628, 0.8629, 0.8632, 0.8653, 0.8666, 0.8667, 0.8687, 0.872, 0.8758, 0.8822]),np.mean([0.8478, 0.8519, 0.8582, 0.8637, 0.8645, 0.8672, 0.8688, 0.8745, 0.8783, 0.8882])]
elasticnet = [0.452, 0.322, 0.6291, 0.1312]
fardeep = [ 0.2357, 0.003, 0.3579, -0.1905]
epic = [0.4287, 0.5452, 0.3454, -0.1405]
cibersort = [0.2478, 0.2896, np.nan, -0.0617]
dcq = [0.3095, 0.0092, -0.0754, -0.4384]

data = [liverdeconv,elasticnet,fardeep,epic,cibersort,dcq]
methods = ['LiverDeconv','Elascit Net','FARDEEP','EPIC','CIBERSORT','DCQ']
cells = ['Neutrophils', 'Monocytes', 'NK', 'Eosinophils']

hatch_list = ['/', '|', '-', '+', 'x', 'o', 'O', '.', '*']
fig,ax = plt.subplots(figsize=(10,6),dpi=300)
x = [k*len(data) for k in range(len(cells))]
for i in range(len(data)):
    v = data[i]
    x2 = [t+i*0.8 for t in x]
    plt.bar(x2,v,label=methods[i],width=0.7,hatch=hatch_list[i]*3)
plt.xticks([t+1.6 for t in x],cells,fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel('Correlation',fontsize=18)
plt.legend(loc='upper right',shadow=True,fontsize=12,bbox_to_anchor=(1.1,1))
#plt.legend(loc='upper center',shadow=True,fontsize=13,ncol=6,bbox_to_anchor=(.5, 1.12))
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=1.0, axis='y')
plt.show()
# %%
