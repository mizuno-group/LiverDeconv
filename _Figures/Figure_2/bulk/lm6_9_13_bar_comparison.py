# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 13:03:14 2023

@author: I.Azuma
"""
import pandas as pd
import matplotlib.pyplot as plt

#%%
lm6 = [0.592, 0.516, 0.250]
lm9 = [0.299, 0.586, 0.555]
lm13 = [0.332, 0.486, 0.626]
data = [lm6,lm9,lm13]
methods = ['6 cells', '9 cells', '13 cells']

cells = ['Neutrophils', 'Monocytes', 'NK']

hatch_list = ['/', '|', '-', '+', 'x', 'o', 'O', '.', '*']
fig,ax = plt.subplots(figsize=(10,6),dpi=300)
x = [k*len(data) for k in range(len(cells))]
for i in range(len(data)):
    v = data[i]
    x2 = [t+i*0.8 for t in x]
    plt.bar(x2,v,label=methods[i],width=0.7,hatch=hatch_list[i]*3)
plt.xticks([t+0.8 for t in x],cells,fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel('Correlation',fontsize=18)
plt.legend(loc='upper right',shadow=True,fontsize=20,bbox_to_anchor=(1.5,1))
#plt.legend(loc='upper center',shadow=True,fontsize=13,ncol=6,bbox_to_anchor=(.5, 1.12))
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=1.0, axis='y')
plt.show()