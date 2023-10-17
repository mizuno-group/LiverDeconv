# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 13:03:14 2023

@author: I.Azuma
"""
# %%
import pandas as pd
import matplotlib.pyplot as plt

#%%
lm5 = [0.573, 0.334, 0.451]
lm10 = [0.558, 0.419, 0.266]

data = [lm5,lm10]
methods = ['5 cells', '10 cells']

cells = ['Neutrophils', 'Monocytes', 'NK']

hatch_list = ['/', '|', '-', '+', 'x', 'o', 'O', '.', '*']
fig,ax = plt.subplots(figsize=(10,6),dpi=300)
x = [k*len(data) for k in range(len(cells))]
for i in range(len(data)):
    v = data[i]
    x2 = [t+i*0.8 for t in x]
    plt.bar(x2,v,label=methods[i],width=0.7,hatch=hatch_list[i]*3)
plt.xticks([t+0.4 for t in x],cells,fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel('Correlation',fontsize=18)
plt.legend(loc='upper right',shadow=True,fontsize=20,bbox_to_anchor=(1.1,1))
#plt.legend(loc='upper center',shadow=True,fontsize=13,ncol=6,bbox_to_anchor=(.5, 1.12))
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
ax.grid(color="#ababab",linewidth=1.0, axis='y')
plt.title('Elastic Net',fontsize=20)
plt.show()
# %%
