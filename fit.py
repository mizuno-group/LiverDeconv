# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 06:04:20 2020

@author: I.Azuma, K.Morita
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet,Ridge,Lasso,HuberRegressor,LinearRegression
from sklearn.model_selection import GridSearchCV
from sklearn.svm import NuSVR
import scipy as sp
import statsmodels.api as sm
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

### fitting modules ###

"""
input : reference dataframe / mix dataframe / parameters

output : one result matrix

"""

# ElasticNet
def fit_ElasticNet(ref,dat,alpha=1,l1_ratio=0.05,max_iter=1e5,**kwargs):
    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=max_iter, tol=1e-5, random_state=None, fit_intercept=True)
    model.fit(ref,dat)
    #print(model.score(ref,dat))
    res_mat = pd.DataFrame(model.coef_,index=dat.columns, columns=ref.columns)
    return res_mat

# Ridge
def fit_Ridge(ref,dat,alpha=1,max_iter=1e5,**kwargs):
    """alpha is equal to ElasticNet's alpha * (1-l1_ratio)"""
    model = Ridge(alpha=alpha, max_iter=max_iter, tol=1e-5, random_state=None, fit_intercept=True)
    model.fit(ref,dat)
    #print(model.score(ref,dat))
    res_mat = pd.DataFrame(model.coef_,index=dat.columns, columns=ref.columns)
    return res_mat

# Lasso
def fit_Lasso(ref,dat,alpha=1,max_iter=1e5,**kwargs):
    """alpha is equal to ElasticNet's alpha * l1_ratio"""
    model = Lasso(alpha=alpha, max_iter=max_iter, tol=1e-5, random_state=None, fit_intercept=True)
    model.fit(ref,dat)
    #print(model.score(ref,dat))
    res_mat = pd.DataFrame(model.coef_,index=dat.columns, columns=ref.columns)
    return res_mat

# NuSVR
def fit_NuSVR(ref,dat,nu=0.5,C=1.0,max_iter=1e5,**kwargs):
    """
    nufloat, default=0.5
        An upper bound on the fraction of training errors and a lower bound of the fraction of support vectors. 
        Should be in the interval (0, 1]. By default 0.5 will be taken.
    """
    res_mat=[]
    ap = res_mat.append
    x = np.array(ref)
    for i in range(len(dat.columns)):
        y = np.array(list(dat.iloc[:,i]))
        model = NuSVR(nu=nu, C=C, kernel='linear')
        model.fit(x, y)
        ap(model.coef_[0])
    res_mat = pd.DataFrame(res_mat,index=dat.columns, columns=ref.columns)
    return res_mat   

# OLS (Ordinary Least Square)
def fit_OLS(ref,dat,**kwargs):
    A = np.array(ref)
    res_mat = pd.DataFrame(index=dat.columns, columns=ref.columns)
    for i in range(len(list(dat.columns))):
        b = np.array(dat.iloc[:,i])
        ols_model = sm.OLS(b,A)
        result = ols_model.fit()
        res_mat.iloc[i,:] = result.params
    return res_mat

# NNLS (Non Negative Least Square)
def fit_NNLS(ref,dat,**kwargs):
    A = np.array(ref)
    res_mat = pd.DataFrame(index=dat.columns, columns=ref.columns)
    for i in range(len(list(dat.columns))):
        b = np.array(dat.iloc[:,i])
        res_mat.iloc[i,:] = sp.optimize.nnls(A,b)[0]
    return res_mat

# LR (Linear Regression)
def fit_LR(ref,dat,**kwargs):
    res_mat = pd.DataFrame(index=dat.columns, columns=ref.columns)
    for i in  range(len(dat.T)):
        linear = LinearRegression()
        tmp = dat.iloc[:,i]
        linear.fit(ref,tmp)
        res_mat.iloc[i,:] = linear.coef_
    res_mat = res_mat.astype(float)
    return res_mat

# RLR (Huber Robut Linear Regression)
def fit_RLR(ref,dat,epsilon=1.35, max_iter=1e5, alpha=0.0001,**kwargs):
    """
    Linear regression model that is robust to outliers.
    epsilonfloat, greater than 1.0, default=1.35
        The parameter epsilon controls the number of samples that should be classified as outliers. 
        The smaller the epsilon, the more robust it is to outliers.
    alphafloat, default=0.0001
    """
    res_mat = pd.DataFrame(index=dat.columns, columns=ref.columns)
    for i in  range(len(dat.T)):
        huber = HuberRegressor(epsilon=epsilon, max_iter=max_iter, alpha=alpha,warm_start=False, fit_intercept=True)
        tmp = dat.iloc[:,i]
        huber.fit(ref,tmp)
        res_mat.iloc[i,:] = huber.coef_
    res_mat = res_mat.astype(float)
    return res_mat