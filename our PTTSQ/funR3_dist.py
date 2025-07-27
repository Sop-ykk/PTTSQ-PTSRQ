# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 20:53:38 2022

@author: YiKelai
"""

import numpy as np
from scipy.spatial.distance import cdist


#### 标准二维空间坐标系
# 1) 查询与多条轨迹之间的标准欧式距离
# trnE = data
# trQE = data_Q

def traj_dist_Eu(trnE,trQE):
    l = np.shape(trnE)[1]   # 点数
    distn = []
    if len(trQE) == l:
        for ui in range(len(trnE)):     
            disti = np.linalg.norm(trnE[ui].astype(int) - trQE.astype(int))     # 查询Q与第ui条轨迹数据的距离
            distn.append(disti)
    return distn

def traj_dist_SEu(trnE,trQE):
    l = np.shape(trnE)[1]   # 点数
    distn = []
    if len(trQE) == l:
        for ui in range(len(trnE)):     
            disti = np.sum((trnE[ui] - trQE)**2)     # 查询Q与第ui条轨迹数据的距离
            distn.append(disti)
    return distn

#### 初筛Hil法 根据table求两条H轨迹之间的距离：（假设对齐）欧式距离
# 1) 两条轨迹之间求距离
def tr2tr_dist_H1(trPH,trQH,tb1):
    l = len(trQH)
    dTr = []
    if len(trPH) == l:
        for h in range(l):
            d0 = abs(trPH[h].astype(int) - trQH[h].astype(int))     # 原始H0键的距离
            PHrs = tb1[trPH[h]][2:].astype(int)
            QHrs = tb1[trQH[h]][2:].astype(int)
            dP = min(d0,min(abs(PHrs - QHrs)))   # 修正d0
            dTr.append(dP)
        dist1 = np.mean(dTr)
    return dist1

# 2) 查询与多条轨迹之间的距离
# trnH[条数,点数]
def traj_dist_H1(trnH,trQH,tb1):
    l = np.shape(trnH)[1]   # 点数
    distn = []

    if len(trQH) == l:
        for ui in range(len(trnH)):   
            disti = tr2tr_dist_H1(trnH[ui],trQH,tb1)    # 查询Q与第ui条轨迹数据的距离
            distn.append(disti)
    return distn














