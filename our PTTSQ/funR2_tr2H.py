# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 15:24:06 2022

@author: YiKelai
"""
from hilbert import decode, encode
import numpy as np

#### 将各条轨迹坐标编码成Hil值
# 所有轨迹点位于order平面内，data元素为array

def traj_enc(data,order):
    H = encode(data, 2, order)
    Ha = H.reshape(np.shape(data)[:-1])
    return Ha
