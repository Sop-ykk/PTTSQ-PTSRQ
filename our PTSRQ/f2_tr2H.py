# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 15:24:06 2022

@author: YiKelai
"""
from hilbert import decode, encode
import numpy as np

#### 将各条轨迹坐标编码成Hil值
# 所有轨迹点位于ordN平面内，data元素为array

def traj_enc(data, ordN):
    H = encode(data, 2, ordN)
    Ha = H.reshape(np.shape(data)[:-1])

    return Ha

def traj_dec(Hdata, ordN):
    Ploc = decode(Hdata, 2, ordN)
    return Ploc
