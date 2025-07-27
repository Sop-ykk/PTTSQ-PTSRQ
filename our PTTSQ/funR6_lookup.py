# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:10:56 2022

@author: YiKelai
"""

#### 根据H0值，查表得到密文

# 一条轨迹trH (l,)

import numpy as np

def look_up_H0(trH, tb):
    trEnc = []
    for i in range(len(trH)):
        trEnc.append(tb[trH[i]][:2])
    trEnc = np.asarray(trEnc)
    return trEnc
    
        
    