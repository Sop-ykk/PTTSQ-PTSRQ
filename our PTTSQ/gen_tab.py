# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 09:21:48 2022

@author: YiKelai
"""

## 生成映射表
import funR1_table as fR1
# import funR2_tr2H as fR2
import funR5_enc as fR5

# import random
import numpy as np

#### step 1：DO 建立映射表
# order = 11      # order from 8 to 11
# os = 1          # os \in [0,order-1]

# k0=2048
# k1=24       # plaintext 介于 (-2**k1,2**k1)
# k2=160
# [pp, sk] = fR5.keygen(k0,k1,k2)
# E01 = fR5.enc(0, pp, sk)
# E02 = fR5.enc(0, pp, sk)
# pk = pp + [E01, E02]

## 生成映射表
# tb1 = fR1.gen_table2E(order,os,pk)
# np.save('tab2048111.npy', tb1)
# np.save('key2048111.npy', [pk,sk])

tab111 = np.load('tab111.npy', allow_pickle=True).item()
key111 = np.load('key111.npy', allow_pickle=True)
[pk, sk] = key111
