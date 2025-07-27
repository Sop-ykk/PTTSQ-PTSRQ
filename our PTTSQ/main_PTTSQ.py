# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:48:47 2022

@author: YiKelai
"""

#%%
#### step 0：Data Owner 数据预处理

# 使轨迹点分布在[2^n,2^n]平面内
# 下面采用 order = 2n 的希尔伯特曲线来编码
import json
import numpy as np

with open('../traj_data/GT11.json') as file:
    data0 = json.load(file)
data0 = np.array(data0)
# data0[0][0][0]


# 等间隔抽稀
un = np.shape(data0)[0]     # user_num = 1100
tl = np.shape(data0)[1]     # traj_len = 50

n = 1     # un间隔n取值
l = 50      # tl间隔l取值
data = data0[0:un:n,0:tl:l]


#%%

import funR1_table as fR1
import funR2_tr2H as fR2
import funR5_enc as fR5

import random
import numpy as np

#### step 1：DO 建立映射表
order = 11
# os = 1

# k0=2048
# k1=24       # plaintext 介于 (-2**k1,2**k1)
# k2=160
# [pp, sk] = fR5.keygen(k0,k1,k2)
# E01 = fR5.enc(0, pp, sk)
# E02 = fR5.enc(0, pp, sk)
# pk = pp + [E01, E02]

## 生成映射表
# tb111 = fR1.gen_table2E(order,os,pk)
# np.save('tab111.npy', tb111)
# np.save('key111.npy', [pk,sk])

## 加载映射表，需要同时存下pk和sk才行
tab111 = np.load('../data/tab_11_1.npy', allow_pickle=True).item()
pk,sk = np.load('../data/key_11_1.npy', allow_pickle=True)
# tb12 = np.load('tb12.npy', allow_pickle=True).item()
#%%
#### step 2：DO 数据编码
trnH = fR2.traj_enc(data[0:,0:,0:2],order)


#### step 3：QU 查询编码
j = random.randint(0,len(data)-1)
data_Q = data[j]


trQH = fR2.traj_enc(data_Q[0:,0:2],order)



#%% CLOUD_I
#### step 4: 查表

import funR3_dist as fR3
import heapq

ke = 20
kh = 50

trnE = data[0:,0:,0:2]
trQE = data_Q[0:,0:2]

dne = fR3.traj_dist_Eu(trnE,trQE)   # 欧式距离 for verify
dnse = fR3.traj_dist_SEu(trnE,trQE)   # 欧氏平方距离 for verify
dnh = fR3.traj_dist_H1(trnH,trQH,tab111)


#### step 5: 初筛kh
import funR4_filter as fR4

ie = heapq.nsmallest(ke, range(len(dne)), dne.__getitem__) # 获取前ke个Eu距离最近的id
ise = heapq.nsmallest(ke, range(len(dnse)), dnse.__getitem__) # 获取前ke个Eu距离最近的id
ih = heapq.nsmallest(kh, range(len(dnh)), dnh.__getitem__) # 获取前kh个H1距离最近的id

pr = fR4.topk_rate1(ie, ih)      # 初筛查准率


#%% CLOUD_I  --  精确搜索
#### step 6: 获取密文轨迹坐标

import funR6_lookup as fR6
import numpy as np

# 这里原始索引改变了
trKEnc = []
for i in range(0,len(ih)):
    trKEnc.append(fR6.look_up_H0(trnH[ih[i]], tab111))
trKEnc = np.asarray(trKEnc)
    
trQEnc = fR6.look_up_H0(trQH, tab111)
    

#### step 7: 密文计算欧式距离、密文比较排序

import funR7_secpro as fR7


En1 = fR5.encp(-1, pk)      # E(-1)
Kseu = fR7.traj_dist_SumEu(trKEnc,trQEnc,pk[:4],En1)    # 各轨迹点距之和的密文list，没有排序
Klen = [np.shape(trKEnc)[1]] * len(Kseu)        # trKEnc中每一条轨迹的长度（可不一）
# Kdis = [Kseu,Klen]            # 合并列表，[轨迹点距之和,轨迹长度(点数)]

# [ls1,Kss] = fR7.sort_rrk(Kseu,Klen,pk[:4],sk) 
lss, Kss = fR7.sort_rrk(Kseu,Klen,pk[:4],sk)      # 排序后

#%% Test
# import funR5_enc as fR5

dcro = []
for i in range(10):
    dcro.append(fR5.dec(Kseu[i], pk[:4],sk))

drrk = []
for i in range(10):
    drrk.append(fR5.dec(Kss[i], pk[:4],sk))

deacr = []
for i in range(10):
    deacr.append(dne[ie[i]]) 
    
dneacr = []
for i in range(10):
    dneacr.append(dnse[ie[i]]) 
    
    
dhacr = []
for i in range(10):
    dhacr.append(dnh[ie[i]]) 
    
trQDec = []
for i in range(3):
    trQDec.append([fR5.dec(trQEnc[i][0], pk[:4],sk),fR5.dec(trQEnc[i][1], pk[:4],sk)]) 

trKDec = []
for i in range(3):
    trKDec.append([fR5.dec(trKEnc[1][i][0], pk[:4],sk),fR5.dec(trKEnc[1][i][1], pk[:4],sk)]) 
    
print('dcro:'+str(dcro)+',\ndrrk:'+str(drrk)+',\ndneacr:'+str(dneacr)+',\ndhacr:'+str(dhacr))
    

