# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 09:59:06 2022

@author: YiKelai
"""
import numpy as np
import math
from hilbert import decode, encode
import funR5_enc as fR5



#### 枚举每个数据点
def gen_allp(order):
    PH = []
    Ploc = []
    for i in range(2**(2*order)):
        PH.append(i)
    Ploc = decode(np.array(PH), 2, order)
    # PH = np.array(PH)
    # Ploc = np.array(Ploc,dtype = int)
    return PH,Ploc
# PH,Ploc = gen_allp(3)


#### 旋转函数
## P(x,y)绕Hil中心点, 顺时针旋转i*90°, 得到Pr(xr,yr)

# 1）Ploc to Ploc的旋转
def xy_rotate(P, order, i):
    [px,py] = P
    angle = i*math.pi/2
    oxy = np.array(2**(order-1)-0.5)
    Prx = round((px - oxy) * math.cos(angle) + (py - oxy) * math.sin(angle) + oxy)
    Pry = round((py - oxy) * math.cos(angle) - (px - oxy) * math.sin(angle) + oxy)
    Pr = [Prx, Pry]
    return Pr

# 3) Ploc to H的旋转
def nP2H_rotate(nP, order, ri):  # P是array
    Pr = []
    for j in range(len(nP)):
        Pr.append(xy_rotate(nP[j], order, ri))
    Hr = encode(np.array(Pr), 2, order)
    return Hr



####平移函数
## P(x,y)向右上方平移2^i个单位，得到Ps(xs,ys).
#  Ploc to H的平移，Ploc在order平面内，H在order+1平面内。0 <= si <= order
def nP2H_shift(nP, order, si):
    if si < order:
        Ps = nP + 2**si
        Hs = encode(np.array(Ps), 2, order+1)
        return Hs
    

#### 生成映射表
#0) ['H0','（明文）坐标']
def gen_table0(order):
    H0,Ploc = gen_allp(order)   # H0作为键的列表, Ploc作为值[0]的列表
    dic0 = dict(zip(H0,Ploc))                # 转换为字典
    return dic0
    
#1) ['H0','坐标','Hr1','Hr2','Hr3']
def gen_table1(order):
    H0,Ploc = gen_allp(order)   # H0作为键的列表, Ploc作为值[0]的列表
    Hr = []
    for i in range(1,4):        # 旋转90，180，270
        H0_ri = nP2H_rotate(Ploc, order, i)
        Hr.append(H0_ri)
    Hr = np.array(Hr)
    values = [Ploc,Hr]
    dic1 = dict(zip(H0,values))                # 转换为字典
    return dic1

# dic1 = gen_table1(3)

#2) 'H0': ['x坐标','y坐标','Hr1','Hr2','Hr3','Hs',...]
## P(x,y)向右上方平移1，2，4，2^os个单位，得到Ps(xs,ys).
#  Ploc to H的平移，Ploc在order平面内，H在order+1平面内。0 <= os是最大平移量<order
#  当os = order，成为了原型H0

# 键: H0
# 值[0,1]: Ploc
# 值[2,.]: Hrs

def gen_table2E(order,os,pk):      # to array
    H0,Ploc = gen_allp(order)   
    Hrs = []
    for ri in range(1,4):        # 旋转90，180，270
        H0_ri = nP2H_rotate(Ploc, order, ri)
        Hrs.append(H0_ri)
    for si in range(os+1):      # 向右上方平移1，2，4，2^os个单位
        H0_si = nP2H_shift(Ploc, order, si)
        Hrs.append(H0_si)  
    Hrs = np.array(Hrs)
    Hrs = np.transpose(Hrs)
    EPloc = fR5.enc_Ploc(Ploc, pk)
    values = np.hstack((EPloc,Hrs))             
    dic2 = dict(zip(H0,values))                # 转换为字典
    return dic2

# dic2 = gen_table2(4,2,pk)
