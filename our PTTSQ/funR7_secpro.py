# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 17:28:59 2022

@author: YiKelai
"""


import funR5_enc as fR5
import random
import numpy as np

#### 密文计算欧式距离 ESED
# def ESED(EP1, EP2, pp, En1):
#     Edx = fR5.hsub(EP1[0], EP2[0], pp, En1)
#     Edy = fR5.hsub(EP1[1], EP2[1], pp, En1)
#     Edx2 = fR5.hmul(Edx, Edx, pp)
#     Edy2 = fR5.hmul(Edy, Edy, pp)
#     ESED = fR5.hadd(Edx2, Edy2, pp)
#     return ESED

#### 标准二维空间坐标系 - 密文

# 1) 两条轨迹之间求距离（所有点距平方之和，后面比较再取平均）
def Etr2tr_dist_SumEu(trPEnc,trQEnc,pp,En1):
    l = len(trQEnc)
    sdTr = 0
    if len(trPEnc) == l:
        for i in range(l):
            EP1 = trPEnc[i]
            EP2 = trQEnc[i]
            Edx = fR5.hsub(EP1[0], EP2[0], pp, En1)
            Edy = fR5.hsub(EP1[1], EP2[1], pp, En1)
            Edx2 = fR5.hmul(Edx, Edx, pp)
            Edy2 = fR5.hmul(Edy, Edy, pp)
            dP = fR5.hadd(Edx2, Edy2, pp)
            sdTr = fR5.hadd(sdTr, dP, pp)     # 点距之和
    return sdTr


# 2) 查询与多条轨迹之间的距离
# trKEnc[条数K,点数,2]
def traj_dist_SumEu(trKEnc,trQEnc,pp,En1):
    l = np.shape(trKEnc)[1]   # 点数
    distn = []
    if len(trQEnc) == l:
        for ui in range(len(trKEnc)):   
            disti = Etr2tr_dist_SumEu(trKEnc[ui],trQEnc,pp,En1)    # 查询Q与第ui条轨迹数据的距离
            distn.append(disti)
    return distn


#### 安全比较排序 s/n
# 密文E1,E2,E3,... (对应s1,s2), 另外轨迹长度明文n1,n2,n3,...; 安全排序s/n
# 交叉相乘比较密文, 若1<2, 则输出True，否则False
# Kseu, Klen为list
import gmpy2

def sort_rrk(Kseu,Klen,pp,sk):
    if len(Kseu) == len(Klen):
        # ls = list(range(len(Kseu)))
        n = pp[0]
        ## 1) 密文扰动
        ra = random.randint(0,2**pp[2])
        rb = random.randint(-2**pp[2],2**pp[2])
        Kseur = []
        for i in range(len(Kseu)):
            x1 = gmpy2.mod(gmpy2.mul(ra, Kseu[i]), n)
            x2 = gmpy2.mod(gmpy2.mul(rb, Klen[0]), n)
            Kseur.append(gmpy2.mod(gmpy2.add(x1, x2), n))
        ## 2) 解密扰动密文
        Ksr = []
        for j in range(len(Kseu)):
            Ksr.append(fR5.dec(Kseur[j], pp, sk))    # 解密所有Kseur
        avgKsr = np.asarray(Ksr)/Klen[0]
            
        ## 3) 所有扰动明文的排序
        # for i in range(0,len(Kseu)):
        #     ci = i      # current_index
        #     while avgKsr[ci] < avgKsr[ci-1] and ci - 1 >= 0:
        #         ls[ci], ls[ci-1] = ls[ci-1], ls[ci]
        #         ci -= 1
        ass = np.sort(avgKsr)
        lss = np.argsort(avgKsr)
        
        ## 4) 用户消除扰动，按ls的index排序
        # Kss = [0] * len(Kseu)
        # for i in range(0,len(Kseu)):
        #     Kss[ls[i]] = (np.asarray(Ksr[i])-rb)/ra
        Kss = (np.asarray(ass)-rb)/ra
        return lss,Kss


#%%

# En1 = fR5.encp(-1, pk)      # E(-1)

# EP1 = trQEnc[0]
# EP2 = trQEnc[1]
# esed = ESED(EP1, EP2, pp, En1)
# sed = fR5.dec(esed, pp, sk)

# etr2trd = traj_dist_Eu(trKEnc,trQEnc,pp,En1)
# sed = fR5.dec(etr2trd[0],pp,sk)

#%%   Test
# pp = pk[:4]
# e3=fR5.enc(3, pp, sk)
# e5=fR5.enc(5, pp, sk)
# e7=fR5.enc(7, pp, sk)
# e9=fR5.enc(9, pp, sk)
# e11=fR5.enc(11, pp, sk)

# d3=fR5.dec(3, pp, sk)
# d5=fR5.dec(5, pp, sk)
# d7=fR5.dec(7, pp, sk)
# d9=fR5.dec(9, pp, sk)
# d11=fR5.dec(11, pp, sk)

# TKseu = [e3,e7,e11,e9,e5]
# TKlen = [3,3,3,3,3]

# Tls1 = sort_cro(TKseu,TKlen,pp,sk)
# # Tls2,TKsseu2 = sort_rrk(TKseu,TKlen,pp,sk)

# # print(Tls1,Tls2)
# # print(str(TKsseu1)+'\n'+str(TKsseu2))