# 已经【等长】的DO、QU轨迹: 一条p_qu, 一条p_do

import numpy as np
from CKKS_HE import *
from matplotlib import pyplot as plt

# 单条轨迹之间每个点对的距离密文向量
def sCKKS(ep1_qu, ep1_do):  # 此处为平方欧式距离，最好需要拟合一下开方
    edv = []
    for i in range(len(ep1_qu)):
        evxdi = Bootstrap_vec(ep1_qu[i][0] - ep1_do[i][0]) 
        evydi = Bootstrap_vec(ep1_qu[i][1] - ep1_do[i][1]) 
        evdi2 = Bootstrap_vec(evxdi ** 2 + evydi ** 2)
        edv.append(evdi2)
        # evdi2 = Bootstrap(evdi2)
        # edv.append(evdi2 @ np.ones((2,1),int))
    return edv

# 单条轨迹之间的DSED距离加权求和的密文值（未取时间\delta平均）
def DSED_ckks(ep1_qu, ep1_do, t1_qu):   # 因为后面只需比较大小，总时间间隔均相同，所以eps0同时乘以分母。
    es = sCKKS(ep1_qu, ep1_do)
    h = len(t1_qu) - 1
    esd2 = (t1_qu[1] - t1_qu[0]) * es[0] + (t1_qu[h] - t1_qu[h - 1]) * es[h]     # 首尾的和式
    for k in range(1, h):
        # print(k)
        esd2 += (t1_qu[k + 1] - t1_qu[k - 1]) * es[k]
    # dsed = esd2 / (2 * (t1_qu[h] - t1_qu[0]))
    esd2 = Bootstrap_vec(esd2) 
    return esd2


# K条轨迹集合与某条查询轨迹的sumed_DSED距离向量
def KE_SD(ep1_qu, epK_do, t1_qu):
    edK = []
    for i in range(len(epK_do)):
        edK.append(DSED_ckks(ep1_qu, epK_do[i], t1_qu))
    return edK


# -------------------------
### For verify, 解密某条轨迹密文list
# -------------------------
def Dec_tr1(ep1):
    p = []
    for i in range(len(ep1)):
        p.append(Dec_vec(ep1[i])[0])
    p = np.array(p).astype(int)
    return p

def Plot_etr1(etr):
    tr1 = Dec_tr1(etr)
    plt.plot(tr1[:,0],tr1[:,1],'o-',color='red',label='Q')

def Plot_etrn(ep1_qu, epn_do):
    tr1 = Dec_tr1(ep1_qu)
    plt.plot(tr1[:,0],tr1[:,1],'o-',color='red',label='Q')
    for i in range(len(epn_do)):
        tri = Dec_tr1(epn_do[i])
        plt.plot(tri[:,0],tri[:,1],'o-',label='T'+str(i))
    plt.legend()

# -------------------------
### 密文比较排序
# -------------------------
def f1(ex):
    ex = Bootstrap_vec(ex)
    ex2 = Bootstrap_vec(ex ** 2)
    ex3 = Bootstrap_vec(ex * ex2)
    ef1 = (3 * ex - ex3) * 0.5
    return ef1


def g1(ex):
    ex = Bootstrap_vec(ex)
    ex2 = Bootstrap_vec(ex ** 2)
    ex3 = Bootstrap_vec(ex * ex2)
    eg1 = (-1359 * ex3 + 2126 * ex) * (2 ** (-10))
    return eg1


def leq_e1e2(ex):
    '''
    密文判定符号 ex = e1 - e2 <= 0 ?
    :param ex: 待判定符号的密文，其对应明文介于 (-1,1)，是Bootstrap过的新鲜密文
    :return: sgn(x)
    '''
    # nf: int=1, ng: int=1, df: int=8, dg: int=7
    # g1(x) = (-1359 * x ** 3 + 2126 * x) * (1/(2 ** 10))
    # f1(x) = (3*x - x**3) * 0.5
    # E(sgn(x)) = f1^8_g1^7(Ex)
    em = g1(ex)
    for i in range(6):
        em = g1(em)
    for j in range(8):
        em = f1(em)
    return em

def scl(edK):
    '''

    :param edK: K个相似度密文的[列表]
    :return: 1/最大差值, s.t. (edi-edj) * Δ 均介于[-1,1]
    '''
    dK = Dec_tr1(edK)
    delta = 1/(np.ptp(dK)*1.001)       # ptp: 最大最小值之差
    return delta


def scl_vec(evK):
    '''

    :param evK: 一个相似度密文的K维[向量]
    :return: 1/最大差值, s.t. (edi-edj) * Δ 均介于[-1,1]
    '''
    dK = Dec_vec(evK)
    delta_v = 1/(np.ptp(dK)*1.001)       # ptp: 最大最小值之差
    return delta_v


def rk_esd2(edK):
    '''

    :param edK: K个相似度密文的列表
    :return: lsi: K维关系矩阵（上三角、不包含对角线）的压缩列表; rk: K条轨迹的相似度排名
    '''
    delta = scl(edK)
    lsi, rk = [], []
    rkn = 0     # 最后一位，前面项之和的相反数
    for i in range(len(edK)-1):
        rki, lsj = 0, []
        for j in range(i+1, len(edK)):
            # print(i,j)
            edj = (edK[i] - edK[j]) * delta
            es = leq_e1e2(edj)      # 介于[-1, 1]
            lsj.append(es)
            rki += es
        lsi.append(lsj)
        for k in range(1,i+1):
            rki -= lsi[k-1][i-k]
        rk.append(rki)
        rkn -= rki
    rk.append(rkn)
    return rk
    # return lsi, rk


def rk_esd(edK, eps0):
    '''
    
    :param edK: K个相似度密文的[列表]
    :return: lsf: 长度为K的符号密文列表, E2代表符合范围, E0代表相似度在范围之外
    '''
    delta = scl(edK)
    lsf = []
    for i in range(len(edK)):
        edj = (edK[i] - eps0) * delta
        ef = (leq_e1e2(edj) + 1)      # 介于[0, 2], 本应再除以2, 但是报错需要bootstrap, 因此省去。事实上，在恢复结果时再除也无影响。
        lsf.append(ef)
    return lsf


def rv_esd(edK, eps0):
    '''
    
    :param evK: 一个相似度密文的K维[向量]
    :return: lsf: 一个符号密文的K维向量, E2代表符合范围, E0代表相似度在范围之外
    '''
    evK = Epack_ls2vec(edK)
    delta = scl_vec(evK)
    edv = (evK - eps0) * delta
    vecf = (leq_e1e2(edv) + 1)
    return vecf


def HE_vecmul3(vecf, hK_do, EidK_do, vm1, tl):
    '''
    前面的编码还是按照每条轨迹内部独立进行打包
    用于同态乘法vecf, 返回掩码后的轨迹编码结果rs_trH, ID结果rs_id
    密文向量vecf首先进行拆分
    :param hK_do: 初筛之后的轨迹编码[K, tl]维矩阵, 即, hn_do_[ir_Pv], 将做H明文矩阵的同态乘法
    :param EidK_do: 初筛之后的轨迹ID密文列表, 长度为K, 将做K次Eid密文向量的同态乘法
    :param vm1: 初筛轨迹的条数K
    :param vm2: ID字符串长度, 默认50
    :return: rs_trH: 
    :return: rs_Eid: 

    '''
    rs_trH, rs_Eid = [], []
    for i in range(tl):
        rs_trH.append(vecf * (hK_do[:,i]))
    I = np.eye(vm1)
    # lsf = Esubs_vec2ls(vecf, vm1)
    for i in range(vm1):
        Ei = Bootstrap_vec(EidK_do[i])
        vecf_Ie = vecf @ (I[i].reshape((-1,1)))
        vecf_Ie_ = Bootstrap_vec(vecf_Ie)
        rs_Eid.append(vecf_Ie_ * Ei)
    return rs_trH, rs_Eid

    lsf = []
    I = np.eye(vm)
    for i in range(vm):
        lsf.append(vecf @ (I[i].reshape((-1,1))))


if __name__ == '__main__':
    edK = KE_SD(ep1_qu, epK_do, t1_qu)
    Dec_tr1(edK)
    rk = rk_esd2(edK)
    Dec_tr1(rk)
    lsf = rk_esd(edK, eps0)
    Dec_tr1(lsf)
    vecf = rv_esd(edK, eps0)
    Dec_vec(lsf)
