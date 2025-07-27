# 一条轨迹trH (l,)
import time
from CKKS_HE import *
import numpy as np
import math
import os
from hilbert import decode, encode
from s1_table import load_arHrs_, gen_allp
from s2_CKKSvec2ten import load_EPloc_sr, load_context


def loc_xy(hi, vmax:int = 4096):
    '''

    :param hi: 从0开始，0,1,2,...
    :param ordN:
    :param vmax:
    :return: hi所在的exi,eyi密文向量，vri为hi在该向量的位置
    '''
    # e_num = np.shape(Ev_sr)[1]
    vqi = hi // vmax       # 商数：hi位于第几个向量，从0开始，0,1,2,...
    vri = hi % vmax        # 余数：hi位于第vqi个向量的第几位，从0开始，0,1,2,...
    return vqi, vri


def read_xy(Ev_sr, ctx, hi, vmax:int = 4096):
    '''

    :param Ev_sr:
    :param hi: 从0开始，0,1,2,...
    :param ordN:
    :param vmax:
    :return: hi所在的exi,eyi密文向量，vri为hi在该向量的位置
    '''
    # e_num = np.shape(Ev_sr)[1]
    vqi, vri = loc_xy(hi, vmax)
    evxi = ts.ckks_vector_from(ctx, Ev_sr[0][vqi])
    evyi = ts.ckks_vector_from(ctx, Ev_sr[1][vqi])
    # li = one_hot(vri, vmax)
    li = [0] * vmax
    li[vri] = 1     # 仅vri位为1，独热编码
    exi = evxi.dot(li)
    eyi = evyi.dot(li)
    return exi, eyi


def look_up_tr1(h1, Ev_sr, ctx, vmax:int = 4096):
    '''

    :param h1:一条轨迹的Hil向量
    :param Ev_sr:
    :param ctx:
    :param vmax:
    :return: 长为tlq的list，每个元素为(exi, eyi)
    '''
    trEnc = []
    for i in range(len(h1)):
        trEnc.append(read_xy(Ev_sr, ctx, h1[i], vmax))
    return trEnc


# 只截取对齐后的DO轨迹长度
def look_up_trK(hn_do, ih, Ev_sr, ctx, vmax:int = 4096):
    trKEnc = []
    for i in range(0, len(ih)):
        trKEnc.append(look_up_tr1(hn_do[ih[i]], Ev_sr, ctx, vmax))
    return trKEnc


if __name__ == '__main__':
    # t1 = time.time()
    # trKEnc = look_up_trK(hn_do, ih, EPloc)
    # t2 = time.time()
    # print(t2-t1)

    from s0_data_pre import data_pre

    ordN = 8
    vmax = 4096
    all = False
    un = 200
    tl = 100
    tlq = 15
    tmax = 600
    # hi = h1_qu
    hi = 14
    ctx = load_context()
    Ev_sr = load_EPloc_sr(ordN, ctx)

    # vqi, vri = loc_xy(hi, vmax)
    # exi, eyi = read_xy(Ev_sr, ctx, hi)
    #
    # x1 = np.round(Dec_vec(exi)).astype(int)[vri]
    # y1 = np.round(Dec_vec(eyi)).astype(int)[vri]

    newfile = 'GT8_tall600'
    hn_do, h1_qu, t1_qu, pn_do, p1_qu = data_pre(newfile, ordN, tlq, all, un, tl)

    ih = [1,3,5,4,2]

    t1 = time.time()
    trKEnc = look_up_trK(hn_do, ih, Ev_sr, ctx, vmax)
    t2 = time.time()
    print(t2-t1)