# coding:utf-8

from CKKS_HE import *
import numpy as np
import math
import time
import os
# from hilbert import decode, encode
from s1_table import load_arHrs_, gen_allp


def save_context():
    # with open('./map_data/context.txt', 'wb') as f:
    #     f.write(gencontext().serialize(True,True,True,True))
    context = gencontext()
    np.save('./map_data/context.npy', context.serialize(True,True,True,True))
    return context


def load_context():
    # with open('./map_data/context.txt', 'rb') as f:
    #     ctx_sr = f.read()
    ctx_sr = np.load('./map_data/context.npy', allow_pickle=True).item()
    ctx = ts.context_from(ctx_sr)
    return ctx


def gen_EPloc(Ploc, ctx, vmax:int = 4096):
    '''
    :return: 密文矩阵 构建时间、空间(MB)
    '''
    t1 = time.time()
    l = len(Ploc)
    e_num = math.ceil(l / vmax)     # 向上取整，e_num个vx, e_num个vy
    evx, evy = [], []
    for ki in range(e_num):
        vxi = Ploc[:, 0][ki * vmax:(ki + 1) * vmax]
        vyi = Ploc[:, 1][ki * vmax:(ki + 1) * vmax]
        evxi = Enc_vec(ctx, vxi)
        evyi = Enc_vec(ctx, vyi)
        evx.append(evxi.serialize())
        evy.append(evyi.serialize())
    ev_sr = [evx, evy]
    np.save('./map_data/Evec_' + str(l) + '_' + str(e_num) + '.npy', ev_sr)
    t2 = time.time()
    te = t2 - t1
    seb = os.path.getsize('./map_data/Evec_' + str(l) + '_' + str(e_num) + '.npy')   # 以字节为单位
    sem = round(seb / (1024 * 1024), 2)
    with open('./map_data/enc_time_size.csv', 'a') as f:
        # f.write('映射字典耗时\n' + str(td) + '\n加密矩阵耗时\n' + str(te))
        f.write('加密向量组,' + str(l) + ',' + str(e_num)+',time_s,'+str(te)+',size_mb,'+str(sem)+'\n')
    # if (Ploc == np.round(Dec_ten(EPloc))).all():    # 检查Ploc加密取整是否正确
    return te, sem


def load_EPloc_sr(ordN, ctx, vmax:int = 4096):
    '''
    :return: 密文矩阵
    '''
    l = 2 ** (2*ordN)
    e_num = math.ceil(l / vmax)     # 向上取整，e_num个vx, e_num个vy
    evcpath = './map_data/Evec_' + str(l) + '_' + str(e_num) + '.npy'
    if not os.path.exists(evcpath):
        print("... [1] gen_EPloc ...")
        Ploc = gen_allp(ordN)[1]
        gen_EPloc(Ploc, ctx, vmax)
    Ev_sr = np.load('./map_data/Evec_' + str(l) + '_' + str(e_num) + '.npy', allow_pickle=True)
    # [Evx, Evy] = Ev_sr
    # for ki in range(e_num):
    #     ts.ckks_vector_from(ctx, Evxi_sr)
    return Ev_sr




if __name__ == "__main__":
    ordN, dro, dso = 11, 3, 3
    hi = 2001
    arHrs_, Ploc = load_arHrs_(ordN, dro, dso) 

    # ctx = load_context()
    # [Evx_sr, Evy_sr] = load_EPloc_sr(ordN)
    # Evx1 = ts.ckks_vector_from(ctx, Evx_sr[0])
    # vx1 = np.round(Dec_vec(Evx1)).astype(int)
