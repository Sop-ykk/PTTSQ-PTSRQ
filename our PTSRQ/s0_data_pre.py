import json
import numpy as np
import random
from f2_tr2H import traj_enc

def file_data(filename: str, all: bool = True, un=0, tl=0):
    # 读取文件中的轨迹信息 (x,y), un条数, tl长度，默认读全
    # return 三维数组：(条数, 长度, 2/4)
    with open("./traj_data/" + filename + ".json") as f1:
        data0 = json.load(f1)
    data = np.array(data0)
    if all:
        un, tl = np.shape(data)[:2]
    data = data[: un, :tl]
    return data

def gen_tsp(data, asyn: int = 3, tmax: int = 24):
    '''
    :param data: 轨迹数据库 (条数, 长度, 2)
    :param tmax: 时间最大值
    :param asyn: 1 异步采样 / 2 同步采样 / 3 同步全部时间戳
    :return: 每个个体返回一列递增的时间戳: (条数, 长度, 1)
    '''
    un, tl = np.shape(data)[:2]
    if asyn == 1:
        # tsp = np.random.randint(tmax, size=(un, tl))      # 时间戳可重复
        # tsp = np.sort(tsp)
        tsp = np.empty((un, tl), int)
        for i in range(0, un):
            tsi = random.sample(range(1,tmax+1), tl)      # 时间戳无重复
            tsp[i] = np.sort(tsi)
    elif asyn == 2:
        # ts = random.randint(tmax, size=tl)        # 时间戳可重复
        ts = random.sample(range(1, tmax+1), tl)      # 时间戳无重复
        ts = np.sort(ts)
        tsp = np.tile(ts, un)
    elif asyn == 3:
        ts = np.arange(1, tl+1)
        tsp = np.tile(ts, un)
    tsp1 = np.reshape(tsp, (un, tl, 1))
    return tsp1


def gen_id(data):
    '''
    :param data: 轨迹数据库 (条数, 长度, *)
    :return: 每个个体返回一列相同id: (条数, 长度, 1)
    '''
    un, tl = np.shape(data)[:2]
    num = np.arange(1, un+1)
    gid = np.repeat(num, tl)
    gid1 = np.reshape(gid, (un, tl, 1))
    return gid1


def data_tsp(filename: str, tmax, all: bool = True, un=0, tl=0, asyn: int = 3):
    data = file_data(filename)
    if all:
        un, tl = np.shape(data)[:2]
    data = data[: un, :tl]
    tsp = gen_tsp(data, asyn, tmax)
    gid = gen_id(data)
    dt = np.c_[data, tsp, gid]
    return dt

# def load_tab(order=11, os=1):
#     ## 加载映射表，需要同时存下pk和sk
#     tab = np.load('./tab_data/tab_'+str(order)+'_'+str(os)+'.npy', allow_pickle=True).item()
#     pk, sk = np.load('./tab_data/key_'+str(order)+'_'+str(os)+'.npy', allow_pickle=True)
#     return tab, pk, sk


def gen_qry(tlq, ordN, tmax: int = 600):
    trQE = np.random.randint(0,2**ordN-1, size=[tlq,2])
    trQE = np.reshape(trQE, (1, tlq, 2))
    tsq = gen_tsp(trQE, 2, tmax)
    id0 = np.repeat(0, tlq)
    id = np.reshape(id0, (1, tlq, 1))
    trQ = np.c_[trQE, tsq, id][0]
    return trQ


def save_data(rawfile:str, tmax, all: bool = True, un=0, tl=0, asyn: int = 3):
    data_new = data_tsp(rawfile, tmax, all, un, tl, asyn)
    if asyn == 1:
        tstype = '_tasyn'
    elif asyn == 2:
        tstype = '_tsyn'
    elif asyn == 3:
        tstype = '_tall'
    outpath = "./traj_data/" + rawfile + tstype + str(tmax) + ".json"
    with open(outpath, 'w') as f:
        json.dump(data_new.tolist(), f)


def read_gen(newfile: str, tlq, ordN, all: bool = True, un=0, tl=0):
    data_new = file_data(newfile, all, un, tl)
    if all:
        tl = np.shape(data_new[0:, 0:, 2])[1]
    assert tl >= tlq, "YKL sets \'tl>=tlq\', please adjust tl or tlq"
    tlq = min(tlq, np.shape(data_new)[1])
    tmax = np.max(data_new[0:, 0:, 2])
    trQ = gen_qry(tlq, ordN, tmax)
    return data_new, trQ


def data_pre(newfile, ordN, tlq, all: bool = True, un=0, tl=0):
    '''

    :param newfile: 数据文件名
    :param ordN: 曲线阶数
    :param un: 轨迹条数 user number（截取数据集）
    :param tl: 轨迹长度 traj length（截取数据集）
    :param tlq: 查询长度
    :return: hn_do, h1_qu, t1_qu
    '''
    # Step 1: 数据准备
    data, trQ = read_gen(newfile, tlq, ordN, all, un, tl)
    pn_do_ = data[:, :, 0:2]
    t1_do = data[:, :, 2][0]  # 假设DO时间戳都齐全且一样
    p1_qu = trQ[0:, 0:2]
    t1_qu = trQ[0:, 2]

    # step 23：DO 数据编码、QU 查询编码
    hn_do_ = traj_enc(pn_do_, ordN).astype(int)
    h1_qu = traj_enc(p1_qu, ordN).astype(int)

    # step 45: CS1 查表、初筛kh
    tid = np.nonzero(np.in1d(t1_do, t1_qu))[0]      # 匹配的tlq个时间戳，[tlq,]
    hn_do = hn_do_[:, tid]      # 取出所有轨迹在所查询的tlq个时间戳下的轨迹点h，[un, tlq]
    pn_do = pn_do_[:, tid]      # 取出所有轨迹在所查询的tlq个时间戳下的轨迹点p，[un, tlq,2]
    return hn_do, h1_qu, t1_qu, pn_do, p1_qu


def gen_ID(un):
    '''
    :param un: 轨迹数据库的条数
    :return gID: 一个形状为(len(ID_num), un)的身份num矩阵
    //一般而言, len(ID_num) < un < 4096, 因此加密时按un方向打包为向量, 即, gID[0], gID[1], ...
    但是, 过滤步骤需要按un拆分, 提取出其中的K个用户, 因此还是只能按上述转置方向打包向量
    '''
    ID_plain = 'Kelai YI, 510108190000000000, 18000000000, xxxxxxx'
    # 将字母字符串转换为数字表示  
    ID_num = [ord(char) for char in ID_plain]
    # gID1 = np.repeat(ID_num, un)
    # gID = np.reshape(gID1, (len(ID_num), un))
    gID1 = np.tile(ID_num, un)
    gID = np.reshape(gID1, (un, len(ID_num)))
    return gID


def rec_ID(ID_arr):
    '''
    :param ID_num: 一位用户的ID信息ASCII码的向量
    :return ID_plain: 一条un维的身份列表
    '''
    ID_plain = []
    # 将数字表示还原为字母字符串
    for i in range(len(ID_arr)):
        # ID_letters_i = [chr(num) for num in ID_ls[i]]
        ID_letters_i = ''.join(chr(x) if x else "" for x in ID_arr[i])
        ID_plain.append(ID_letters_i)
    return ID_plain


if __name__=='__main__':
    # rawfiles = ['GT8','GT9']
    # # newfiles = ['GT11_tall600']
    # tmax = 600
    # tlq = 18
    # ordN = 11
    # for rawfile in rawfiles:
    #     save_data(rawfile, tmax)

    # for newfile in files:
        # data = data_tsp(file, tmax)
        # data = data_tsp(file, tmax, True, 0, 0, 1)
        # data = data_tsp(file, 600, False, 10, 5, 1)
        # data = data_tsp(file, 600, False, 10, 24)
        # data_new, trQ = read_gen(newfile, tlq, ordN, all: bool = True, un=0, tl=0)
        # data_new, trQ = read_gen(newfile, tlq, ordN)

    ordN = 8
    all = True
    # un = 200
    # tl = 100
    tlq = 15
    tmax = 600
    save_data('GT'+str(ordN), tmax, True, un = 0, tl = 0)

    # newfile = 'GT11_tall600'
    # hn_do, h1_qu, t1_qu, pn_do, p1_qu = data_pre(newfile, ordN, tlq, all, un, tl)

