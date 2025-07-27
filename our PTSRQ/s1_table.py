import numpy as np
import math
from hilbert import decode, encode
from CKKS_HE import *
import time
import os

#### 枚举每个数据点
def gen_allp(ordN):
    PH = []
    Ploc = []
    for i in range(2 ** (2 * ordN)):
        PH.append(i)
    Ploc = decode(np.array(PH), 2, ordN).astype(int)
    # PH = np.array(PH)
    # Ploc = np.array(Ploc,dtype = int)
    return PH, Ploc

#### 旋转函数
## P(x,y)绕Hil中心点, 顺时针旋转i*90°, 得到Pr(xr,yr)

# 1）Ploc to Ploc的旋转
def xy_rotate(P, ordN, i):
    [px, py] = P
    angle = i * math.pi / 2
    oxy = np.array(2 ** (ordN - 1) - 0.5)
    Prx = round((px - oxy) * math.cos(angle) + (py - oxy) * math.sin(angle) + oxy)
    Pry = round((py - oxy) * math.cos(angle) - (px - oxy) * math.sin(angle) + oxy)
    Pr = [Prx, Pry]
    return Pr


# 3) Ploc to H的旋转
def nP2H_rotate(nP, ordN, ri):  # P是array
    Pr = []
    for j in range(len(nP)):
        Pr.append(xy_rotate(nP[j], ordN, ri))
    Hr = encode(np.array(Pr), 2, ordN)
    return Hr


####平移函数
# P(x,y)向右上方平移(2^i-1)个单位: s=0,1,3,7...，得到Ps(xs,ys)

def nP2H_shift(nP, ordN, si):
    if si < ordN:
        Ps = nP + 2 ** si - 1
        Hs = encode(np.array(Ps), 2, ordN + 1)
        return Hs


def gen_arrayH(ordN, dro, dso):
    '''
    键: H0
    值: H1, H2, ..., Hr
    :return: Hil字典的矩阵，列号即为H0值，节省一列的空间占用
    '''
    arrpath = './map_data/arrayH_' + str(ordN) + '_ro' + str(dro) + '_so' + str(dso) + '.npy'
    plcpath = './map_data/Ploc_' + str(ordN) + '.npy'
    if not os.path.exists(arrpath):
        t1 = time.time()
        H0, Ploc = gen_allp(ordN)
        Hrs = []
        for ri in range(1, dro + 1):  # 旋转0, 90, 180, 270: ro = 0,1,2,3
            H0_ri = nP2H_rotate(Ploc, ordN, ri)
            Hrs.append(H0_ri)
        for si in range(1, dso + 1):  # 向右上方平移2^so-1个单位: so = 0,1,2,...
            H0_si = nP2H_shift(Ploc, ordN, si)
            Hrs.append(H0_si)
        arrayH = np.transpose(Hrs).astype(int)
        # dicH = dict(zip(H0, Hrs))  # 转换为字典
        t2 = time.time()
        ta = t2-t1
        np.save('./map_data/arrayH_' + str(ordN) + '_ro' + str(dro) + '_so' + str(dso) + '.npy', arrayH)
        if not os.path.exists(plcpath):
            np.save('./map_data/Ploc_' + str(ordN) + '.npy', Ploc)
        sab = os.path.getsize('./map_data/arrayH_' + str(ordN) + '_ro' + str(dro) + '_so' + str(dso) + '.npy')  # 以字节为单位
        sam = round(sab / (1024 * 1024), 2)
        with open('./map_data/map_time_size.csv', 'a') as f:
            # f.write('映射字典耗时\n' + str(td) + '\n加密矩阵耗时\n' + str(te))
            f.write('映射字典耗时,ordN,' + str(ordN) + ',dro' + str(dro) + ',dso,' + str(dso) + ',time_s,' + str(ta) + ',size_mb,' + str(sam) + '\n')
        return arrayH, ta, Ploc


def load_arHrs_(ordN, dro: int = 3, dso: int = 3):
    '''
    根据全集矩阵（省去第1列），取出所需ro, so的几列矩阵
    :param arrayH: 全集矩阵_ordN_dro_dso
    :param ro: 旋转次数
    :param so: 平移次数
    :return: 该ro-so情形的部分映射表矩阵，补上第一列为H0:(0,2^2n-1)
    '''
    arrpath = './map_data/arrayH_' + str(ordN) + '_ro' + str(dro) + '_so' + str(dso) + '.npy'
    plcpath = './map_data/Ploc_' + str(ordN) + '.npy'
    if not os.path.exists(arrpath):
        arrayH, ta, Ploc = gen_arrayH(ordN, dro, dso)
    else:
        arrayH = np.load(arrpath, allow_pickle=True)
        Ploc = np.load(plcpath, allow_pickle=True)
    # arHrs = np.hstack([arrayH[:, :ro], arrayH[:, 3:(3+so)]])
    H0 = np.arange(2**(2*ordN)).reshape((2**(2*ordN),1))      # H0
    arHrs_ = np.hstack([H0, arrayH])
    return arHrs_, Ploc


def take_arHrs(arHrs_, ro, so):
    '''
    根据全集矩阵（包含第1列），取出所需ro, so的几列矩阵
    '''
    arHrs = np.hstack([arHrs_[:, :(ro+1)], arHrs_[:, 4:(4+so)]])
    return arHrs




if __name__ == "__main__":
    # ordN = 5
    # so = 1
    # context = gencontext()
    # dicH, EPloc, Ploc, td, te = gen_HP(ordN, so, context)
    # print('映射字典耗时', td, '\n加密矩阵耗时', te)

    # ordN1 = 6
    # ordN2 = 8
    # so1 = 3
    # so2 = 3
    # ro1 = 3
    # ro2 = 3
    # # En = True
    # En = False
    # td, te = save_all(ordN1, ordN2, so2, ro2, En)
    # tab_time(ordN1, ordN2, so1, so2, ro1, ro2, En)
    # context = load_context()
    # dicH, EPloc = load_cipher(ordN1, ro2, so2, context, En)

    ordN = 7
    dro = 3
    dso = 3

    gen_arrayH(ordN, dro, dso)

    # print(dicH)
    # print(Dec_ten(EPloc))

    # np.save('ctx_' + str(ordN) + '_' + str(so) + '.npy', context)
    # ## 生成映射表，存储序列化后的table、context
    # t1 = time.time()
    # tab = gen_table2E(ordN, so, context)
    # t2 = time.time()
    # print(t2-t1)
    # tab[0][0].decrypt()
    # tab[1][0].decrypt()
