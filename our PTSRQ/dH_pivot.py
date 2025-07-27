# 假设：
# 已经【等长】的DO、QU轨迹: 一条p_qu, 一条p_do
# 为避免重复生成table，假设统一固定为全集 ro = 3, so = 3，计算时再摘取实际所需的列(r<=7)


import numpy as np
from s0_data_pre import file_data
from f2_tr2H import traj_enc
from s1_table import load_arHrs_, take_arHrs


'''
-- dMh -- 曼哈顿距离
-- dEu -- 欧氏距离
-- sEu -- 欧式平方距离
-- dH  -- Hil距离
-- dHP -- 引入支点的Hil距离

-- CKKS -- 密文操作，另写
'''
def dMh(v1, v2):
    '''
    :param v1: 已对齐轨迹向量 v = [(x1,y1), (x2,y2), ...]: (len, 2)
    :return: 已对齐轨迹【欧氏距离】的向量 dv = [d1,d2,...]
    '''
    vd1 = np.abs(v1 - v2)
    dv = vd1[:, 0] + vd1[:, 1]
    # dv = np.sum(vd1)
    return dv


def dEu(v1, v2):
    '''
    :param v1: 已对齐轨迹向量 v = [(x1,y1), (x2,y2), ...]: (len, 2)
    :return: 已对齐轨迹【欧氏距离】的向量 dv = [d1,d2,...]
    '''
    vd2 = np.square(v1 - v2)
    dv = np.sqrt(vd2[:, 0] + vd2[:, 1])
    return dv


def sEu(v1, v2):
    '''
    :param v1: 已对齐轨迹向量 v = [(x1,y1), (x2,y2), ...]: (len, 2)
    :return: 已对齐轨迹【平方欧氏距离】的向量 dv = [d1,d2,...]
    '''
    vd2 = np.square(v1 - v2)
    dv2 = vd2[:, 0] + vd2[:, 1]
    return dv2



def dH(h1, h2, arHrs):
    '''
    for comparison without transformation
    :param h1: 已对齐轨迹向量 hv = [h1, h2, ...]: (len, 1)
    :param arHrs: 映射表矩阵
    :return: 已对齐轨迹【Hil距离】的向量 dv = [d1,d2,...]
    '''
    # if np.shape(arHrs)[1] == 1:
    #     # dh = abs(h1 - h2)       # origin hil dist
    #     hm1 = arHrs[h1]
    #     hm2 = arHrs[h2]
    #     dh = abs(hm1 - hm2)
    # else:
    hm1 = arHrs[h1]
    hm2 = arHrs[h2]
    dh = abs(hm1 - hm2).min(1)     # modified hil dist
    # dh = np.minimum(d0, dm)
    return dh


# version2: 若支点出界，则截断（clip）
def dHP(h1, h2, arHrs, radius, bi:bool=False):
    '''
    for comparison without transformation
    :param h1: 已对齐轨迹向量 hv = [h1, h2, ...]: (len, 1)
    :param arHrs: 映射表矩阵
    :return: 已对齐轨迹【Hil-pivot距离】的向量 dv = [d1,d2,...]
    '''
    # hm1 = np.array([arHrs[(i-radius):(i+radius+1)] for i in h1])
    pv = np.array([i-radius for i in range(2*radius+1)]).reshape(-1, 1)
    hm1 = h1 + pv
    hm1_clip = np.clip(hm1, 0, len(arHrs)-1)
    dhp2 = np.array([dH(hm1_clip[i,:], h2, arHrs) for i in range(2*radius+1)])           # 支点到h2的距离矩阵
    apv = abs(pv)      # 支点到h1的距离矩阵
    dh12 = dhp2 + apv
    dh = np.amin(dh12, axis=0)
    # print(np.argmin(dh12, axis=0))
    if bi:
        hm2 = h2 + pv     # 将h2每个元素上下扩充radius，得到2r+1列，最中间h1[:,r]即本身
        hm2_clip = np.clip(hm2, 0, len(arHrs) - 1)
        dhp1 = np.array([dH(h1, hm2_clip[i,:], arHrs) for i in range(2 * radius + 1)])          # 支点到h1的距离矩阵
        dh21 = dhp1 + apv
        dh_ = np.amin(dh21, axis=0)
        dh = [min(x) for x in zip(dh, dh_)]
    return dh


# version3: 考虑多个中间支点, 2023.11.13
def dHPm(h1, h2, arHrs, radius, bi:bool=False):
    '''
    for comparison without transformation
    :param h1: 已对齐轨迹向量 hv = [h1, h2, ...]: (len, 1)
    :param arHrs: 映射表矩阵
    :return: 已对齐轨迹【Hil-pivot距离】的向量 dv = [d1,d2,...]
    '''
    # # hm1 = np.array([arHrs[(i-radius):(i+radius+1)] for i in h1])
    # pv = np.array([i-radius for i in range(2*radius+1)]).reshape(-1, 1)
    # hm1 = h1 + pv
    # hm1_clip = np.clip(hm1, 0, len(arHrs)-1)
    # dhp2 = np.array([dH(hm1_clip[i,:], h2, arHrs) for i in range(2*radius+1)])           # 支点到h2的距离矩阵
    # apv = abs(pv)      # 支点到h1的距离矩阵
    # dh12 = dhp2 + apv
    # dh = np.amin(dh12, axis=0)
    # # print(np.argmin(dh12, axis=0))
    # if bi:
    #     hm2 = h2 + pv     # 将h2每个元素上下扩充radius，得到2r+1列，最中间h1[:,r]即本身
    #     hm2_clip = np.clip(hm2, 0, len(arHrs) - 1)
    #     dhp1 = np.array([dH(h1, hm2_clip[i,:], arHrs) for i in range(2 * radius + 1)])          # 支点到h1的距离矩阵
    #     dh21 = dhp1 + apv
    #     dh_ = np.amin(dh21, axis=0)
    #     dh = [min(x) for x in zip(dh, dh_)]
    return dh


def switch_dist(dist_type: str, v1, v2):
    if dist_type == 'dEu':
        return dEu(v1, v2)
    if dist_type == 'sEu':
        return sEu(v1, v2)
    if dist_type == 'dMh':
        return dMh(v1, v2)


def DSED(dist_type:str, p1_qu, p1_do, t_qu):
    # tid = np.nonzero(np.in1d(t_do, t_qu))[0]
    # p_do = p_do_[tid]
    s = switch_dist(dist_type, p1_qu, p1_do)
    h = len(t_qu) - 1
    sd2 = (t_qu[1] - t_qu[0]) * s[0] + (t_qu[h] - t_qu[h-1]) * s[h]
    for k in range(1, h):
        # print(k)
        sd2 += (t_qu[k + 1] - t_qu[k - 1]) * s[k]
    dsed = sd2 / (2 * (t_qu[h] - t_qu[0]))
    return dsed


def DSED_hil(h1_qu, h1_do, t_qu, arHrs):
    # tid = np.nonzero(np.in1d(t_do, t_qu))[0]
    # h_do = h_do_[tid]
    s = dH(h1_qu, h1_do, arHrs)
    h = len(t_qu) - 1
    sd2 = (t_qu[1] - t_qu[0]) * s[0] + (t_qu[h] - t_qu[h - 1]) * s[h]
    for k in range(1, h):
        # print(k)
        sd2 += (t_qu[k + 1] - t_qu[k - 1]) * s[k]
    dsed = sd2 / (2 * (t_qu[h] - t_qu[0]))
    return dsed


def DSED_hilpv(h1_qu, h1_do, t_qu, arHrs, radius, bi:bool=False):
    # tid = np.nonzero(np.in1d(t_do, t_qu))[0]
    # h_do = h_do_[tid]
    s = dHP(h1_qu, h1_do, arHrs, radius, bi)
    h = len(t_qu) - 1
    sd2 = (t_qu[1] - t_qu[0]) * s[0] + (t_qu[h] - t_qu[h - 1]) * s[h]
    for k in range(1, h):
        # print(k)
        sd2 += (t_qu[k + 1] - t_qu[k - 1]) * s[k]
    dsed = sd2 / (2 * (t_qu[h] - t_qu[0]))
    return dsed


if __name__ == '__main__':
    from s1_table import load_context, load_cipher
    import os
    import random
    ordN = 8
    dro, dso = 3, 3
    ro, so = 3, 1
    radius = 2

    # ctx = load_context()
    # dicH, EPloc = load_cipher(ordN, dro, dso, ctx, False)
    # arrayH = np.load('./map_data/arrayH_' + str(ordN) + '_ro' + str(dro) + '_so' + str(dso) + '.npy', allow_pickle=True)
    # os.chdir('./extension4j')
    arHrs_, Ploc = load_arHrs_(ordN, dro, dso)
    arHrs = take_arHrs(arHrs_, ro, so)
    # os.chdir('..')
    data = file_data('GT8_tall600')
    # data = file_data('GT11_tsyn600')
    # dist_type = "dEu"
    ri = random.randint(0, len(data)-1)
    rj = random.randint(0, len(data)-1)
    p1_do = data[0:, :15:3, 0:2][ri]
    # t_do = data[0:, 0:, 2][0]
    p1_qu = data[0:, :20:4, 0:2][rj]
    t_qu = data[0:, 0:20:4, 2][rj]

    h1_do = traj_enc(p1_do, ordN).astype(int)
    h1_qu = traj_enc(p1_qu, ordN).astype(int)

    dh1 = dH(h1_do, h1_qu, arHrs)
    dh2 = dHP(h1_do, h1_qu, arHrs, radius)
    mh = dMh(p1_do, p1_qu)
    eu = dEu(p1_do, p1_qu)

    print(dh1,dh2,mh,eu)

    DSED('dEu', p1_qu, p1_do, t_qu)
    DSED_hil(h1_qu, h1_do, t_qu, arHrs)
    DSED_hilpv(h1_qu, h1_do, t_qu, arHrs, radius)