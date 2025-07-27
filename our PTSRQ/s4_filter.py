#### 除了topk查询以外，增加范围查询的功能
from dH_pivot import DSED, DSED_hil, DSED_hilpv
import heapq
import numpy as np


def nDSED(dist_type:str, p1_qu, pn_do, t_qu):
    '''
    查询与多条轨迹之间的距离
    :param p1_qu: (点数,2)
    :param t_qu: (点数,)
    :param pn_do: [条数,点数]
    :param nt_do: [条数,点数]
    :return: 列表 [d1, d2, ..., d_un]
    '''
    distn = []
    for ui in range(len(pn_do)):
        disti = DSED(dist_type, p1_qu, pn_do[ui], t_qu)  # 查询Q与第ui条轨迹数据的距离
        distn.append(disti)
    return distn


def nDSED_Hil(h1_qu, hn_do, t_qu, arHrs):
    '''
    查询与多条轨迹之间的距离
    :param h1_qu: (点数,)
    :param t_qu: (点数,)
    :param hn_do: [条数,点数]
    :param nt_do: [条数,点数]
    :return: 列表 [d1, d2, ..., d_un]
    '''
    distn = []
    for ui in range(len(hn_do)):
        disti = DSED_hil(h1_qu, hn_do[ui], t_qu, arHrs)  # 查询Q与第ui条轨迹数据的距离
        distn.append(disti)
    return distn


def nDSED_HilPv(h1_qu, hn_do, t_qu, arHrs, radius, bi:bool=False):
    '''
    查询与多条轨迹之间的距离
    :param h1_qu: (点数,)
    :param t_qu: (点数,)
    :param hn_do: [条数,点数]
    :param nt_do: [条数,点数]
    :return: 列表 [d1, d2, ..., d_un]
    '''
    distn = []
    for ui in range(len(hn_do)):
        disti = DSED_hilpv(h1_qu, hn_do[ui], t_qu, arHrs, radius, bi)  # 查询Q与第ui条轨迹数据的距离
        distn.append(disti)
    return distn


## Top-k query
def idxh(kh, h1_qu, hn_do, t_qu, arHrs, radius, bi:bool=False):
    # 返回top_kh的索引下标
    dnh1 = nDSED_Hil(h1_qu, hn_do, t_qu, arHrs)
    ih1 = heapq.nsmallest(kh, range(len(dnh1)), dnh1.__getitem__)  # 获取前kh个H1距离最近的id
    dnh2 = nDSED_HilPv(h1_qu, hn_do, t_qu, arHrs, radius, bi)
    ih2 = heapq.nsmallest(kh, range(len(dnh2)), dnh2.__getitem__)  # 获取前kh个H1_pv距离最近的id
    return ih1, ih2


def idKh(kh, h1_qu, hn_do, t_qu, arHrs, radius, bi:bool=False):
    # 返回top_kh的索引下标
    dnh2 = nDSED_HilPv(h1_qu, hn_do, t_qu, arHrs, radius, bi)
    ih2 = heapq.nsmallest(kh, range(len(dnh2)), dnh2.__getitem__)  # 获取前kh个H1距离最近的id
    return ih2


def ids(dist_type, ke, p1_qu, pn_do, t_qu):
    # for verify
    dne = nDSED(dist_type, p1_qu, pn_do, t_qu)  # 欧式(平方)距离 for verify
    ie = heapq.nsmallest(ke, range(len(dne)), dne.__getitem__)  # 获取前ke个Eu距离最近的id
    return ie


# topk的正确率，ie,ih分别是前ke,kh个距离最近的轨迹下标
def topk_rate(ie, ih):
    cr = set(ie).intersection(set(ih))     # 旋转r正确集合
    if len(ie) == 0:
        r1 = 0
        print('no results  --ie')
    else:
        r1 = len(cr)/len(ie)
    return r1


def prec(dist_type, ke, kh, p1_qu, pn_do, t1_qu, h1_qu, hn_do, arHrs, radius, bi:bool=False):
    ih1, ih2 = idxh(kh, h1_qu, hn_do, t1_qu, arHrs, radius, bi)
    ie = ids(dist_type, ke, p1_qu, pn_do, t1_qu)
    pr1 = topk_rate(ie, ih1)
    pr2 = topk_rate(ie, ih2)
    return pr1, pr2


## Range query
def idRh(per2, h1_qu, hn_do, t1_qu, arHrs, radius, bi):
    dnh = nDSED_HilPv(h1_qu, hn_do, t1_qu, arHrs, radius, bi)
    epsh = np.quantile(dnh, per2)   # 根据假设per2, 用于HilPv初筛,(不计时
    irh2 = np.argwhere(dnh < epsh)     # 获取H2距离满足范围条件的id
    ir_Pv = list(map(int, irh2))     # 变为list类型,(不计时
    return ir_Pv


def ids_range(dist_type, eps0, p1_qu, pn_do, t1_qu, rk:bool = False, per = 0.25):
    '''

    :param eps0: 范围查询的epsilon阈值: = np.float64(100.0)
    :param rk: 是否按百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :param per: 0到1之间的小数，百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :return: 满足该eps范围或rk_per范围的轨迹下标
    '''
    # for verify
    dne = nDSED(dist_type, p1_qu, pn_do, t1_qu)  # 欧式(平方)距离 for verify
    if rk:
        eps0 = np.quantile(dne, per)
    ir = np.argwhere(dne < eps0)
    ir = list(map(int, ir))  # 变为list类型
    print(len(ir))
    return ir


def idh_range(eps0, err, h1_qu, hn_do, t1_qu, arHrs, radius, bi:bool=False, rk:bool = False, per = 0.25):
    '''

    :param eps0: 范围查询的epsilon阈值
    :param err: 用Hilbert曲线过滤时，对应的阈值选取为eps0的(1+err)倍
    :param rk: 是否按百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :param per: 0到1之间的小数，百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :return: 满足该eps范围或rk_per范围的轨迹下标
    '''
    dnh1 = nDSED_Hil(h1_qu, hn_do, t1_qu, arHrs)
    dnh2 = nDSED_HilPv(h1_qu, hn_do, t1_qu, arHrs, radius, bi)
    if rk:
        epsh1 = np.quantile(dnh1, per)
        epsh2 = np.quantile(dnh2, per)
        print(epsh1/eps0, epsh2/eps0)
    else:
        epsh1 = eps0 * (1 + err)    # 误差err>=0
        epsh2 = epsh1
    irh1 = np.argwhere(dnh1 < epsh1)     # 获取H1距离满足范围条件的id
    irh1 = list(map(int, irh1))     # 变为list类型
    irh2 = np.argwhere(dnh2 < epsh2)     # 获取H1距离满足范围条件的id
    irh2 = list(map(int, irh2))     # 变为list类型
    print('len(irh1), len(irh2)', len(irh1), len(irh2))
    return irh1, irh2


def idh_range2(dist_type, err, p1_qu, pn_do, h1_qu, hn_do, t1_qu, arHrs, radius, bi:bool=False, per = 0.25):
    '''
    :param eps0: 范围查询的epsilon阈值
    :param err: 用Hilbert曲线过滤时，对应的阈值选取为eps0的(1+err)倍
    :param rk: 是否按百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :param per: 0到1之间的小数，百分比位序（分位数），即，阈值eps0定为dne排序在第per的距离值
    :return: 满足该eps范围或rk_per范围的轨迹下标
    '''
    dne = nDSED(dist_type, p1_qu, pn_do, t1_qu)  # 欧式(平方)距离 for verify
    dnh1 = nDSED_Hil(h1_qu, hn_do, t1_qu, arHrs)
    dnh2 = nDSED_HilPv(h1_qu, hn_do, t1_qu, arHrs, radius, bi)

    eps0 = np.quantile(dne, per)
    epsh1 = eps0 * (1 + err)  # 误差err>=0
    epsh2 = epsh1

    irh1 = np.argwhere(dnh1 < epsh1)     # 获取H1距离满足范围条件的id
    irh1 = list(map(int, irh1))     # 变为list类型
    irh2 = np.argwhere(dnh2 < epsh2)     # 获取H1距离满足范围条件的id
    irh2 = list(map(int, irh2))     # 变为list类型
    print('len(irh1), len(irh2)', len(irh1), len(irh2))
    return irh1, irh2


# range_search的正确率，ie,ih分别是前ke,kh个距离最近的轨迹下标

# def range_prec(ie, ih):
## 与topk_rate函数相同
#     cr = set(ie).intersection(set(ih))     # 旋转r正确集合
#     r1 = len(cr)/len(ie)
#     return r1


def range_recall(ie, ih):
    cr = set(ie).intersection(set(ih))     # 旋转r正确集合
    if len(ih) == 0:
        r2 = 0
        print('no results -- ih')
    else:
        r2 = len(cr)/len(ih)
    return r2


def prrc(dist_type, eps0, err, p1_qu, pn_do, t1_qu, h1_qu, hn_do, arHrs, radius, bi:bool=False, rk:bool = False, per = 0.25):
    irh1, irh2 = idh_range(eps0, err, h1_qu, hn_do, t1_qu, arHrs, radius, bi, rk, per)
    ir = ids_range(dist_type, eps0, p1_qu, pn_do, t1_qu, rk, per)
    pr1 = topk_rate(ir, irh1)
    pr2 = topk_rate(ir, irh2)
    rc1 = range_recall(ir, irh1)
    rc2 = range_recall(ir, irh2)
    return pr1, pr2, rc1, rc2


def prrc2(dist_type, eps0, err, p1_qu, pn_do, t1_qu, h1_qu, hn_do, arHrs, radius, bi:bool=False, rk:bool = False, per = 0.25):
    irh1, irh2 = idh_range2(dist_type, err, p1_qu, pn_do, h1_qu, hn_do, t1_qu, arHrs, radius, bi, per)
    ir = ids_range(dist_type, eps0, p1_qu, pn_do, t1_qu, rk, per)
    pr1 = topk_rate(ir, irh1)
    pr2 = topk_rate(ir, irh2)
    rc1 = range_recall(ir, irh1)
    rc2 = range_recall(ir, irh2)
    return pr1, pr2, rc1, rc2

# eps0 = np.array(1000)
# err=0
# rk=True
# per=0.25
# bi=False
# pr1, pr2, rc1, rc2 = prec2('dMh', eps0, 0, p1_qu, pn_do, t1_qu, h1_qu, hn_do, arHrs, radius, rk)