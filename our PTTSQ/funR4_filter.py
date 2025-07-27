# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 08:38:34 2022

@author: YiKelai
"""


# topk的正确率，ie,ih分别是前ke,kh个距离最近的轨迹下标
def topk_rate1(ie,ih):
    cr = set(ie).intersection(set(ih))     # 旋转r正确集合
    r1 = len(cr)/len(ie)
    return r1

