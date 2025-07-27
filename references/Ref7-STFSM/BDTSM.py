#encoding=utf-8
#Definition 1: Bi-Directional Trajectory Similarity Measure
from math import inf
from typing import List


import numpy as np
import numpy.linalg as LA

# T_len() is used to calculate the length of trajectory T
# not used
def T_len(T: List[List[float]]) ->float:
    length=0
    for i in range(0,len(T)-1):
        length+=np.sqrt((T[i][0]-T[i+1][0])**2+(T[i][1]-T[i+1][1])**2)
    return length

# Dist_ppp() is used to calculate the minimum distance between p and p1->p2
def Dist_ppp(p: List[float],p1: List[float],p2: List[float])->float:
    a_sq=(p[0]-p1[0])**2+(p[1]-p1[1])**2
    b_sq=(p[0]-p2[0])**2+(p[1]-p2[1])**2
    c_sq=(p1[0]-p2[0])**2+(p1[1]-p2[1])**2
    a=np.sqrt(a_sq)
    b=np.sqrt(b_sq)
    if b_sq>=a_sq+c_sq:
        return a
    elif a_sq>=b_sq+c_sq:
        return b
    else:    
        #S=np.array([[p[0],p[1],1],[p1[0],p1[1],1],[p2[0],p2[1],1]])
        #h=LA.det(S)/np.sqrt(c_sq)
        c=np.sqrt(c_sq)
        s=(a+b+c)/2
        A=np.sqrt(s*(s-a)*(s-b)*(s-c))
        h=2*A/c
        return h

# Dist_PT() is used to calculate the minimum distance between point p and trajectory T
# in paper it's Dist_{PT}(p^k_{T_i},T_j)
def Dist_PT(p: List[float],T: List[List[float]]) ->float:
    result=[]
    for i in range(0,len(T)-1):
        result.append(Dist_ppp(p,T[i],T[i+1]))
    return min(result)

# apply D_max threshold
# in paper it's d^k_{T_i \rightarrow T_j}
def d_TiTj(p: List[float],T: List[List[float]],D_max: float) ->float:
    real_dist=Dist_PT(p,T)
    if real_dist > D_max:
        return inf
    else:
        return real_dist

# Bi-Directional Trajectory Similaroty Measure
# in paper it's SIM(T_i,T_j)
def SIM(T1: List[List[float]], T2: List[List[float]],D_max: float) ->float:
    sum1=0
    sum2=0
    for each in T1:
        sum1+=d_TiTj(each,T2,D_max)
    for each in T2:
        sum2+=d_TiTj(each,T1,D_max)
    div1=len(T1)
    div2=len(T2)
    result=(sum1+sum2)/(div1+div2)
    return result