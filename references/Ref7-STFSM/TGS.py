from typing import List
import numpy as np

from convert_2D_to_1D import *

from phe.paillier import PaillierPublicKey as PK
from phe.paillier import EncryptedNumber as EN

def seg_line_set(a,b):
    dx = int(b[0]-a[0])
    dy = int(b[1]-a[1])
    adx = abs(dx)
    ady = abs(dy)
    l = max(adx,ady)
    if adx >= ady:
        ddx = np.ones(l)*np.sign(dx)
        ddy = np.zeros(l)
        d0 = round((l-ady)/(ady+1))
        for i in range(int(l/(d0+1))):
            ddy[(i+1)*(d0+1)-1] = np.sign(dy)
    else:
        ddy = np.ones(l)*np.sign(dy)
        ddx = np.zeros(l)
        d0 = round((l-adx)/(adx+1))
        for i in range(int(l/(d0+1))):
            ddx[(i+1)*(d0+1)-1] = np.sign(dx)
    sx,sy = np.zeros(l+1),  np.zeros(l+1)
    sx[0] = a[0]
    sy[0] = a[1]
    for i in range(l-1):
        sx[i+1] = sx[i] + ddx[i]
        sy[i+1] = sy[i] + ddy[i]
    sx[l] = b[0]
    sy[l] = b[1]
    st = np.transpose([sx, sy]).astype(int)
    return st

def TGS_2point(a: List[float],b: List[float])->List[List[int]]:
    return seg_line_set(a,b)

def TGS(line: List[List[float]]):
    TGS_result=[]
    for i in range(1,len(line)):
        TGS_temp=TGS_2point(line[i-1],line[i])
        for each in TGS_temp:
            TGS_result.append([int(each[0]),int(each[1])])
    return TGS_result

def ETGS(pk: PK,line: List[List[float]], max_coord: int)->List[EN]:
    TGS_res=TGS(line)
    ETGS_res=[]
    for each in TGS_res:
        ETGS_res.append(pk.encrypt(convert(each[0],each[1],max_coord)))
    return ETGS_res

if __name__=='__main__':
    a=seg_line_set([0,0],[1.9,5.9])
    print(a)