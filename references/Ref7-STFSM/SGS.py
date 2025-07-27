import numpy as np
from typing import List

from convert_2D_to_1D import *

from phe.paillier import PaillierPublicKey as PK
from phe.paillier import EncryptedNumber as EN

def SGS_1point(center: List[float],radius: float)->List[List[int]]:
    x=center[0]
    y=center[1]
    r=radius
    result=[]
    for temp_x in range(int(x-r),int(x+r)+1):
        for temp_y in range(int(y-r),int(y+r)+1):
            if (temp_x-x)**2+(temp_y-y)**2<=r**2:
                result.append([temp_x,temp_y])
    return result

def SGS(line: List[List[float]],radius: float)->List[List[List[int]]]:
    SGS_result=[]
    for each in line:
        SGS_result.append(SGS_1point(each,radius))
    return SGS_result

def ESGS(pk: PK,line: List[List[float]],radius: float, max_coord: int):
    SGS_res=SGS(line,radius)
    ESGS_res=[]
    for each_circle in SGS_res:
        ESGS_circle=[]
        for each_point in each_circle:
            ESGS_circle.append(pk.encrypt(convert(each_point[0],each_point[1],max_coord)))
        ESGS_res.append(ESGS_circle)
    return ESGS_res


if __name__=='__main__':
    x=2.1
    y=2
    r=2.1
    res=SGS([x,y],r)
    print(res)