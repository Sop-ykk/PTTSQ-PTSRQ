#encoding=utf-8
#Algorithm 2: Secure Trajectory Filtering based on Signature Matching
import random
import gmpy2
import numpy as np
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

from SGS import *
from TGS import *
from convert_2D_to_1D import *

#running on c1 without sk
#E_TT needs to be encoded from (x,y) to y
def STFSM_PART1(E_TT: List[EN], E_SQ: List[List[EN]], E_TQ: List[EN], E_ST: List[List[EN]]):
    pk=E_TT[0].public_key
    n=pk.n
    if n>1024:
        rand_max=1024
    else:
        rand_max=n
    rand_1=[]
    rand_2=[]
    omega_1=[]
    omega_2=[]
    for i in range(0,len(E_TT)):
        rand_1.append(random.randint(0,rand_max))
    for i in range(0,len(E_TQ)):
        rand_2.append(random.randint(0,rand_max))
    for i in range(0,len(E_SQ)):
        omega_1_sub=[]
        for j in range(0,len(E_SQ[i])):
            omega_1_subsub=[]
            for k in range(0,len(E_TT)):
                temp=(E_SQ[i][j]-E_TT[k])*rand_1[k]
                omega_1_subsub.append(temp)
            omega_1_sub.append(omega_1_subsub)
        omega_1.append(omega_1_sub)
    for i in range(0,len(E_ST)):
        omega_2_sub=[]
        for j in range(0,len(E_ST[i])):
            omega_2_subsub=[]
            for k in range(0,len(E_TQ)):
                temp=(E_ST[i][j]-E_TQ[k])*rand_2[k]
                omega_2_subsub.append(temp)
            omega_2_sub.append(omega_2_subsub)
        omega_2.append(omega_2_sub)
    return omega_1,omega_2

#running on c2 with sk
def STFSM_PART2(sk: SK, omega_1: List[List[List[EN]]], omega_2: List[List[List[EN]]]):
    pk=omega_1[0][0][0].public_key
    e_1=[]
    for i in range(0,len(omega_1)):
        e_i=1
        for j in range(0,len(omega_1[i])):
            for k in range(0,len(omega_1[i][j])):
                if sk.decrypt(omega_1[i][j][k])==0:
                    e_i=0
                    break
            if e_i==0:
                break
        e_1.append(pk.encrypt(e_i))
    e_2=[]
    for i in range(0,len(omega_2)):
        e_i=1
        for j in range(0,len(omega_2[i])):
            for k in range(0,len(omega_2[i][j])):
                if sk.decrypt(omega_2[i][j][k])==0:
                    e_i=0
                    break
            if e_i==0:
                break
        e_2.append(pk.encrypt(e_i))
    return e_1,e_2

#running on c1 without sk
def STFSM_PART3(e_1: List[EN], e_2: List[EN])->EN:
    res_1=e_1[0]
    for i in range(1,len(e_1)):
        res_1=res_1+e_1[i]
    res_2=e_2[0]
    for i in range(1,len(e_1)):
        res_2=res_2+e_2[i]
    res=res_1+res_2
    return res

def STFSM(sk: SK, E_TT: List[EN], E_SQ: List[List[EN]], E_TQ: List[EN], E_ST: List[List[EN]])->EN:
    omega_1,omega_2=STFSM_PART1(E_TT,E_SQ,E_TQ,E_ST)
    e_1,e_2=STFSM_PART2(sk,omega_1,omega_2)
    res=STFSM_PART3(e_1,e_2)
    return res

if __name__=='__main__':
    T=[[1,1],[5,5],[5,10],[3,3]]
    Q=[[0,4],[8,4],[5,8],[5,3]]
    r=2
    max_coord=16
    SGS_T=[]
    SGS_Q=[]
    TGS_T=[]
    TGS_Q=[]
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    '''
    for i in range(0,len(T)):
        SGS_T_temp=SGS_1point(T[i],r)
        SGS_T_1D=[]
        for each in SGS_T_temp:
            SGS_T_1D.append(pk.encrypt(convert(each[0],each[1],max_coord)))
        SGS_T.append(SGS_T_1D)
        SGS_Q_temp=SGS_1point(Q[i],r)
        SGS_Q_1D=[]
        for each in SGS_Q_temp:
            SGS_Q_1D.append(pk.encrypt(convert(each[0],each[1],max_coord)))
        SGS_Q.append(SGS_Q_1D)
    for i in range(1,len(T)):
        TGS_T_temp=TGS_2point(T[i-1],T[i])
        for each in TGS_T_temp:
            TGS_T.append(pk.encrypt(convert(int(each[0]),int(each[1]),max_coord)))
        TGS_Q_temp=TGS_2point(Q[i-1],Q[i])
        for each in TGS_Q_temp:
            TGS_Q.append(pk.encrypt(convert(int(each[0]),int(each[1]),max_coord)))
    '''
    TGS_T=ETGS(pk,T,max_coord)
    TGS_Q=ETGS(pk,Q,max_coord)
    SGS_T=ESGS(pk,T,r,max_coord)
    SGS_Q=ESGS(pk,Q,r,max_coord)
    a=STFSM(sk,TGS_T,SGS_Q,TGS_Q,SGS_T)
    D_a=sk.decrypt(a)
    print(D_a)
