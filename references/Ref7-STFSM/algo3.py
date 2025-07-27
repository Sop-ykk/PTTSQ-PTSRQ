#encoding=utf-8
#Algorithm 3: Signature-Based Secure Trajectory Similarity Search
import random
import gmpy2
import numpy as np
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

from SM import *
from SGS import *
from TGS import *
from SMC import *
from SDC import *
from STFSM import *
from SSPLD import *

large_enough=1E20+9

#this part is used to generate ETGS and ESGS of Q and T
def SBSTSS_INIT(pk: PK, T: List[List[List[int]]], Q: List[List[int]], radius: int, max_coord: int):
    E_TQ=ETGS(pk,Q,max_coord)
    E_SQ=ESGS(pk,Q,radius,max_coord)
    E_TTs=[]
    E_STs=[]
    E_Ts=[]
    E_Q=[]
    v=[]
    for each_line in T:
        v_i=[]
        E_TTs.append(ETGS(pk,each_line,max_coord))
        E_STs.append(ESGS(pk,each_line,radius,max_coord))
        E_Ti=[]
        for each_point in each_line:
            E_Ti.append([pk.encrypt(each_point[0]),pk.encrypt(each_point[1])])
            v_i.append(pk.encrypt(1))
            #E_Ti.append(pk.encrypt(convert(each_point[0],each_point[1],max_coord)))
        E_Ts.append(E_Ti)
        v.append(v_i)
    for each_point in Q:
        E_Q.append([pk.encrypt(each_point[0]),pk.encrypt(each_point[1])])
        #E_Q.append(pk.encrypt(convert(each_point[0],each_point[1],max_coord)))
    return E_TQ,E_SQ,E_TTs,E_STs,E_Ts,E_Q,v

#running on C1 and C2
#line: 1-2
def SBSTSS_PART1(sk: SK,E_TTs,E_SQ,E_TQ,E_STs)->List[EN]:
    sigma=[]
    for i in range(0,len(E_TTs)):
        sigma.append(STFSM(sk,E_TTs[i],E_SQ,E_TQ,E_STs[i]))
    return sigma

#running on C1 without sk
#line: 3-4
def SBSTSS_PART2(sigma: List[EN]):
    eta=[]
    rand=[]
    pk=sigma[0].public_key
    n=pk.n
    if n<1024:
        rand_max=n
    else:
        rand_max=1024
    for i in range(0,len(sigma)):
        rand.append(random.randint(0,rand_max))
        eta.append(sigma[i]*rand[i])
    #require a random permutation function for eta here
    eta_=eta
    return eta_

#running on C2 with sk
#line: 6 - 10
def SBSTSS_PART3(sk: SK, eta_: List[EN]):
    pk=eta_[0].public_key
    L_=0
    M_=[]
    for i in range(0,len(eta_)):
        if sk.decrypt(eta_[i])==0:
            A=[]
            L_+=1
            for i in range(0,len(eta_)):
                A.append(pk.encrypt(0))
            A[i]=pk.encrypt(1)
            M_.append(A)
    M_=np.array(M_).T.tolist()
    return M_,L_

#running on C1
#line: 11
def SBSTSS_PART4(E_Ts, M_: List[List[EN]],L_: int):
    #require a reverse permutation function here
    M=M_
    #FILTER() function here,confusing
    E_Ts_=[]
    for i in range(0,L_):
        E_T_=[]
        for j in range(0,len(E_Ts)):
            temp_x=SM(sk,E_Ts[0][j][0],M[j][i])
            temp_y=SM(sk,E_Ts[0][j][1],M[j][i])
            for k in range(1,len(E_Ts)):
                temp_x+=SM(sk,E_Ts[k][j][0],M[k][i])
                temp_y+=SM(sk,E_Ts[k][j][1],M[k][i])
            E_T_.append([temp_x,temp_y])
        E_Ts_.append(E_T_)
    '''
    for i in range(0,len(E_Ts[0])):
        E_T_=[]
        for j in range(0,len(E_Ts)):
            E_T_.append(E_Ts[i][j]+M[j][i])
        E_Ts_.append(E_T_)
    E_Ts_=np.array(E_Ts_).T.tolist()
    '''
    return E_Ts_

#running on C1 and C2
#line: 12 - 24
def SBSTSS_PART5(sk: SK, L_: int, E_Ts_, E_Q,v):
    phi=[]
    delta=[]
    E_d2_TQ=[]
    E_d2_QT=[]
    E_T2Q=[]
    for i in range(0,L_):
        E_d2_TQ_i=[]
        E_d2_QT_i=[]
        for j in range(0,len(E_Ts_[i])):
            phi_j=pk.encrypt(large_enough)
            for k in range(0,len(E_Q)-1):
                delta_k=SSPLD_(sk,E_Ts_[i][j][0],E_Ts_[i][j][1],E_Q[k][0],E_Q[k][1],E_Q[k+1][0],E_Q[k+1][1])
                delta.append(delta_k)
                phi_j=SMC_(sk,phi_j,delta_k)
            phi.append(phi_j)
        for j in range(0,len(E_Ts_[i])):
            E_d2_TQ_i.append(SM(sk,phi[j],v[i][j]))
        for j in range(0,len(E_Q)):
            E_d2_QT_temp=pk.encrypt(large_enough)
            for k in range(0,len(E_Ts_[i])-1):
                delta_k=SSPLD_(sk,E_Q[j][0],E_Q[j][1],E_Ts_[i][k][0],E_Ts_[i][k][1],E_Ts_[i][k+1][0],E_Ts_[i][k+1][1])
                E_d2_QT_temp=SMC_(sk,E_d2_QT_temp,delta_k)
            E_d2_QT_i.append(E_d2_QT_temp)
        E_d2_TQ.append(E_d2_TQ_i)
        E_d2_QT.append(E_d2_QT_i)
    return E_d2_TQ,E_d2_QT
            

#running on C1
#line: 25 - 31
def SBSTSS_PART6(sk: SK,E_d2_TQ: List[List[EN]],E_d2_QT: List[List[EN]],E_Q,v,L_: int):
    pk=E_d2_TQ[0][0].public_key
    n=pk.n
    if n<1024:
        rand_max=n
    else:
        rand_max=1024
    Phi=[]
    for i in range(0,L_):
        Phi.append(random.randint(0,rand_max))
    Gamma=[]
    Lambda=[]
    E_SIM_=[]
    for i in range(0,L_):
        gamma_temp=E_d2_TQ[i][0]
        for j in range(1,len(E_d2_TQ[i])):
            gamma_temp+=E_d2_TQ[i][j]
        for j in range(0,len(E_d2_QT[i])):
            gamma_temp+=E_d2_QT[i][j]
        Gamma.append(gamma_temp)
        lambda_temp=v[0][i]+pk.encrypt(len(E_Q))
        Lambda.append(lambda_temp)
        E_SIM_temp=SDC(sk,gamma_temp,lambda_temp)+pk.encrypt(Phi[i])
        E_SIM_.append(E_SIM_temp)
    return E_SIM_,Phi

#running on C2 with sk
#line: 32
def SBSTSS_PART7(sk: SK,E_SIM_: List[EN])-> List[int]:
    SIM_=[]
    for each in E_SIM_:
        SIM_.append(sk.decrypt(each))
    return SIM_

#running on U
#line: 33
def SBSTSS_PART8(SIM_: List[int],Phi: List[int]):
    res=[]
    for i in range(0,len(SIM_)):
        res.append(SIM_[i]-Phi[i])
    return res

def SBSTSS(sk: SK,pk: PK, T: List[List[List[int]]], Q: List[List[int]], radius: int, max_coord: int)-> List[int]:
    E_TQ,E_SQ,E_TTs,E_STs,E_Ts,E_Q,v=SBSTSS_INIT(pk,T,Q,radius,max_coord)
    sigma=SBSTSS_PART1(sk,E_TTs,E_SQ,E_TQ,E_STs)
    eta_=SBSTSS_PART2(sigma)
    M_,L_=SBSTSS_PART3(sk,eta_)
    E_Ts_=SBSTSS_PART4(E_Ts,M_,L_)
    E_d2_TQ,E_d2_QT=SBSTSS_PART5(sk,L_,E_Ts_,E_Q,v)
    E_SIM_,Phi=SBSTSS_PART6(sk,E_d2_TQ,E_d2_QT,E_Q,v,L_)
    SIM_=SBSTSS_PART7(sk,E_SIM_)
    res=SBSTSS_PART8(SIM_,Phi)
    return res

if __name__=='__main__':
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    Q=[[1,1],[5,3],[3,6],[5,8]]
    T1=[[3,1],[3,4],[5,6],[4,8]]
    T2=[[8,1],[10,3],[9,6],[11,4]]
    T3=[[5,1],[2,2],[6,4],[3,7]]
    T=[T1,T2,T3]
    r=3
    max_coord=16
    res=SBSTSS(sk,pk,T,Q,r,max_coord)
    print(res)





