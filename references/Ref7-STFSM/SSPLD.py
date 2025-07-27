#encoding=utf-8
#Algorithm 1: SSPLD Computation Protocol
import random
import gmpy2
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

from SM import *
from SDC import *
from SSED import *

# running on both C1 and C2
# C1 only has pk, C2 has both pk and sk
def SSPLD_PART1(sk: SK,E_p0x: EN,E_p0y: EN,E_p1x: EN,E_p1y: EN,E_p2x: EN,E_p2y: EN,E_c2: EN):
    E_a2=SSED(sk,E_p0x,E_p0y,E_p1x,E_p1y)
    E_b2=SSED(sk,E_p0x,E_p0y,E_p2x,E_p2y)
    E_s_1=SM(sk,E_p0x,E_p1y)
    E_s_2=SM(sk,E_p1x,E_p2y)
    E_s_3=SM(sk,E_p2x,E_p0y)
    E_s_p1=E_s_1+E_s_2+E_s_3
    E_s_4=SM(sk,E_p2x,E_p1y)
    E_s_5=SM(sk,E_p1x,E_p0y)
    E_s_6=SM(sk,E_p0x,E_p2y)
    E_s_p2=E_s_4+E_s_5+E_s_6
    E_s=E_s_p1-E_s_p2
    E_s2=SM(sk,E_s,E_s)
    E_h2=SDC(sk,E_s2,E_c2)
    return E_a2,E_b2,E_h2

# running on C1
def SSPLD_PART2(E_a2: EN,E_b2: EN,E_c2: EN):
    pk=E_a2.public_key
    n=pk.n
    if n>1024:
        rand_max=1024
    else:
        rand_max=n
    r1=random.randint(0,rand_max)
    r2=random.randint(0,rand_max)
    alpha=(E_b2-E_a2-E_c2)*r1
    beta=(E_a2-E_b2-E_c2)*r2
    return alpha,beta

# running on C2 with sk
def SSPLD_PART3(sk: SK,alpha: EN,beta: EN):
    pk=alpha.public_key
    D_alpha=sk.decrypt(alpha)
    D_beta=sk.decrypt(beta)
    if D_alpha>=0:
        theta_0=pk.encrypt(1)
        theta_1=pk.encrypt(0)
        theta_2=pk.encrypt(0)
    elif D_beta>=0:
        theta_0=pk.encrypt(0)
        theta_1=pk.encrypt(1)
        theta_2=pk.encrypt(0)
    else:
        theta_0=pk.encrypt(0)
        theta_1=pk.encrypt(0)
        theta_2=pk.encrypt(1)
    return theta_0,theta_1,theta_2

#running on both C1 and C2, actually on C1
def SSPLD_PART4(sk: SK,E_a2: EN,E_b2: EN,E_h2: EN,theta_0: EN,theta_1: EN,theta_2: EN)->EN:
    E_d_1=SM(sk,E_a2,theta_0)
    E_d_2=SM(sk,E_b2,theta_1)
    E_d_3=SM(sk,E_h2,theta_2)
    E_d=E_d_1+E_d_2+E_d_3
    return E_d

def SSPLD(sk: SK,E_p0x: EN,E_p0y: EN,E_p1x: EN,E_p1y: EN,E_p2x: EN,E_p2y: EN,E_c2: EN)->EN:
    E_a2,E_b2,E_h2=SSPLD_PART1(sk,E_p0x,E_p0y,E_p1x,E_p1y,E_p2x,E_p2y,E_c2)
    alpha,beta=SSPLD_PART2(E_a2,E_b2,E_c2)
    theta_0,theta_1,theta_2=SSPLD_PART3(sk,alpha,beta)
    E_result=SSPLD_PART4(sk,E_a2,E_b2,E_h2,theta_0,theta_1,theta_2)
    return E_result

def SSPLD_(sk: SK,E_p0x: EN,E_p0y: EN,E_p1x: EN,E_p1y: EN,E_p2x: EN,E_p2y: EN)->EN:
    E_c2=SSED(sk,E_p1x,E_p1y,E_p2x,E_p2y)
    E_a2,E_b2,E_h2=SSPLD_PART1(sk,E_p0x,E_p0y,E_p1x,E_p1y,E_p2x,E_p2y,E_c2)
    alpha,beta=SSPLD_PART2(E_a2,E_b2,E_c2)
    theta_0,theta_1,theta_2=SSPLD_PART3(sk,alpha,beta)
    E_result=SSPLD_PART4(sk,E_a2,E_b2,E_h2,theta_0,theta_1,theta_2)
    return E_result

if __name__=='__main__':
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    E_p0x=pk.encrypt(1)
    E_p0y=pk.encrypt(3)
    E_p1x=pk.encrypt(0)
    E_p1y=pk.encrypt(0)
    E_p2x=pk.encrypt(3)
    E_p2y=pk.encrypt(3)
    E_c2=SSED(sk,E_p1x,E_p1y,E_p2x,E_p2y)
    E_res=SSPLD(sk,E_p0x,E_p0y,E_p1x,E_p1y,E_p2x,E_p2y,E_c2)
    print(sk.decrypt(E_res))