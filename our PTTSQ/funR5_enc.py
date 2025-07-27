# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:16:16 2022

SHE加密

@author: YiKelai
"""

#%%
#### SHE加密

import gmpy2
import random
import time
import libnum

### 1) Key Generation

def get_prime(rs,bitsize):
    p = gmpy2.mpz_urandomb(rs, bitsize)
    while not gmpy2.is_prime(p):    # 直到为素数
        p = p + 1
    return p

def keygen(k0,k1,k2):
    rs = gmpy2.random_state(int(time.time()))
    p = get_prime(rs,k0)
    q = get_prime(rs,k0)
    #print(p, q)
    n = p * q
    L = gmpy2.mpz_urandomb(rs, k2)

    pp = [n,k0,k1,k2]
    sk = [p,L]
    return pp, sk

### 2.1) sk加密
# 𝐸𝑛𝑐(𝑚): 𝑐 = (𝑟𝐿+𝑚)(1+𝑟'𝑝) 𝑚𝑜𝑑 n
# plaintext 介于 (-2**k1,2**k1)
def enc(plaintext, pp, sk):
    rs = gmpy2.random_state(int(time.time()))
    m = gmpy2.mpz(int(plaintext))
    [n,k0,k1,k2] = pp
    [p,L] = sk
    r = gmpy2.mpz_urandomb(rs, k2)
    r1 = gmpy2.mpz_urandomb(rs, k0)
    # c = ((r*L+m) % n)*((1+r1*p) % n) % n
    a = gmpy2.mod(gmpy2.add(gmpy2.mul(r, L), m),n)
    b = gmpy2.mod(gmpy2.add(gmpy2.mul(r1, p), 1),n)
    c = gmpy2.mod(gmpy2.mul(a,b),n)
    return c

### 2.2) pk加密
# 𝐸𝑛𝑐(𝑚): 𝑐 = 𝑚 + 𝑟1·𝐸1(0) + 𝑟2·𝐸2(0) 𝑚𝑜𝑑 𝑁

def encp(plaintext, pk):
    rs = gmpy2.random_state(int(time.time()))
    m = gmpy2.mpz(int(plaintext))
    [n,k0,k1,k2,E01,E02] = pk
    
    r1 = gmpy2.mpz_urandomb(rs, k2)
    r2 = gmpy2.mpz_urandomb(rs, k2)

    a = gmpy2.mod(gmpy2.mul(r1,E01),n)
    b = gmpy2.mod(gmpy2.mul(r2,E02),n)
    c = gmpy2.mod(gmpy2.add(m,gmpy2.add(a,b)),n)
    return c

# 3) sk解密
# 𝐷𝑒𝑐(𝑐): 𝑚 = (𝑐 𝑚𝑜𝑑 𝑝) 𝑚𝑜𝑑 𝐿.
def dec(c, pp, sk):
    [n,k0,k1,k2] = pp
    [p,L] = sk
    m1 = (c % p) % L
    if m1 <= (L/2):
        m = m1
    else:
        m = m1-L
    return m

# 4.1) 同态加法
def hadd(c1, c2, pp):
    n = pp[0]
    ha1 = gmpy2.mod(gmpy2.add(c1, c2), n)
    return ha1

# 4.2) 同态乘法
def hmul(c1, c2, pp):
    n = pp[0]
    hm1 = gmpy2.mod(gmpy2.mul(c1,c2), n)
    return hm1

# 4.3) 同态减法: c1-c2
# En1 = encp(-1, pk)      # E(-1)

def hsub(c1, c2, pp, En1):
    n = pp[0]
    nc2 = gmpy2.mod(gmpy2.mul(En1, c2), n)       # E(-c2)
    hs1 = gmpy2.mod(gmpy2.add(c1,nc2), n)
    # nc2 = hmul(En1, c2, pk[:4])       # E(-c2)
    # hs1 = hadd(c1, nc2, pk[:4])
    return hs1





#%%
## SHE加密Ploc数组(2^2n,2)的每个元素
import numpy as np

def enc_Ploc(Ploc, pk):
    EP = []
    for i in range(len(Ploc)):
        EPx = encp(Ploc[i][0], pk)
        EPy = encp(Ploc[i][1], pk)
        EP.append([EPx, EPy])
    EP = np.asarray(EP)
    return EP
        
#%%

# k0=2048
# k1=24
# k2=160

# [pp, sk] = keygen(k0,k1,k2)
# E01 = enc(0, pp, sk)
# E02 = enc(0, pp, sk)
# pk = pp + [E01, E02]

# m = 2**22
# Em = encp(m, pk)
# DEm = dec(Em, pp, sk)

# order = 3
# H0,Ploc = gen_allp(order)
# EPloc = enc_Ploc(Ploc, pk)

#%% --------------------------- Script -----------------------------




