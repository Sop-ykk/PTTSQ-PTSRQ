#encoding=utf-8
#Secure Division Computation
import random
from unittest import result
import gmpy2
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

# running on P1 with pk
def SDC_PART1(E_x: EN,E_d: EN):
    pk=E_x.public_key
    n=pk.n
    if pk.n>1024:
        rand_max=1024
    else:
        rand_max=pk.n
    r1=random.randint(0,rand_max)
    r2=random.randint(0,rand_max)
    x_=(E_x*r1)+(E_d*(r1*r2))
    d_=E_d*r1
    E_r2=pk.encrypt(r2)
    # x_1=gmpy2.powmod(gmpy2.mpz(E_x.ciphertext()),r1,n**2)
    # x_2=gmpy2.powmod(gmpy2.mpz(E_d.ciphertext()),r1*r2,n**2)
    # x_=gmpy2.mul(x_1,x_2)
    # #x_=gmpy2.t_mod(x_,n**4)
    # x_=EN(pk,int(x_))
    # d_=gmpy2.powmod(gmpy2.mpz(E_d.ciphertext()),r1,n**2)
    # d_=EN(pk,int(d_))
    # E_r2=pk.encrypt(r2)
    return x_,d_,E_r2

# running on P2 with pk and sk
def SDC_PART2(sk: SK,x_: EN,d_: EN)->EN:
    pk=x_.public_key
    D_x=sk.decrypt(x_)
    D_d=sk.decrypt(d_)
    h=int(D_x/D_d)
    return pk.encrypt(h)

# running on P1 with pk
def SDC_PART3(E_h: EN,E_r2: EN)->EN:
    pk=E_r2.public_key
    n=pk.n
    #result=E_h+(E_r2*(n-1))
    result=gmpy2.mul(gmpy2.mpz(E_h.ciphertext()),gmpy2.powmod(gmpy2.mpz(E_r2.ciphertext()),n-1,n**2))
    #result=gmpy2.t_mod(result,n**4)
    return pa.EncryptedNumber(pk,int(result))

# the SDC is designed to return INT!!! result
def SDC(sk: SK,e0: EN,e1: EN):
    x_,d_,E_r2=SDC_PART1(e0,e1)
    E_h=SDC_PART2(sk,x_,d_)
    e2=SDC_PART3(E_h,E_r2)
    return e2

if __name__=='__main__':
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    e0=pk.encrypt(101000)
    e1=pk.encrypt(101)
    e2=SDC(sk,e0,e1)
    # E_ara,E_brb,ra,rb=SM_PART1(e0,e1)
    # d0=sk.decrypt(E_ara)
    # d1=sk.decrypt(E_brb)
    # h=SM_PART2(sk,E_ara,E_brb)
    # e2=SM_PART3(h,e0,e1,ra,rb)
    d2=sk.decrypt(e2)
    print(d2)
