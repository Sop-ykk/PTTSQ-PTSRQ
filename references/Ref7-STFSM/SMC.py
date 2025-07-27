#encoding=utf-8
#Secure Minimum Computation
import random
import gmpy2
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

# running on P1 with pk
def SMC_PART1(E_d0: EN,E_d1: EN,E_id0: EN,E_id1: EN):
    pk=E_d0.public_key
    n=pk.n
    if pk.n>1024:
        rand_max=1024
    else:
        rand_max=pk.n
    r1=random.randint(0,rand_max)
    r2=random.randint(0,rand_max)
    r3=random.randint(0,rand_max)
    r4=random.randint(0,rand_max)
    r5=random.randint(0,rand_max)
    E_r2=pk.encrypt(r2)
    E_r3=pk.encrypt(r3)
    E_r4=pk.encrypt(r4)
    E_r5=pk.encrypt(r5)
    F=random.randint(0,1)
    if F==0:
        E_a=(E_d0-E_d1)*r1
        E_b=E_d0+E_r2
        E_b_=E_d1+E_r3
        E_p=E_id0+E_r4
        E_p_=E_id1+E_r5
    else:
        E_a=(E_d1-E_d0)*r1
        E_b=E_d1+E_r2
        E_b_=E_d0+E_r3
        E_p=E_id1+E_r4
        E_p_=E_id0+E_r5
    return E_a,E_b,E_b_,E_p,E_p_,r2,r3,r4,r5

# running on P2 with pk and sk
def SMC_PART2(sk: SK,E_a: EN,E_b: EN,E_b_: EN,E_p: EN,E_p_: EN):
    a=sk.decrypt(E_a)
    pk=E_a.public_key
    n=pk.n
    if a>0:
        e=0
    else:
        e=1
    o=((E_b*(1-e))+(E_b_*e))
    o=gmpy2.powmod(o.ciphertext(),n-1,n**2)
    o=EN(pk,o)
    o_=((E_p*(1-e))+(E_p_*e))
    o_=gmpy2.powmod(o_.ciphertext(),n-1,n**2)
    o_=EN(pk,o_)
    E_e=pk.encrypt(e)
    return E_e,o,o_

# running on P1 with pk
def SMC_PART3(E_d0: EN,E_d1: EN,E_id0: EN,E_id1: EN,o: EN,o_: EN,E_e: EN,r2: int,r3: int,r4: int,r5: int)->EN:
    pk=E_d0.public_key
    n=pk.n
    E_1=pk.encrypt(1)
    E_dmin=((E_1-E_e)*r2)+(E_e*r3)+o+E_d1+E_d0
    E_idmin=((E_1-E_e)*r4)+(E_e*r5)+o_+E_id1+E_id0
    return E_dmin,E_idmin

def SMC(sk: SK,E_d0: EN,E_d1: EN,E_id0: EN,E_id1: EN):
    E_a,E_b,E_b_,E_p,E_p_,r2,r3,r4,r5=SMC_PART1(E_d0,E_d1,E_id0,E_id1)
    E_e,o,o_=SMC_PART2(sk,E_a,E_b,E_b_,E_p,E_p_)
    E_dmin,E_idmin=SMC_PART3(E_d0,E_d1,E_id0,E_id1,o,o_,E_e,r2,r3,r4,r5)
    return E_dmin,E_idmin

def SMC_(sk: SK,E_0: EN,E_1: EN)->EN:
    pk=E_0.public_key
    E_id0=pk.encrypt(0)
    E_id1=pk.encrypt(1)
    E_d,E_id=SMC(sk,E_0,E_1,E_id0,E_id1)
    return E_d

if __name__=='__main__':
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    e0=pk.encrypt(101000)
    e1=pk.encrypt(1010000)
    e2=pk.encrypt(1)
    e3=pk.encrypt(2)
    E_dmin,E_idmin=SMC(sk,e0,e1,e2,e3)
    # E_ara,E_brb,ra,rb=SM_PART1(e0,e1)
    # d0=sk.decrypt(E_ara)
    # d1=sk.decrypt(E_brb)
    # h=SM_PART2(sk,E_ara,E_brb)
    # e2=SM_PART3(h,e0,e1,ra,rb)
    d0=sk.decrypt(E_dmin)
    d1=sk.decrypt(E_idmin)
    print(d0)
    print(d1)

