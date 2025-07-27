#encoding=utf-8
#Secure Multiplication
import random
import gmpy2
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK

# running on P1 with pk
def SM_PART1(E_a: EN,E_b: EN):
    pk=E_a.public_key
    if pk.n>1024:
        rand_max=1024
    else:
        rand_max=pk.n
    ra=random.randint(0,rand_max)
    rb=random.randint(0,rand_max)
    E_ra=pk.encrypt(ra)
    E_rb=pk.encrypt(rb)
    E_ara=E_a+E_ra
    E_brb=E_b+E_rb
    return E_ara,E_brb,ra,rb

# running on P2 with pk and sk
def SM_PART2(sk: SK,E_ara: EN,E_brb: EN)->EN:
    ha=sk.decrypt(E_ara)
    hb=sk.decrypt(E_brb)
    pk=E_ara.public_key
    n=pk.n
    #h=(ha*hb)%n
    h=ha*hb
    return pk.encrypt(h)

# running on P1 with pk
def SM_PART3(h: EN,E_a: EN,E_b: EN,ra: int,rb: int)->EN:
    pk=E_a.public_key
    n=pk.n
    temp1=gmpy2.mpz(E_a.ciphertext())
    s=gmpy2.mpz(h.ciphertext())*gmpy2.powmod(temp1,n-rb,n**2)
    temp2=gmpy2.mpz(E_b.ciphertext())
    s_=s*gmpy2.powmod(temp2,n-ra,n**2)
    E_rarb=pk.encrypt(ra*rb)
    result=s_*gmpy2.powmod(gmpy2.mpz(E_rarb.ciphertext()),n-1,n**2)
    E_ab=pa.EncryptedNumber(pk,int(result))
    return E_ab

def SM(sk: SK,e0: EN,e1: EN):
    E_ara,E_brb,ra,rb=SM_PART1(e0,e1)
    h=SM_PART2(sk,E_ara,E_brb)
    e2=SM_PART3(h,e0,e1,ra,rb)
    return e2

if __name__=='__main__':
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    e0=pk.encrypt(101)
    e1=pk.encrypt(101)
    e2=SM(sk,e0,e1)
    # E_ara,E_brb,ra,rb=SM_PART1(e0,e1)
    # d0=sk.decrypt(E_ara)
    # d1=sk.decrypt(E_brb)
    # h=SM_PART2(sk,E_ara,E_brb)
    # e2=SM_PART3(h,e0,e1,ra,rb)
    d2=sk.decrypt(e2)
    print(d2)
