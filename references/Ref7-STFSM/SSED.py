#encoding=utf-8
#Basic Secure Squared Euclidean Distance Protocol
#Although BSSED support multiple points
#But I think 2 is enough
import random
import gmpy2
from typing import List
from phe import paillier as pa
from phe.paillier import EncryptedNumber as EN
from phe.paillier import PaillierPublicKey as PK
from phe.paillier import PaillierPrivateKey as SK
from SM import *

def SSED(sk: SK,E_xx: EN,E_xy: EN,E_yx: EN,E_yy: EN)-> EN:
    E_xxyx=E_xx-E_yx
    E_xyyy=E_xy-E_yy
    E_1=SM(sk,E_xxyx,E_xxyx)
    E_2=SM(sk,E_xyyy,E_xyyy)
    E=E_1+E_2
    return E

if __name__=="__main__":
    pk,sk=pa.generate_paillier_keypair(n_length=512)
    e0x=pk.encrypt(0)
    e0y=pk.encrypt(6)
    e1x=pk.encrypt(4)
    e1y=pk.encrypt(0)
    res=SSED(sk,e0x,e0y,e1x,e1y)
    print(sk.decrypt(res))