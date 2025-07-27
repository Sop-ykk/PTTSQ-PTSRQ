from phe import paillier as pa
pk,sk=pa.generate_paillier_keypair(n_length=512)
e0=pk.encrypt(10000)
e1=pk.encrypt(2000)
e2=pk.encrypt(12000)
e3=e1+e0
e4=e0*10
e5=e0.public_key.encrypt(3000)
d0=sk.decrypt(e3)
d1=sk.decrypt(e4)
d2=sk.decrypt(e5)
print('a')