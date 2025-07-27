# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 21:06:30 2022
update version by ykl
@author: YiKelai
"""

#%%
## 1) Key Generation

from umbral import pre, keys, signing
import time

for traj_num in [50,100,1200,1600,2000]:
# traj_num = [500,1000,1500,2000]
    for ke in [5,10,15,20]:

        start = time.time()

        # Generate Umbral keys for Alice.
        alices_private_key = keys.SecretKey.random()
        alices_public_key = alices_private_key.public_key()

        alices_signing_key = keys.SecretKey.random()
        alices_verifying_key = alices_signing_key.public_key()
        alices_signer = signing.Signer(secret_key=alices_signing_key)

        # Generate Umbral keys for Bob.
        bobs_private_key = keys.SecretKey.random()
        bobs_public_key = bobs_private_key.public_key()

        stop1 = time.time()

        #%%
        ## 2) Enc \times traj_num

        for num1 in range(traj_num):
            # Encrypt data with Alice's public key.
            plaintext = b'Kelai YI, 510108190000000000, 18000000000, xxxxxxx'
            CCapsule, ciphertext = pre.encrypt(alices_public_key, plaintext)

        # Decrypt data with Alice's private key.
        # cleartext = pre.decrypt_original(ciphertext=ciphertext,
        #                         capsule=CCapsule,
        #                         delegating_sk=alices_private_key)

        stop2 = time.time()
            
        #%%
        ## 3) ReKey Generation

        # Alice generates "M of N" re-encryption key fragments (or "KFrags") for Bob.
        # In this example, 10 out of 20.
        kfrags = pre.generate_kfrags(delegating_sk=alices_private_key,
                                    signer=alices_signer,
                                    receiving_pk=bobs_public_key,
                                    threshold=10,
                                    shares=20)

        stop3 = time.time()

        # Several Ursulas perform re-encryption, and Bob collects the resulting `cfrags`.
        # He must gather at least `threshold` `cfrags` in order to activate the capsule.

        # Capsule = CapsuleFrag.verify(Capsule,capsule=Capsule,
        #                verifying_pk=alices_verifying_key,
        #                delegating_pk=alices_public_key,
        #                receiving_pk=bobs_public_key)

        #%%
        ## 4) ReEnc \times ke

        for num2 in range(ke):
            cfrags = list()           # Bob's cfrag collection
            for kfrag in kfrags[:10]:
                cfrag = pre.reencrypt(kfrag=kfrag, capsule=CCapsule)
                cfrags.append(cfrag)    # Bob collects a cfrag
                
        stop4 = time.time()
            
        #%% 
        ## 5) ReDec \times ke

        # Bob activates and opens the capsule
        # for cfrag in cfrags:
        #     capsule.attach_cfrag(cfrag)

        for num2 in range(ke):
                bob_cleartext = pre.decrypt_reencrypted(ciphertext=ciphertext,
                                            capsule=CCapsule,
                                            receiving_sk=bobs_private_key,
                                            verified_cfrags=cfrags,
                                            delegating_pk=alices_public_key)
                
                assert bob_cleartext == plaintext

        end = time.time()

        t1 = stop1 - start
        t2 = stop2 - stop1
        t3 = stop3 - stop2
        t4 = stop4 - stop3
        t5 = end - stop4
        duration = end-start

        print('n,'+str(traj_num)+',ke,'+str(ke)+',t1,'+str(t1)+',t2,'+str(t2)+',t3,'+str(t3)+',t4,'+str(t4)+',t5,'+str(t5)+',all,'+str(duration))
        # print(str(plaintext) +'\n'+str(cleartext)+'\n'+str(bob_cleartext))