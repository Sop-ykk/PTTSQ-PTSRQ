#### 注意tenseal库的版本==0.3.12，用0.3.14不行
import sys
import numpy as np
import os

from CKKS_HE import *

from s0_data_pre import read_gen, data_pre, gen_ID, rec_ID
from s1_table import gen_allp, load_arHrs_, take_arHrs
from s2_CKKSvec2ten import save_context, load_context, load_EPloc_sr
from f2_tr2H import traj_enc, traj_dec
from s4_filter import nDSED, idRh
# from dH_pivot import get_arHrs
from s6_lookup import *
from s7_fhe import Dec_tr1, KE_SD, rv_esd, rk_esd, HE_vecmul3
import csv



if __name__ == '__main__':
    ordN = 11
    ro = 3
    so = 1
    dro = 3
    dso = 3
    all = True
    un, tl = 100, 60
    tlq = 5
    # tmax = 600
    # ke, kh = 5, 20
    per1 = 0.005     # 不同数据集的初始eps0依据百分位来设定，确保结果有效
    per2 = 0.05      # 过滤范围K，设定依据：扩大后的百分位 / 扩大后的eps

    radius = 5
    bi = True
    dist_type = 'sEu'
    file = 'GT11_tall600'
    files = ['GT11_tall600']
    vmax = 4096


    # -------------------------
    #### Step 0: DO 生成映射表
    # -------------------------
    # ordN1 = 2
    # ordN2 = 3
    # so1 = 0
    # so2 = 1
    # td, te = save_all(ordN1, ordN2, so1, so2)
    # -------------------------
    #### Step 0: 咱们读取映射表
    # -------------------------
    ctx = load_context()
    arHrs_, Ploc = load_arHrs_(ordN, dro, dso)
    arHrs = take_arHrs(arHrs_, ro, so)
    Ev_sr = load_EPloc_sr(ordN, ctx)

    for newfile in files:
    # newfile = files[0]
        #### Step 1: 数据准备
        data, trQ = read_gen(newfile, tlq, ordN, all, un, tl)
        pn_do_ = data[0:, 0:, 0:2]
        t1_do = data[0:, 0:, 2][0]      # 假设DO时间戳都齐全且一样
        p1_qu = trQ[0:, 0:2]
        t1_qu = trQ[0:, 2]
        un, tl = np.shape(pn_do_)[0], np.shape(pn_do_)[1]
        gID = gen_ID(un)

        #### step 2：DO 数据编码和ID加密
        time2s = time.time()
        hn_do_ = traj_enc(pn_do_, ordN).astype(int)
        Eidn_do = Enc_nvec(ctx, gID)
        time2e = time.time()

        #### step 3：QU 查询编码
        h1_qu = traj_enc(p1_qu, ordN).astype(int)
        time3e = time.time()

        #### step 4: CS 查表
        tid = np.nonzero(np.in1d(t1_do, t1_qu))[0]      # 匹配的tlq个时间戳，[tlq,]
        hn_do = hn_do_[:, tid]      # 取出所有轨迹在所查询的tlq个时间戳下的轨迹点h，[un, tlq]
        time4e = time.time()

        # -------------------------
        ##### (eps值的预估，是场景模拟，不计入运行时间)
        pn_do = pn_do_[:, tid]
        dne = nDSED(dist_type, p1_qu, pn_do, t1_qu)
        eps0 = np.quantile(dne, per1)   # 条件预假设, 用于后续精确验证
        # -------------------------

        #### step 5: CS 初筛per2
        time5s = time.time()
        ir_Pv = idRh(per2, h1_qu, hn_do, t1_qu, arHrs, radius, bi)
        time5e = time.time()

        # -------------------------
        # 从top-kh之中，精确搜索top-ke
        # -------------------------

        #### step 6: CS 获取密文轨迹坐标
        ep1_qu = look_up_tr1(h1_qu, Ev_sr, ctx, vmax)      # ls长度为traj_len
        epK_do = look_up_trK(hn_do, ir_Pv, Ev_sr, ctx, vmax)  # ls长度为kh, 每个元素lsi也是列表，长度traj_len
        # time6e = time.time()

        #### step 7: CS 密文计算欧式距离(通分)
        edK = KE_SD(ep1_qu, epK_do, t1_qu) 
        # Dec_tr1(edK)

        #### step 8: CS 密文范围验证        
        h = len(t1_qu) - 1
        vecf = rv_esd(edK, eps0 * (2 * (t1_qu[h] - t1_qu[0])))
        # Dec_vec(vecf)
        # lsf = rk_esd(edK, eps0)
        # Dec_tr1(lsf)
        time8e = time.time()

        #### step 9: CS 同态乘法
        vecf_ = Bootstrap_vec(vecf)
        # lsf = Esubs_vec2ls(vecf_, len(edK))
        hK_do = hn_do_[ir_Pv]
        ## 假设Eidn_do是轨迹ID_num密文向量的列表, 有len(ID_num)个元素, 每个元素是un维向量, 则如下 2023.12.22以后
        # Eidn_do = Enc_nvec(ctx, gID)
        EidK_do = [Eidn_do[j] for j in ir_Pv]
        rs_trH, rs_Eid = HE_vecmul3(vecf_, hK_do, EidK_do, len(edK), tl)

        time9e = time.time()

        #### step 10: CS 代理重加密
        # 用openfhe的库

        #### step 11: QU 解密
        rs_tr = Dec_vls_half(rs_trH)        # ckks解密, 恢复Hilbert编码值
        rs_trP = traj_dec(rs_tr, ordN)      # Hilbert解码, 恢复轨迹坐标明文
        rs_id = Dec_vls_half(rs_Eid)        # ckks解密, 恢复ASCII码值
        rs_idP = rec_ID(rs_id)              # ASCII解码, 恢复轨迹坐标ID明文
        time11e = time.time()

        # ------------ Timing -------------
        t_DOencode = time2e - time2s
        t_QUrequest = time3e - time2e
        t_CSfilter = time4e - time3e + time5e - time5s
        t_CSverify = time8e - time5e
        t_CSmultiply = time9e - time8e
        t_QUrecover = time11e - time9e


        print('DO_encode', t_DOencode, '\nQU_request', t_QUrequest, '\nCS_filter', t_CSfilter, '\nCS_verify', t_CSverify, '\nCS_multiply', t_CSmultiply, '\nQU_recover', t_QUrecover)

        with open('./Figures/time_res.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([str(newfile),'per1',str(per1),'esp',str(eps0),'per2',str(per2),'tlq',str(tlq),'ordN',str(ordN),'rs',str(ro),str(so),'r',str(radius),'bi',str(bi)[0],dist_type,'all',str(all)[0],np.shape(data)[0], np.shape(data)[1]])  # 写入表头
            writer.writerow(['TIME', 'DO_encode', t_DOencode, 'QU_request', t_QUrequest, 'CS_filter', t_CSfilter, 'CS_verify', t_CSverify, 'CS_multiply', t_CSmultiply, 'QU_recover', t_QUrecover])