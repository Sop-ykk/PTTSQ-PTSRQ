import tenseal as ts
import numpy as np
from scipy.special import comb

# -------------------------
### CKKS同态加密
# -------------------------
def gencontext():
    # context = ts.context(ts.SCHEME_TYPE.CKKS, 8192, coeff_mod_bit_sizes=[22 ,21, 21, 21, 21, 21, 21, 21, 21, 21])
    # context.global_scale = pow(2, 21)
    context = ts.context(ts.SCHEME_TYPE.CKKS, 8192, coeff_mod_bit_sizes=[60, 40, 40, 60])


    context.global_scale = pow(2, 40)

    context.generate_galois_keys()
    return context

def Enc_vec(context, ls_vector):
    ## ls_vector 是 列表
    return ts.ckks_vector(context, ls_vector)

def Dec_vec(enc_vector):
    return enc_vector.decrypt()

def Enc_ten(context, np_tensor):
    ## np_tensor 是 数组
    return ts.ckks_tensor(context, np_tensor)

def Dec_ten(enc_tensor):
    return np.array(enc_tensor.decrypt().tolist())

def Dec_vls_half(evls):
    # evls是一个列表，其中每个元素为一个ckks.vector
    # 返回的明文值/2, numpy_array取整
    vls = []
    for i in range(len(evls)):
        vls.append(Dec_vec(evls[i]))
    intarr_half = (np.array(vls)/2).astype(int)
    return intarr_half

def Bootstrap_vec(enc_vector):
    ctx = enc_vector.context()
    new_encvec = Enc_vec(ctx, Dec_vec(enc_vector))
    return new_encvec

def Enc_nvec(context, np_arr):
    ## np_tensor 是 数组, 按行向量加密为列数长的列表, 解密用Dec_vec(a[i])
    els = []
    for i in range(len(np_arr)):
        els.append(Enc_vec(context, np_arr[i]))
    return els

# -------------------------
### CKKS同态比较，符号函数
# -------------------------
def f_n(x, n: int=1):
    s = 0
    # if n > 0 and -1 <= x <= 1:
    for i in range(n+1):
        s += (comb(2*i, i) * x * (1-x**2)**i)/(4**i)
    return s

def f_n_d(x, n: int=1, d: int=8):
    m = f_n(x, n)
    for i in range(d):
        m = f_n(m, n)
    return m

def g_n(x, n: int=1):
    if n == 1:
        g = (-1359 * x**3 + 2126 * x)/(2**10)
    if n == 2:
        g = (3796 * x ** 5 - 6108 * x ** 3 + 3334 * x) / (2 ** 10)
    if n == 3:
        g = (-12860 * x ** 7 + 25614 * x ** 5 - 16577 * x ** 3 + 4589 * x) / (2 ** 10)
    return g

def fg(x, nf: int=1, ng: int=1, df: int=8, dg: int=6):
    m = g_n(x, ng)
    for i in range(dg):
        m = g_n(m, ng)
    for j in range(df):
        m = f_n(m, nf)
    return m

# def sigmoid(x):
#     return 1/(1+np.exp(-x))


# -------------------------
### CKKS密文打包 ls->vec
# -------------------------
def Epack_ls2vec(edK):
    '''
    对于edK密文列表首先进行打包, 得到一个向量evK, 再与eps批量化符号比较
    
    '''
    I = np.eye(len(edK))
    evK = edK[0] @ [I[0]]       # @的作用等价于enc_vec.matmul(plain_matrix)
    for i in range(1, len(edK)):
        evK = evK + edK[i] @ [I[i]]
    return evK

def Esubs_vec2ls(vecf, vm):
    '''
    密文向量vecf打包后不可直接拆分, 影响矩阵运算, 此函数用于将vecf拆分为单个密文的列表lsf
    vm: 向量长度, 即, len(Dec_vec(vecf))
    '''
    lsf = []
    I = np.eye(vm)
    for i in range(vm):
        lsf.append(vecf @ (I[i].reshape((-1,1))))
    return lsf



if __name__ == "__main__":
    context = gencontext()

### tensor
    a = np.array([[1.,2.,3.,4.], [1.,2.,5.,4.]])
    b = np.array([[1.,3.,5.,7.], [2.,4.,6.,8.]])
    enc_a = Enc_ten(context, a)
    enc_at = Enc_ten(context, a.T)
    enc_b = Enc_ten(context, b)
    enc_bt = Enc_ten(context, b.T)
    # res = enc_a + enc_b
    # res = enc_a - enc_b
    res = enc_a * enc_b      # 矩阵乘法（点乘）
    # res = enc_a @ enc_bt   # 矩阵乘法（叉乘）
    print(Dec_ten(res))