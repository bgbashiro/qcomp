from qcomp.utils import *
import numpy as np
from scipy.sparse import random

TOL = 1e-20

def test_scalar_sparse():
    for _ in range(10):
        n = np.random.randint(10,20)
        m = np.random.randint(10,20)
        A = random(m,n, density=0.3).toarray()
        sA = arr_to_sparse(A)
        res = A * 5
        sres = sA.mult_scalar(5)
        assert ((sres.to_np_array() - res)**2).sum() < TOL

def test_mult_sparse():
    for _ in range(10):
        n = np.random.randint(10,30)
        m = np.random.randint(10,30)
        k = np.random.randint(10,30)
        A = random(n,k,density=0.3).toarray()
        B = random(k,m,density=0.3).toarray()
        sA = arr_to_sparse(A)
        sB = arr_to_sparse(B)
        res = np.matmul(A,B)
        sres = sA.mult_sparse(sB)
        assert ((sres.to_np_array() - res)**2).sum() < TOL