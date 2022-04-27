from qcomp.gate import *
from qcomp.qregister import *
from qcomp.qbit import *
from qcomp.utils import *

def test_complex():
    c = ComplexMatrixGate([Hadamard, ID, Hadamard])
    assert c.qreg_size == 3

def test_matrix():
    m = MatrixGate([[1,2],[-1,0]], 1)
    assert m.qreg_size == 1

def test_smatrix():

    # construct matrices
    sN = sparseMatrix(3,3)
    sM = sparseMatrix(3,4)


    for i in range(3):
        for j in range(3):
            if i%2 == 0:
                sN.elements.append(element(i,j,1))
            if j%2 == 0:
                sM.elements.append(element(i,j,1))

    N = sN.convertToMatrix()
    M = sM.convertToMatrix()
    assert (np.matmul(N,M) == sN.multiplySparseMatrix(sM).convertToMatrix()).all()
